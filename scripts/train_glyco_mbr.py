#!/usr/bin/env python3
"""Train andes' intact-glycopeptide MBR match scorer.

The scorer is a 10 -> 64 -> 32 -> 8 -> 2 MLP with BatchNorm1d and ReLU on
every layer, trained on labels andes generates itself with
`--glyco-mbr-match-train`. It reimplements the network and training recipe of
the published method (Zeng et al., Nat. Commun. 2022; focal loss gamma=2,
Adam lr 1e-3 weight decay 1e-3, StepLR(step=2, gamma=0.99) stepped per batch,
100 epochs, batch 512) because the authors' pretrained weights are not
distributed under a redistributable license.

Input: the TSV written by `--glyco-mbr-match-train`:
    label  f0 .. f9  peptide  donor>target
Labels are 1 for a scored window that covers the peak the target run itself
identified, 0 otherwise (including every window of the +11 Th decoy twin).

The peptide column is the split group: a peptide's rows go entirely to train or
entirely to validation, so the reported AUC is not inflated by the same peptide
appearing on both sides.

Output: the text model read by `glyco_match::Model::parse`:
    arch 10 64 32 8 2
    fc1.weight <640 values> ...
"""

import argparse
import hashlib
import os
import sys

import torch
import torch.nn as nn
import torch.nn.functional as F

N_FEATURES = 10


class MyNetNN(nn.Module):
    """The published forward pass, unchanged."""

    def __init__(self):
        super().__init__()
        self.fc1 = nn.Linear(10, 64)
        self.drop1 = nn.Dropout(0.5)
        self.bn1 = nn.BatchNorm1d(64)
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 8)
        self.fc4 = nn.Linear(8, 2)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.drop1(x)
        x = self.bn1(x)
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        return F.relu(self.fc4(x))


def focal_loss(preds, labels, gamma=2.0):
    log_soft = F.log_softmax(preds, dim=1)
    soft = torch.exp(log_soft).gather(1, labels.view(-1, 1))
    log_soft = log_soft.gather(1, labels.view(-1, 1))
    return -torch.mul(torch.pow(1.0 - soft, gamma), log_soft).mean()


def load_rows(path, max_negatives):
    """Read the dump, keeping every positive and up to `max_negatives`
    negatives per (peptide, run pair) candidate set."""
    by_set = {}
    labels = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != N_FEATURES + 3:
                raise SystemExit(f"{path}: bad row with {len(cols)} columns: {line[:120]!r}")
            label = int(cols[0])
            feats = [float(v) for v in cols[1 : 1 + N_FEATURES]]
            group, pair = cols[1 + N_FEATURES], cols[2 + N_FEATURES]
            key = (group, pair)
            by_set.setdefault(key, []).append((label, feats))
            labels[key] = group

    rows, groups = [], []
    kept = {0: 0, 1: 0}
    for key, items in by_set.items():
        pos = [x for x in items if x[0] == 1]
        neg = [x for x in items if x[0] == 0]
        # Deterministic subsample of the negatives of this candidate set.
        neg.sort(key=lambda r: hashlib.blake2b(
            repr(r[1]).encode(), digest_size=8).digest())
        for label, feats in pos + neg[:max_negatives]:
            rows.append((label, feats))
            groups.append(labels[key])
            kept[label] += 1
    return rows, groups, kept[1], kept[0]


def split_by_group(groups, val_fraction):
    """Deterministic group split: a peptide never straddles the two sides."""
    train, val = [], []
    for i, g in enumerate(groups):
        h = hashlib.blake2b(g.encode(), digest_size=8).digest()
        frac = int.from_bytes(h[:8], "little") / 2**64
        (val if frac < val_fraction else train).append(i)
    return train, val


def roc_auc(scores, labels):
    pairs = sorted(zip(scores, labels))
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    rank_sum = 0.0
    i = 0
    while i < len(pairs):
        j = i
        while j + 1 < len(pairs) and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        rank_sum += avg_rank * sum(lab for _, lab in pairs[i : j + 1])
        i = j + 1
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def write_model(model, path, provenance):
    sd = {k: v.detach().to(torch.float64).contiguous() for k, v in model.state_dict().items()}
    with open(path, "w", encoding="utf-8") as fh:
        for line in provenance:
            fh.write(f"# {line}\n")
        fh.write("arch 10 64 32 8 2\n")
        for key in (
            "fc1.weight", "fc1.bias",
            "bn1.weight", "bn1.bias", "bn1.running_mean", "bn1.running_var",
            "fc2.weight", "fc2.bias",
            "fc3.weight", "fc3.bias",
            "fc4.weight", "fc4.bias",
        ):
            vals = sd[key].reshape(-1).tolist()
            fh.write(key + " " + " ".join(f"{v:.17g}" for v in vals) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", help="TSV from --glyco-mbr-match-train")
    ap.add_argument("--out", required=True, help="model text file to write")
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--lr", type=float, default=0.001)
    ap.add_argument("--weight-decay", type=float, default=0.001)
    ap.add_argument("--max-negatives", type=int, default=5,
                    help="negatives kept per (peptide, run pair) candidate set")
    ap.add_argument("--val-fraction", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--source", default="", help="provenance string recorded in the header")
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    rows, groups, n_pos, n_neg = load_rows(args.table, args.max_negatives)
    print(f"rows: {len(rows)} ({n_pos} positive, {n_neg} negative) from {args.table}",
          file=sys.stderr)
    if n_pos == 0 or n_neg == 0:
        raise SystemExit("training table needs both classes")

    train_idx, val_idx = split_by_group(groups, args.val_fraction)
    if not val_idx:
        raise SystemExit("no validation groups; lower --val-fraction or add data")

    def tensor(idx):
        x = torch.tensor([rows[i][1] for i in idx], dtype=torch.float32)
        y = torch.tensor([rows[i][0] for i in idx], dtype=torch.long)
        return x, y

    x_tr, y_tr = tensor(train_idx)
    x_va, y_va = tensor(val_idx)
    print(f"train: {len(train_idx)} rows; validation: {len(val_idx)} rows "
          f"({int(y_va.sum())} positive)", file=sys.stderr)

    model = MyNetNN()
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    sched = torch.optim.lr_scheduler.StepLR(opt, step_size=2, gamma=0.99)
    gen = torch.Generator().manual_seed(args.seed)
    n_batches = (len(train_idx) + args.batch_size - 1) // args.batch_size

    best_auc, best_state = -1.0, None
    for epoch in range(args.epochs):
        model.train()
        order = torch.randperm(len(train_idx), generator=gen)
        total, correct, loss_sum = 0, 0, 0.0
        for b in range(n_batches):
            idx = order[b * args.batch_size : (b + 1) * args.batch_size]
            if len(idx) < 2:
                # BatchNorm1d cannot compute a batch statistic from one row.
                continue
            xb, yb = x_tr[idx], y_tr[idx]
            opt.zero_grad()
            pred = model(xb)
            loss = focal_loss(pred, yb)
            loss.backward()
            opt.step()
            sched.step()
            total += len(idx)
            correct += int((pred.argmax(1) == yb).sum())
            loss_sum += float(loss) * len(idx)
        model.eval()
        with torch.no_grad():
            logits = model(x_va)
            prob = torch.softmax(logits, dim=1)[:, 1].tolist()
        auc = roc_auc(prob, y_va.tolist())
        if auc == auc and auc > best_auc:
            best_auc = auc
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        if epoch % 10 == 0 or epoch == args.epochs - 1:
            print(f"epoch {epoch:3d}: loss {loss_sum / total:.4f}, "
                  f"train acc {correct / total:.4f}, val AUC {auc:.4f}", file=sys.stderr)

    if best_state is not None:
        model.load_state_dict(best_state)
    provenance = [
        "andes intact-glycopeptide MBR match scorer (10 -> 64 -> 32 -> 8 -> 2 MLP,",
        "BatchNorm1d + ReLU on every layer; method of Zeng et al., Nat. Commun. 2022).",
        "Trained by scripts/train_glyco_mbr.py with the published recipe:",
        "focal loss gamma=2, Adam lr 1e-3 weight decay 1e-3,",
        "StepLR(step_size=2, gamma=0.99) per batch, 100 epochs, batch 512.",
        f"Training table: a --glyco-mbr-match-train dump "
        f"({os.path.getsize(args.table)} bytes)",
        f"Rows: {len(rows)} ({n_pos} positive, {n_neg} negative); "
        f"groups split {1 - args.val_fraction:.0%}/{args.val_fraction:.0%} by peptide.",
        f"Validation ROC AUC: {best_auc:.4f} (peptide-disjoint).",
    ]
    if args.source:
        provenance.append(args.source)
    write_model(model, args.out, provenance)
    print(f"wrote {args.out} (validation AUC {best_auc:.4f})", file=sys.stderr)


if __name__ == "__main__":
    main()
