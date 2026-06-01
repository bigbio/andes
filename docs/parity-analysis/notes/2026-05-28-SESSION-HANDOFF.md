# Session handoff — 2026-05-28 (PSM-gain + chimeric)

Saved before a machine restart. Roadmap: `docs/2026-05-28-psm-gain-state-and-roadmap.md`.
Goal: beat Java on PSMs @1% FDR on **PXD001819 + TMT** (the merge gate; speed already wins). Nothing merged.

## State (all VM-bench-validated)

| Item | Result | Where |
|---|---|---|
| **DeltaRawScore (3a)** | ✅ clean win, kept: PXD 14,808→**14,937** (+129), TMT 9,605→**9,617** (+12), Astral 36,715→**36,819** (+104); zero wall cost. PXD now −0.25% vs Java. | `feat/delta-raw-score` `bea5d697` |
| Drop lnEValue (3b) | discarded — noise (−8/+9/+18), costs PXD. Java/Rust E-value formula already matches. | bench note below |
| TMT (−4.9%, gate blocker) | cheap levers **all ruled out** (mods, args, deconvolution=iter30, scoring param=parity). Gap is in the **GF SpecE DP** = Lever 2a (research). Not started. | bench note below |
| Lever 2b (top-n 2) | ❌ FDR inflation (Astral canary, decoy collapse), reverted. | bench note below |
| Chimeric fragment-overlap diagnostic | built (env-gated `MSGF_CHIMERIC_OVERLAP=1`); BSA preview = **low overlap** → tentatively challenges Phase-3 "fragment theft". Astral measurement PENDING. | `feat/chimeric-dda-plus` `59421180` |

Notes: `2026-05-28-delta-raw-score-bench.md`, `2026-05-28-chimeric-fragment-overlap-diagnostic.md`.

## Resume after restart

1. **Reopen the VM socket** (`/tmp/msgfplus-bench.sock` dies on restart): re-establish the ControlMaster to `pride-linux-vm`. Bench root `/srv/data/msgf-bench/`.
2. **Decisive chimeric measurement** (was blocked only by the socket): ship `feat/chimeric-dda-plus` (`59421180`) to the VM, rebuild, run Astral chimeric with `MSGF_CHIMERIC_OVERLAP=1 ... --chimeric 2>astral-overlap.log`, aggregate `CHIM_OVERLAP` (awk recipe in the diagnostic note).
   - High overlap → Phase 3 (shared-fragment competition) validated. Low (BSA pattern) → theft refuted → Phase 3 won't help; chimeric needs per-scan FDR or stays shelved.
3. **Strategic fork (your call):** TMT Lever-2a (GF SpecE-shape GF-trace + fix, ~1-2wk, Rule-2 risk) vs bank 3a (still gate-blocked since TMT unmoved) vs keep pushing chimeric.

Branches: `feat/delta-raw-score` (3a, base `a71553ce` = post-PR#40 mainline) · `feat/chimeric-dda-plus` (roadmap + chimeric, parked) · `dev` = Java (no Rust crates).
Full context also in auto-memory: `project_psm_gain_2026_05_28_session.md`.
</content>
