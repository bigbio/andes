"""fetch.py DATA_DIR: download the raws, FASTAs and truth tables of datasets.json.

PRIDE serves one connection at about 1 MB/s; the raws are fetched in parallel 64 MB ranges.
"""
import concurrent.futures as cf, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PRIDE = 'https://ftp.ebi.ac.uk/pride/data/archive'
IMPX = 'https://raw.githubusercontent.com/vbc-proteomics-org/imp-x-fdr/master/IMP_X_FDR/Resources/libraries'
SMALL = {
    'db/re-analysis_PXD014337.fasta': f'{PRIDE}/2023/08/PXD041955/re-analysis_PXD014337.fasta',
    'db/peplib.fasta': f'{PRIDE}/2023/08/PXD041955/peplib.fasta',
    'db/all_320_proteins_from_first_second_third_fourth_fifth_batch.fasta':
        f'{PRIDE}/2024/08/PXD042173/all_320_proteins_from_first_second_third_fourth_fifth_batch.fasta',
    'db/Ribosome_EColi_SG_AND_human_abundant500excl-ribosomal_671sequences.fasta':
        f'{PRIDE}/2022/05/PXD029252/Ribosome_EColi_SG_AND_human_abundant500excl-ribosomal_671sequences.fasta',
    'truth/protein-construct-sequencesandinteraction-groups-corrected.xlsx':
        f'{PRIDE}/2024/08/PXD042173/protein-construct-sequencesandinteraction-groups-corrected.xlsx',
    'truth/main_library.xlsx': f'{IMPX}/main_library.xlsx',
    'truth/enrichable_library.xlsx': f'{IMPX}/enrichable_library.xlsx',
    'truth/peplib_Cas9_Beveridge.xlsx': f'{IMPX}/peplib_Cas9_Beveridge.xlsx',
}
CHUNK = 64 * 1024 * 1024


def size(url):
    out = subprocess.run(['curl', '-sSIL', url], capture_output=True, text=True).stdout
    return int([l for l in out.lower().splitlines() if l.startswith('content-length')][-1].split(':')[1])


def part(url, path, lo, hi):
    for _ in range(8):
        if os.path.exists(path) and os.path.getsize(path) == hi - lo + 1:
            return
        subprocess.run(['curl', '-sS', '-f', '--max-time', '600', '-r', f'{lo}-{hi}', '-o', path, url])
    raise RuntimeError(f'{path}: download failed')


def fetch_large(url, dest):
    if os.path.exists(dest):
        return
    n = size(url)
    parts = [(f'{dest}.part{i}', lo, min(lo + CHUNK, n) - 1) for i, lo in enumerate(range(0, n, CHUNK))]
    with cf.ThreadPoolExecutor(10) as ex:
        list(ex.map(lambda p: part(url, *p), parts))
    with open(dest + '.tmp', 'wb') as out:
        for p, _, _ in parts:
            with open(p, 'rb') as f:
                out.write(f.read())
            os.remove(p)
    assert os.path.getsize(dest + '.tmp') == n
    os.rename(dest + '.tmp', dest)
    print('fetched', dest, flush=True)


def main(data):
    for sub in ('raw', 'db', 'truth'):
        os.makedirs(f'{data}/{sub}', exist_ok=True)
    for rel, url in SMALL.items():
        if not os.path.exists(f'{data}/{rel}'):
            subprocess.run(['curl', '-sS', '-fL', '-o', f'{data}/{rel}', url], check=True)
    for d in json.load(open(f'{HERE}/datasets.json')):
        fetch_large(d['url'], f'{data}/raw/{d["raw"]}')


if __name__ == '__main__':
    main(sys.argv[1])
