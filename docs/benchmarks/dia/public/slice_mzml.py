"""slice_mzml.py <in.mzML> <out.mzML> <rt_lo_s> <rt_hi_s>: keep the spectra whose scan start time
is in [rt_lo, rt_hi] seconds, unchanged, with the mzML header (instrument, data processing),
as a plain (non-indexed) mzML with a correct spectrumList count and spectrum index 0..n-1.
Prints the number of spectra, of MS1 and of MS2 kept. Same scan selection as slice.py.
"""
import re, sys, tempfile, os, shutil

src, dst, lo, hi = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
rt_re = re.compile(rb'<cvParam[^>]*name="scan start time"[^>]*>')
val_re = re.compile(rb'\svalue="([\d.eE+-]+)"')
lvl_re = re.compile(rb'value="(\d+)" name="ms level"|name="ms level" value="(\d+)"')
idx_re = re.compile(rb'(<spectrum [^>]*?\bindex=")(\d+)(")')
cnt_re = re.compile(rb'<spectrumList count="\d+"')


def scan_seconds(blk):
    m = rt_re.search(blk)
    v = m and val_re.search(m.group(0))
    if not v:
        return -1
    return float(v.group(1)) * (60 if b'unitName="minute"' in m.group(0) else 1)


head = []
n = {1: 0, 2: 0}
tmp = tempfile.NamedTemporaryFile(dir=os.path.dirname(os.path.abspath(dst)), delete=False)
with open(src, 'rb') as f:
    in_head = True
    buf = None
    for line in f:
        s = line.lstrip()
        if in_head:
            if s.startswith(b'<spectrum '):
                in_head = False
                buf = [line]
                continue
            if s.startswith(b'<indexedmzML'):
                continue
            head.append(line)
            continue
        if buf is not None:
            buf.append(line)
            if s.startswith(b'</spectrum>'):
                blk = b''.join(buf)
                buf = None
                rt = scan_seconds(blk)
                if lo <= rt <= hi:
                    lv = lvl_re.search(blk)
                    lv = int(lv.group(1) or lv.group(2)) if lv else 0
                    n[lv] = n.get(lv, 0) + 1
                    blk = idx_re.sub(lambda m: m.group(1) + str(sum(n.values()) - 1).encode() + m.group(3),
                                     blk, count=1)
                    tmp.write(blk)
                elif rt > hi:
                    break
            continue
        if s.startswith(b'<spectrum '):
            buf = [line]
            continue
        if s.startswith(b'</spectrumList>'):
            break
tmp.close()
total = sum(n.values())
with open(dst, 'wb') as o:
    for line in head:
        o.write(cnt_re.sub(b'<spectrumList count="%d"' % total, line))
    with open(tmp.name, 'rb') as t:
        shutil.copyfileobj(t, o, 1 << 24)
    o.write(b'\t\t</spectrumList>\n\t</run>\n</mzML>\n')
os.unlink(tmp.name)
print(f'spectra {total} MS1 {n.get(1, 0)} MS2 {n.get(2, 0)}')
