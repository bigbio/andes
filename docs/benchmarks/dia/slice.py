import sys,re
src,dst,lo,hi=sys.argv[1],sys.argv[2],float(sys.argv[3]),float(sys.argv[4])
rt_re=re.compile(rb'<cvParam[^>]*name="scan start time"[^>]*>')
val_re=re.compile(rb'\svalue="([\d.eE+-]+)"')


def scan_seconds(blk):
    m=rt_re.search(blk)
    v=m and val_re.search(m.group(0))
    if not v:
        return -1
    return float(v.group(1))*(60 if b'unitName="minute"' in m.group(0) else 1)

with open(src,'rb') as f, open(dst,'wb') as o:
    head=True; buf=[]; n=0
    for line in f:
        s=line.lstrip()
        if head:
            if s.startswith(b'<spectrum '):
                head=False; buf=[line]; continue
            o.write(line.replace(b'<spectrumList count="', b'<spectrumList countX="')); continue
        if buf:
            buf.append(line)
            if s.startswith(b'</spectrum>'):
                blk=b''.join(buf); buf=[]
                rt=scan_seconds(blk)
                if lo<=rt<=hi: o.write(blk); n+=1
                elif rt>hi: break
            continue
        if s.startswith(b'<spectrum '): buf=[line]; continue
        if s.startswith(b'</spectrumList>'): break
    o.write(b'\t\t</spectrumList>\n\t</run>\n</mzML>\n')
print(n)
