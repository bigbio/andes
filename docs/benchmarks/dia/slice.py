import sys,re
src,dst,lo,hi=sys.argv[1],sys.argv[2],float(sys.argv[3]),float(sys.argv[4])
rt_re=re.compile(rb'name="scan start time" value="([\d.eE+-]+)"')
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
                m=rt_re.search(blk); rt=float(m.group(1)) if m else -1
                if lo<=rt<=hi: o.write(blk); n+=1
                elif rt>hi: break
            continue
        if s.startswith(b'<spectrum '): buf=[line]; continue
        if s.startswith(b'</spectrumList>'): break
    o.write(b'\t\t</spectrumList>\n\t</run>\n</mzML>\n')
print(n)
