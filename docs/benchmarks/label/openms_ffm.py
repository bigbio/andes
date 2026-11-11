"""Run OpenMS FeatureFinderMultiplex (pyOpenMS) on one run and write its
multiplets as TSV: rt, charge, then per channel m/z and intensity."""
import sys, time
import pyopenms as p
mzml, labels, mc, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
tol = float(sys.argv[5]) if len(sys.argv) > 5 else 6.0
cutoff = float(sys.argv[6]) if len(sys.argv) > 6 else 1000.0
exp = p.MSExperiment()
fh = p.MzMLFile()
o = fh.getOptions(); o.setMSLevels([1]); fh.setOptions(o)
fh.load(mzml, exp)
t = time.time()
ff = p.FeatureFinderMultiplexAlgorithm()
prm = ff.getDefaults()
prm.setValue("algorithm:labels", labels)
prm.setValue("algorithm:max_nr_labelled_aas", mc)
prm.setValue("algorithm:charge", "2:5")
prm.setValue("algorithm:knock_out", "true")
prm.setValue("algorithm:mz_tolerance", tol)
prm.setValue("algorithm:intensity_cutoff", cutoff)
ff.setParameters(prm)
ff.run(exp, True)
cmap = ff.getConsensusMap()
fmap = ff.getFeatureMap()
print(f"FeatureFinderMultiplex: {cmap.size()} multiplets in {time.time()-t:.1f} s", file=sys.stderr)
with open(out, "w") as f:
    f.write("rt\tcharge\tchannels\n")
    for c in cmap:
        hs = sorted(c.getFeatureList(), key=lambda h: h.getMapIndex())
        f.write(f"{c.getRT():.2f}\t{c.getCharge()}\t" + ";".join(f"{h.getMapIndex()}:{h.getMZ():.5f}:{h.getIntensity():.1f}" for h in hs) + "\n")
