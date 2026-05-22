import pathlib, numpy as np, sys
# usage: python compare_tol.py <baseline_dir> <new_dir> [glob]
# Tolerance-based, zero-crossing-aware comparator for CLASS .dat outputs.
# A column entry passes if |a-b| <= atol + RTOL*|a|, with atol = RTOL*column_peak.
# Near a zero-crossing (|a|~0) the absolute floor (0.1% of the column's peak
# magnitude) dominates, so meaningless blown-up relative errors at TE/Ephi
# crossings are ignored; peak-region points are judged at ~0.1% relative.
if len(sys.argv) < 3:
    sys.exit("usage: python compare_tol.py <baseline_dir> <new_dir> [glob]")
base, new = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
glob = sys.argv[3] if len(sys.argv) > 3 else "*.dat"
RTOL = 1e-3
ok = True
for f in sorted(base.glob(glob)):
    n = new / f.name
    if not n.exists(): print(f"MISSING {f.name}"); ok = False; continue
    # ndmin=2 keeps single-row files 2D so colpeak stays per-column.
    a, b = np.loadtxt(f, comments="#", ndmin=2), np.loadtxt(n, comments="#", ndmin=2)
    if a.shape != b.shape: print(f"SHAPE {f.name}: {a.shape} vs {b.shape}"); ok = False; continue
    if a.size == 0: print(f"OK   {f.name}: empty"); continue
    scale = np.maximum(np.abs(a), np.abs(b))
    colpeak = scale.max(axis=0, keepdims=True)          # per-column peak (always 2D via ndmin=2)
    atol = RTOL * np.where(colpeak > 0, colpeak, 1.0)   # absolute floor = 0.1% of column peak
    diff = np.abs(a - b)
    bad = diff > (atol + RTOL * np.abs(a))
    nbad = int(bad.sum())
    worst = float((diff / np.where(colpeak > 0, colpeak, 1.0)).max())  # worst |a-b|/colpeak
    flag = "OK  " if nbad == 0 else "FAIL"
    if nbad: ok = False
    print(f"{flag} {f.name}: n_exceed={nbad} worst_vs_colpeak={worst:.2e}")
sys.exit(0 if ok else 1)
