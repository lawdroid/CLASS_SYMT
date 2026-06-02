import json, math, random
random.seed(20260602)

def randn():
    # Box-Muller
    u1 = random.random()
    u2 = random.random()
    return math.sqrt(-2*math.log(u1+1e-12)) * math.cos(2*math.pi*u2)

def gen_branch(branch, n=5000):
    samples = []
    for _ in range(n):
        if branch == "3a":
            beta = 0.0
            Om   = 0.310 + 0.008*randn()
            h    = 0.682 + 0.005*randn()
            s8   = 0.810 + 0.012*randn()
            wch  = -0.992 + 0.004*randn()
            cD   = 18 + 8*abs(randn())
        elif branch == "3b":
            # skewed unimodal in beta on (0, 2.4]
            while True:
                u1 = random.random()+1e-9
                u2 = random.random()+1e-9
                b = -0.45 * (math.log(u1) + math.log(u2))
                if 0 < b <= 2.4: break
            beta = b
            Om   = 0.310 + 0.008*randn() + 0.012*beta
            h    = 0.682 + 0.006*randn() - 0.006*beta
            s8   = 0.810 + 0.012*randn() - 0.018*beta
            wch  = -0.990 + 0.005*randn() - 0.002*beta
            cD   = 18 + 9*abs(randn()) + 1.5*beta
        else:  # 3c
            mode = 2.78 if random.random() < 0.62 else 3.22
            while True:
                b = mode + 0.055*randn()
                if 2.65 <= b <= 3.40: break
            beta = b
            Om   = 0.310 + 0.013*randn() + 0.020*(beta-3.0)
            h    = 0.685 + 0.010*randn() - 0.018*(beta-3.0)
            s8   = 0.795 + 0.017*randn() - 0.025*(beta-3.0)
            wch  = -0.989 + 0.006*randn()
            cD   = 22 + 11*abs(randn())
        S8 = s8 * math.sqrt(Om/0.3)
        samples.append({
            "beta": beta, "Omega_m": Om, "h": h, "sigma8": s8,
            "S8": S8, "w_chi": wch, "c_D": cD
        })
    return samples

def stdv(arr):
    m = sum(arr)/len(arr)
    return math.sqrt(sum((x-m)**2 for x in arr)/len(arr))

def quantile(arr, q):
    s = sorted(arr)
    i = (len(s)-1)*q
    lo, hi = int(math.floor(i)), int(math.ceil(i))
    return s[lo] + (s[hi]-s[lo])*(i-lo)

def stats_1d(samples, key):
    xs = [s[key] for s in samples]
    return {
        "median": quantile(xs, 0.5),
        "q16":    quantile(xs, 0.16),
        "q84":    quantile(xs, 0.84),
        "mean":   sum(xs)/len(xs),
        "std":    stdv(xs),
        "min":    min(xs),
        "max":    max(xs),
    }

def kde2d(samples, xk, yk, nx=50, ny=50):
    xs = [s[xk] for s in samples]
    ys = [s[yk] for s in samples]
    xmn, xmx = min(xs), max(xs)
    ymn, ymx = min(ys), max(ys)
    if xmx-xmn < 1e-4:
        c = (xmx+xmn)/2; xmn, xmx = c-0.04, c+0.04
    if ymx-ymn < 1e-4:
        c = (ymx+ymn)/2; ymn, ymx = c-0.04, c+0.04
    padx = 0.06*(xmx-xmn); pady = 0.06*(ymx-ymn)
    xmin = xmn - padx; xmax = xmx + padx
    ymin = ymn - pady; ymax = ymx + pady
    sx = max(stdv(xs), (xmax-xmin)*0.04, 1e-4)
    sy = max(stdv(ys), (ymax-ymin)*0.04, 1e-4)
    n = len(xs)
    hx = 1.06 * sx * (n ** -0.2)
    hy = 1.06 * sy * (n ** -0.2)
    dx = (xmax-xmin)/(nx-1); dy = (ymax-ymin)/(ny-1)
    grid = [[0.0]*nx for _ in range(ny)]
    rx = max(2, int(math.ceil(3*hx/dx)))
    ry = max(2, int(math.ceil(3*hy/dy)))
    for s_ in range(n):
        sxv = xs[s_]; syv = ys[s_]
        ix = int(round((sxv-xmin)/dx)); iy = int(round((syv-ymin)/dy))
        i0 = max(0, ix-rx); i1 = min(nx-1, ix+rx)
        j0 = max(0, iy-ry); j1 = min(ny-1, iy+ry)
        for j in range(j0, j1+1):
            gy = ymin + j*dy
            ty = (gy-syv)/hy
            ky = math.exp(-0.5*ty*ty)
            if ky < 0.002: continue
            for i in range(i0, i1+1):
                gx = xmin + i*dx
                tx = (gx-sxv)/hx
                grid[j][i] += math.exp(-0.5*tx*tx) * ky
    # round to 4 sig figs to keep JSON small
    def r4(v):
        if v == 0: return 0
        return float(f"{v:.4g}")
    flat = [r4(grid[j][i]) for j in range(ny) for i in range(nx)]
    return {
        "xmin": round(xmin, 5), "xmax": round(xmax, 5),
        "ymin": round(ymin, 5), "ymax": round(ymax, 5),
        "nx": nx, "ny": ny,
        "grid": flat,
    }

PARAMS = ["beta", "Omega_m", "h", "sigma8", "S8", "w_chi", "c_D"]
PRETTY = {"beta":"β","Omega_m":"Ωₘ","h":"h","sigma8":"σ₈","S8":"S₈","w_chi":"w_χ","c_D":"c_D"}

# Parameter pairs of primary interest for the explorer (keep file size in check)
PAIRS = [
    ("beta", "c_D"),
    ("beta", "Omega_m"),
    ("beta", "sigma8"),
    ("beta", "h"),
    ("beta", "w_chi"),
    ("Omega_m", "sigma8"),
    ("Omega_m", "h"),
    ("h", "sigma8"),
    ("Omega_m", "S8"),
    ("c_D", "sigma8"),
    ("c_D", "Omega_m"),
    ("w_chi", "c_D"),
]

R1 = {"3a": 0.018, "3b": 0.034, "3c": 0.081}
COLORS = {"3a": "#155fa5", "3b": "#0f6e56", "3c": "#993c1d"}
LABELS = {
    "3a": "3a baseline (β = 0, arrested-medium)",
    "3b": "3b non-crossing (β ∈ (0, 2.4])",
    "3c": "3c crossing, B1-narrow (β ∈ [2.65, 2.75])",
}

out = {
    "schema_version": 1,
    "generated_utc": "2026-06-02",
    "data_source": "MOCK SAMPLES — replace with export_chains_summary.py output from real chains.",
    "params": PARAMS,
    "pretty": PRETTY,
    "credible_levels": [0.68, 0.95],
    "branches": {},
}

for br in ["3a", "3b", "3c"]:
    samples = gen_branch(br, n=4000)
    print(f"  {br}: {len(samples)} samples")
    branch_data = {
        "label": LABELS[br],
        "color": COLORS[br],
        "R1":    R1[br],
        "n_samples": len(samples),
        "stats_1d": {p: stats_1d(samples, p) for p in PARAMS},
        "kde_2d": {}
    }
    for xk, yk in PAIRS:
        key = f"{xk}__{yk}"
        branch_data["kde_2d"][key] = kde2d(samples, xk, yk, nx=50, ny=50)
    out["branches"][br] = branch_data

dest = "/sessions/pensive-sleepy-newton/mnt/Tom/CLASS_SYMT/chains_summary.json"
with open(dest, "w") as fh:
    json.dump(out, fh, separators=(",", ":"))

import os
print(f"wrote {dest}  ({os.path.getsize(dest)/1024:.1f} KiB)")
