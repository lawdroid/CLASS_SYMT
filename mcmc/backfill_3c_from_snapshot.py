"""
Stopgap when the Mac mini is offline and the 3c chain can't be re-fetched.

Reads the May 2026 cached 1D-stats snapshot at
  ../kde_three_branch_2026-05-20.json
and injects a 3c entry into chains_summary.json containing:
  - stats_1d for beta, H0/h, Omega_m, sigma8, S8 (propagated), w_chi, c_D
  - n_samples (= n_weighted from snapshot)
  - data_source flag clearly marking 3c as a snapshot, not live

kde_2d entries for 3c are NOT generated (no samples on hand). The explorer
page should fall back to the static PNG at analysis_2026-05-27/kde_2d_3c.png
for 3c until the chain is re-fetched.

Re-run mcmc/export_chains_summary.py after the Mac mini is back to replace
the snapshot block with real live data.

Run from CLASS_SYMT/ root:
  python mcmc/backfill_3c_from_snapshot.py
"""

import json, math, os, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TOM  = REPO.parent  # /mnt/c/Users/ingyu/upwork/tom or C:\Users\ingyu\upwork\Tom
SNAPSHOT = TOM / "kde_three_branch_2026-05-20.json"
SUMMARY  = REPO / "chains_summary.json"

if not SNAPSHOT.exists():
    print(f"ERROR: snapshot not found at {SNAPSHOT}", file=sys.stderr)
    sys.exit(1)
if not SUMMARY.exists():
    print(f"ERROR: chains_summary.json not found — run export_chains_summary.py first",
          file=sys.stderr)
    sys.exit(1)

snap = json.load(SNAPSHOT.open())
s3c  = snap["summary"]["3c"]

# stats_1d block in the export schema expects:
#   median, q16, q84, mean, std, min, max
# The snapshot only gives mean, q16, median, q84.
# Fill std as a robust approx from (q84-q16)/2 (one-sigma half-width under Gaussian).
# Leave min/max as q16/q84 fallbacks — the explorer plot ranges aren't critical
# without samples to draw.
def expand(d, scale=1.0, offset=0.0):
    m   = d["mean"]   * scale + offset
    q16 = d["q16"]    * scale + offset
    q84 = d["q84"]    * scale + offset
    med = d["median"] * scale + offset
    return {
        "median": med, "q16": q16, "q84": q84,
        "mean": m, "std": max((q84-q16)/2.0, 1e-6),
        "min": q16, "max": q84,
    }

stats = {
    "beta":    expand(s3c["beta_aft"]),
    "Omega_m": expand(s3c["Omega_m"]),
    "h":       expand(s3c["H0"], scale=0.01),
    "sigma8":  expand(s3c["sigma8"]),
    "w_chi":   expand(s3c["w_X"]),
    "c_D":     expand(s3c["c_D"]),
}
# Propagate S8 = sigma8 * sqrt(Omega_m/0.3) on the medians/means; the percentile
# bounds are a rough but defensible composition.
def prop_s8(s8_d, Om_d):
    f = lambda s8, Om: s8 * math.sqrt(Om / 0.3)
    med = f(s8_d["median"], Om_d["median"])
    mn  = f(s8_d["mean"],   Om_d["mean"])
    q16 = f(s8_d["q16"],    Om_d["q16"])
    q84 = f(s8_d["q84"],    Om_d["q84"])
    return {"median": med, "q16": q16, "q84": q84,
            "mean": mn, "std": max((q84-q16)/2.0, 1e-6),
            "min": q16, "max": q84}
stats["S8"] = prop_s8(s3c["sigma8"], s3c["Omega_m"])

summary = json.load(SUMMARY.open())

summary["branches"]["3c"] = {
    "label":      "3c crossing, B1-narrow (β ∈ [2.65, 2.75]) — snapshot 2026-05-20",
    "color":      "#993c1d",
    "R1":         None,
    "n_samples":  s3c["n_weighted"],
    "stats_1d":   stats,
    "kde_2d":     {},  # no samples on hand — explorer falls back to static PNG
    "snapshot_source": str(SNAPSHOT.name),
    "snapshot_date":   snap["date"],
}

# Update data_source so the banner is honest about the mixed provenance.
base = summary["data_source"]
if "snapshot" not in base.lower():
    summary["data_source"] = (
        "REAL chains for 3a, 3b; 3c stats from 2026-05-20 snapshot "
        "(Mac mini offline — re-fetch when reachable)"
    )

with open(SUMMARY, "w") as fh:
    json.dump(summary, fh, separators=(",", ":"))

print(f"backfilled 3c stats into {SUMMARY.name}")
print(f"  3c n_weighted: {s3c['n_weighted']:,}")
print(f"  3c snapshot date: {snap['date']}")
print(f"  data_source: {summary['data_source']}")
print()
print("Reload http://localhost:8000/phase3_summary.html to see the updated banner.")
print("When the Mac mini is back, run:")
print("    python mcmc/export_chains_summary.py")
print("to replace the snapshot block with live 3c data.")
