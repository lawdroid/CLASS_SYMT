# Phase 3 MCMC — raw chains + interactive viewer

Folder contents for the Martin–Koh Yang-Mills afterglow framework, Phase 3.
Three pre-registered MCMC branches against Planck NPIPE CamSpec TTTEEE +
Planck 2018 lowl.TT + lowl.EE_sroll2 + lensing.native + DESI DR2 BAO +
Pantheon+ SN.

**Source code:** https://github.com/ingyukoh/CLASS_SYMT (branch
`phase3-interactive`).
**Drive folder (this folder):** https://drive.google.com/drive/folders/1NTaQAoOLeKXfCSWkLLYfu3LSLrkwppud

## Headline (1 minute)

- **β is unconstrained by the data** in both 3b and 3c — the posterior shape
  equals the prior shape in each window. Two independent methods (marginal
  shape, Bayesian evidence) now agree.
- **Bayesian evidence cross-check passes**: Δ log Z (3b − 3a) = −0.18 ± 0.05,
  consistent with the ≈ 0 predicted from a flat-β likelihood.
- **c_D is bounded from above** at c_D ≲ 70 (95%) across all three branches.
- **w_χ ≈ −0.99** in all three branches — effectively a cosmological constant.
- All three branches are fully converged (R−1 ≪ 0.01).

## How to view

Two ways to look at the same data:

1. **Interactive HTML viewer (offline-capable, single file):**
   Download `phase3_summary_standalone.html` and open it in any modern browser.
   No server needed — chain summary is embedded as inline JSON. Banner is
   green, three placeholder slots populate with the evidence numbers, status
   line auto-renders the consistency interpretation.

2. **Live page from the repo:** clone the GitHub repo, `cd CLASS_SYMT`,
   `python -m http.server 8000`, browse to
   `http://localhost:8000/phase3_summary.html`. Identical content, served
   page. The `view_phase3.ps1` script at the repo root does this in one
   command with a confirmation prompt.

## Folder layout

```
Tom3abc/
    Phase3_Drive_README.md                    this file
    phase3_summary_standalone.html             interactive viewer (open in browser)
    chains_summary.json                       265 KB, the data the viewer reads
    evidence_knn_results.json                  raw k-NN evidence output
    Phase3_Evidence_Summary_2026-06-02.md     detailed write-up of the evidence result
    MANIFEST.sha256                            checksums for every raw chain file
    upload_raw_chains.ps1                      run from Windows to upload raw .txt files
                                                                    (see "Raw chains" below)

    phase3a_baseline/                          β = 0, arrested-medium hypothesis
        afterglow_3a.input.yaml                Cobaya input
        afterglow_3a.updated.yaml              Cobaya as-run config (post any patches)
        afterglow_3a.covmat                    proposal covariance from convergence
        afterglow_3a.checkpoint
        afterglow_3a.progress                  R-1 history
        [afterglow_3a.[1-4].txt]               raw chains — see "Raw chains" below

    phase3b_noncrossing/                       β ∈ (0, 2.4]
        ...
        [afterglow_3b.[1-4].txt]

    phase3c_b1narrow/                          β ∈ [2.65, 2.75]
        ...
        [afterglow_3c.[1-4].txt]
```

## Per-branch numbers (extracted from the chains)

| Branch | β prior | n_weighted | R−1 | log Z | β median | c_D median | w_χ median |
|---|---|---|---|---|---|---|---|
| 3a | β = 0 | 505,327 | 0.0013 | −6219.38 ± 0.05 | (fixed) | 22.7 | −0.985 |
| 3b | (0, 2.4] | 764,257 | 0.0002 | −6219.56 ± 0.02 | 1.12 (flat) | 24.4 | −0.986 |
| 3c | [2.65, 2.75] | 1,250,840 | 0.0004 | −6219.90 ± 0.03 | 2.70 | 28.6 | −0.988 |

Cross-checks:
- Δ log Z (3b − 3a) = **−0.18 ± 0.05** — consistent with zero on the
  |ΔlogZ| < 1 Jeffreys "not worth a mention" tier.
- Δ log Z (3c − 3a) = −0.52 ± 0.06 — barely worth a mention, consistent
  with a small Occam penalty for 3c's extra β parameter.

Method: k-NN density evidence estimator (Heavens et al. 2017,
arXiv:1704.03472), clean-room implementation in
`mcmc/evidence_knn.py`. Cross-check absolute log Z values to ±0.5 nat
systematic; deltas are more reliable.

## Raw chains

The Cobaya chain `.txt` files are the actual MCMC samples. Sizes:

| Branch | 4 .txt files total | Avg per file |
|---|---|---|
| 3a | ~540 MB | ~135 MB |
| 3b | ~840 MB | ~210 MB |
| 3c | ~1.4 GB | ~340 MB |

Total ~2.8 GB — too large for normal HTTP uploads. To upload to this folder
from your Windows PC where the chains live:

```powershell
cd C:\Users\ingyu\upwork\Tom\CLASS_SYMT
# Option 1 (easiest): install Google Drive for desktop, sync the Drive folder,
#                     then File Explorer copy the chains/ subdirs in.
# Option 2 (scripted): use rclone (one-time OAuth setup, then):
.\upload_raw_chains.ps1
```

Cobaya `.txt` chain format: each row is
`weight, -log(L·π), <sampled params>, <derived params>`. Header row
(prefixed `#`) names the columns. Read with `np.loadtxt` or `getdist.MCSamples`.

## Reproducing the analysis

```bash
git clone https://github.com/ingyukoh/CLASS_SYMT.git
cd CLASS_SYMT
git checkout phase3-interactive

# To rebuild chains_summary.json from raw .txt files:
python mcmc/export_chains_summary.py

# To rebuild the k-NN evidence numbers:
python mcmc/evidence_knn.py

# To launch the interactive page locally:
.\view_phase3.ps1     # Windows
python -m http.server 8000 && open http://localhost:8000/phase3_summary.html
```

## Pending items

- Δ log Z (3b − ΛCDM) — requires a fresh ΛCDM chain under the matched
  likelihood stack (~12–24 hours wall time on the Ubuntu workstation).
  Not blocking for the consistency claim; needed only if the framing
  sharpens to "preferred/disfavored vs ΛCDM."

## Contact

Ingyu Koh — ingyukoh2@gmail.com. Repo: https://github.com/ingyukoh/CLASS_SYMT.
