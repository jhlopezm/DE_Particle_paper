# DE_Particle_paper

**Reproducible inverse reconstruction of a cubic 2HDM+complex-singlet scalar sector with Differential Evolution.**

This repository contains the **actual 3000 completed DE optimizations** (1000 mass-target triplets × 3 independent restarts), a frozen computational model and public-ready Jupyter notebooks to reproduce the post-processing and inspect representative numerical preimages.

> **Scope:** These are empirical inverse-reconstruction ensembles, **not posterior samples**. The final post-run selection imposes **five strict one-width criteria**. The widths σ are numerical reconstruction scales rather than independently established experimental errors. The source data and their original optimizer status are preserved unchanged.

## Verified results

| Test | Count |
| --- | ---: |
| Runs with all three mass residuals below one width | 2977 / 3000 |
| Charged-scalar residual below one width | 2994 / 3000 |
| Vector-coupling residual below one width | 3000 / 3000 |
| **All five residuals below one width** | **2971 / 3000** |
| Targets with at least one accepted reconstruction | **1000 / 1000** |

The forward model uses `v = 246.22 GeV`, with CP-even eigenvectors as **columns** (`R.T @ M2 @ R = diag(m_i²)`); the H2 sign is fixed by the nonnegative reduced vector coupling. The base optimizer seeds are **11, 23, 37**, and the target-generation seed is **20261020**.

## Start here (VS Code / Jupyter)

1. Clone the repository or download and extract the ZIP. Open its **root folder**, which contains `data/`, `notebooks/`, and `src/`.
2. Use **Python 3.12** (the production environment used Python 3.12.1). Create a virtual environment, activate it, and install the dependencies:

   ```bash
   python -m venv .venv
   # Windows PowerShell: .venv\Scripts\Activate.ps1
   # Linux / macOS: source .venv/bin/activate
   python -m pip install -r requirements.txt
   ```
3. Choose the `.venv` Python kernel in VS Code or Jupyter; start with [`notebooks/00_start_here.ipynb`](notebooks/00_start_here.ipynb).
4. Run the independent audit: `python scripts/verify_dataset.py`. This **does not rerun DE**.

### Google Colab

Create a notebook at [colab.research.google.com](https://colab.research.google.com), and clone this repository using a cell:

```python
!git clone https://github.com/<OWNER>/DE_Particle_paper.git
%cd DE_Particle_paper
!pip install -r requirements.txt
```

Then open an individual notebook from GitHub in Colab (File → Open notebook → GitHub) and ensure that the working directory points to the cloned project root. The public notebooks deliberately **do not download unknown files at runtime** or contain API credentials. Replace `<OWNER>` once the repo URL is known. Colab's runtime Python version may differ from the pinned production version; the notebooks reproduce the *saved-data analysis*, not bit-for-bit DE reruns.

## Notebook guide

| Notebook | What it does |
| --- | --- |
| [`00_start_here.ipynb`](notebooks/00_start_here.ipynb) | Validate environment, data and one real preimage |
| [`01_selection_audit.ipynb`](notebooks/01_selection_audit.ipynb) | Recalculate all five one-width cuts and check target coverage |
| [`02_numerical_preimages.ipynb`](notebooks/02_numerical_preimages.ipynb) | Print parameters, CP-even mass matrix, eigendecomposition and predicted masses |
| [`03_degeneracy_distributions.ipynb`](notebooks/03_degeneracy_distributions.ipynb) | Analyze same-target degeneracy, representative vs. all accepted solutions |
| [`04_correlation_validation.ipynb`](notebooks/04_correlation_validation.ipynb) | Exact identities, eigenvector-composition checks, conditional approximations |
| [`05_protocol_and_demo.ipynb`](notebooks/05_protocol_and_demo.ipynb) | Describe original DE configuration and optional tiny illustrative run |

Saved small outputs make the notebooks legible on GitHub without execution; all notebooks can be **Run All** in a normal Jupyter kernel from either the repo root or its `notebooks/` directory.

## Data and implementation

- [`data/raw`](data/raw): unmodified final `runs.csv` (3000 rows), target draws (1000 rows), original optimizer configuration, full production environment snapshot (text transcoded to UTF-8 for GitHub readability).
- [`data/processed`](data/processed): verified five-one-width subsets, one best accepted run per target, summary JSON.
- [`src/model.py`](src/model.py): canonical mass matrix, diagonalization and numerical observables.
- [`src/run_de.py`](src/run_de.py): resumable 3-seed DE **production script**. A full run is computationally expensive. The script records original optimizer constraints, which differ from the stricter final post-run cuts.
- [`figures`](figures): publication-quality analysis figures calculated from the saved data.
- [`docs/methods.md`](docs/methods.md): detailed method and acceptance definitions; [`docs/data_dictionary.md`](docs/data_dictionary.md): CSV schema.

`python scripts/verify_dataset.py` recomputes every mass, charged scalar and κ from saved parameter vectors, compares against their recorded values, and asserts exact final counts.

## Citing and releasing

This is a **prepublication research dataset**, prepared for coauthor review. The repository contains a provisional [`CITATION.cff`](CITATION.cff), but no article DOI is claimed. Add a persistent Science Data Bank or Zenodo DOI and final publication metadata once they exist. Authors should explicitly approve code/data licensing before publishing; see [`LICENSE_NOTICE.md`](LICENSE_NOTICE.md).
