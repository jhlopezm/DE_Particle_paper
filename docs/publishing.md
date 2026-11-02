# Publishing and archiving DE_Particle_paper

The GitHub repository is already public at **https://github.com/jhlopezm/DE_Particle_paper**. The `main` branch contains the six Jupyter notebooks, the archived 3000-run dataset and analysis code; no new repository needs to be created.

## Licenses

- Code and notebooks: [MIT](../LICENSE).
- Raw and processed datasets: [CC BY 4.0](../data/LICENSE).
- See [../LICENSE_NOTICE.md](../LICENSE_NOTICE.md) for exact scopes and recommended attribution.

## Next publication steps

1. Confirm that all three research authors are comfortable with the public release and proposed attribution and have authority to release their contributions.
2. Keep `python scripts/verify_dataset.py` and GitHub Actions green, and verify the paper's reported 5×1σ counts against the release files.
3. Update `CITATION.cff` with the final manuscript and persistent dataset DOI when available.
4. Create an immutable GitHub release tag such as `v1.0.0` after coauthor approval; archive that **specific tag** (not an evolving branch) in Science Data Bank or Zenodo.
5. Insert the resulting DOI into the paper's data-and-software availability statement.

For later updates from VS Code, use a normal incremental Git commit and `git push origin main`; never force-push over the published history.
