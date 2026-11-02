# Publish this project on GitHub

1. In the GitHub account where this project will live, create a **new empty public repository** named `DE_Particle_paper` (do not initialize it with a README, .gitignore or license: the prepared bundle already has README and .gitignore). Confirm all coauthors agree to sharing the 3000 generated runs and notebooks.
2. Extract `DE_Particle_paper_public.zip` and **open a terminal inside the `DE_Particle_paper_public` folder** (with `README.md` and `notebooks/` at top level).
3. Run the following after checking the connected GitHub account and replacing `<GITHUB-OWNER>` with its login:

   ```bash
   git init
   git add .
   git commit -m "Release reproducible 2HDM+S DE reconstruction analysis"
   git branch -M main
   git remote add origin https://github.com/<GITHUB-OWNER>/DE_Particle_paper.git
   git push -u origin main
   ```

   Git needs an authenticated session for the push. In VS Code you can instead use **Source Control → Publish to GitHub** and choose **Public**; make sure the owner and repository name match.
4. Check the GitHub Actions workflow `Verify published data`. It must be green before publishing a DOI or announcing the release.
5. Approve explicit code/data licenses and a final `CITATION.cff` entry, then tag a release and archive the exact tag (e.g. on ScienceDB or Zenodo).

**Alternative:** upload the project files through the GitHub web interface or a GitHub desktop client. Generating this ZIP does not itself create a GitHub repository.
