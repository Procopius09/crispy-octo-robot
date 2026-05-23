==============================================================
  HOW TO UPLOAD THIS PROJECT TO GITHUB
  Tier 1 Radiation Oncology Citation-Integrity Audit
==============================================================

This file walks you through putting the analysis code on GitHub
so others can find, cite, and reproduce the work.  No programming
experience is required — the web browser method (Option A) takes
about five minutes.


--------------------------------------------------------------
BEFORE YOU START — SAFETY CHECKLIST
--------------------------------------------------------------

The following files contain your personal credentials and must
NEVER be uploaded to GitHub:

  config.sh          (holds your NCBI email and API key)

The .gitignore file in this folder already lists config.sh so
that Git will automatically exclude it.  As long as you follow
the steps below, it will not be uploaded.

The following files do NOT contain credentials and are safe:

  config.example.sh  (blank template — safe to share)
  All .py files      (credentials are read from environment
                      variables, not written in the code)
  journals.json, requirements.txt, README.md, etc.


--------------------------------------------------------------
OPTION A — UPLOAD VIA GITHUB WEBSITE (recommended, no coding)
--------------------------------------------------------------

Step 1: Create a GitHub account (if you do not have one)
  - Go to https://github.com
  - Click "Sign up" and follow the prompts
  - Choose the free tier — it is sufficient for this project

Step 2: Create a new repository
  - After logging in, click the green "New" button
    (top-left, next to your repository list)
  - Repository name: tier1-radonc-citation-audit
    (or any name you prefer)
  - Description: Replication code for Miller & Wrightson,
    "The post-ChatGPT citation fabrication surge is not
    universal across medical specialties", The Lancet 2026
  - Set visibility to Public (so reviewers can access it)
  - Check "Add a README file"
  - Choose MIT License from the dropdown
  - Click "Create repository"

Step 3: Upload the files
  - On your new repository page, click "Add file" → "Upload files"
  - Drag and drop ALL files from this folder EXCEPT config.sh:
      .gitignore
      analysis.py
      analysis_report.md
      audit.py
      build_figure.py
      config.example.sh        ← upload this (it is the blank template)
      eutils.py
      fetch.py
      journals.json
      LICENSE
      README.md
      requirements.txt
      run.sh
      setup.sh
      verify.py
  - DO NOT upload:  config.sh  (your real credentials)
  - DO NOT upload:  data/*.csv  (large data files — use Zenodo)
  - In the "Commit changes" box at the bottom, write:
      "Initial upload of analysis pipeline"
  - Click "Commit changes"

Step 4: Add the Zenodo DOI badge (optional but recommended)
  - Open README.md on GitHub (click the file, then the pencil
    edit icon)
  - At the very top of the file, add this line:
      [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20350450.svg)](https://doi.org/10.5281/zenodo.20350450)
  - Commit the change
  - A clickable DOI badge will now appear on your repository
    front page, linking reviewers directly to the dataset

Step 5: Copy the repository URL
  - The URL will look like:
      https://github.com/YOUR-USERNAME/tier1-radonc-citation-audit
  - Add this URL to the manuscript or supplementary material
    wherever the code is referenced


--------------------------------------------------------------
OPTION B — UPLOAD VIA GIT COMMAND LINE (for users familiar
           with the terminal)
--------------------------------------------------------------

  # One-time setup (if Git is not installed)
  # macOS:  xcode-select --install
  # Ubuntu: sudo apt install git

  cd /path/to/this/folder

  git init
  git remote add origin https://github.com/YOUR-USERNAME/tier1-radonc-citation-audit.git

  # Stage all safe files (.gitignore excludes config.sh and data/*.csv)
  git add .

  # Verify config.sh is NOT staged (this should return nothing)
  git status | grep config.sh

  git commit -m "Initial upload of analysis pipeline"
  git branch -M main
  git push -u origin main


--------------------------------------------------------------
LINKING THE ZENODO DATASET TO THE GITHUB REPOSITORY
--------------------------------------------------------------

The full dataset (328,091 references, ~52 MB) is already
deposited at Zenodo:

  DOI:  10.5281/zenodo.20350450
  URL:  https://doi.org/10.5281/zenodo.20350450

GitHub is for the CODE; Zenodo is for the DATA.  Reviewers
and replicators should:
  1. Clone / download the code from GitHub
  2. Download the data from Zenodo
  3. Place the per-journal CSV files in the data/ folder
  4. Run the analysis stage:  python audit.py --stage analyze


--------------------------------------------------------------
CITATION
--------------------------------------------------------------

When referring to this code in publications or other materials,
use the Zenodo DOI (which covers both the code and data):

  Miller RC, Wrightson T. Tier 1 radiation oncology citation-
  integrity audit — analysis pipeline and dataset. Zenodo, 2026.
  DOI: 10.5281/zenodo.20350450

Or in Vancouver style for The Lancet:

  Miller RC, Wrightson T. Tier 1 radiation oncology citation-
  integrity audit — analysis code and dataset [Internet]. Zenodo;
  2026 [cited 2026 May 23]. Available from:
  https://doi.org/10.5281/zenodo.20350450


--------------------------------------------------------------
QUESTIONS
--------------------------------------------------------------

Contact: Robert C. Miller, MD MBA FRSA
         miller.robert@mayo.edu
         ORCID: 0000-0001-8932-2732

         Tessa Wrightson, DO
         twrightson179@marian.edu
         ORCID: 0009-0006-7191-8830

==============================================================
