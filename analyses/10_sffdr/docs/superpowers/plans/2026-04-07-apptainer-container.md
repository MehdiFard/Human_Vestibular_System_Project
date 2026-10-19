# Apptainer Container Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace conda dependency management with a single Apptainer-compatible Docker image hosted on GHCR, integrated into the Snakemake pipeline via the `container:` directive.

**Architecture:** A single `Dockerfile` based on `rocker/tidyverse:4.4` installs all R (CRAN, Bioconductor, GitHub) and Python dependencies. A GitHub Actions workflow builds and pushes the image to GHCR on changes to the `Dockerfile`. The Snakemake pipeline is updated to reference the image via a global `container:` directive.

**Tech Stack:** Docker, Apptainer, Snakemake 9+, R 4.4 (rocker/tidyverse), BiocManager, remotes, Python 3, pip, GitHub Actions, GHCR

---

### Task 1: Write the Dockerfile

**Files:**
- Create: `Dockerfile`

- [ ] **Step 1: Create `Dockerfile`**

```dockerfile
FROM rocker/tidyverse:4.4

# System dependencies for Python and gql's aiohttp transport
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor packages
RUN Rscript -e "\
    options(repos = BiocManager::repositories()); \
    BiocManager::install(c('AnnotationHub', 'ensembldb', 'GenomicRanges'), \
        ask = FALSE, update = FALSE)"

# Install CRAN packages not already in rocker/tidyverse
RUN Rscript -e "\
    install.packages(c('data.table', 'magrittr', 'ggtext', 'ggrepel'), \
        repos = 'https://cloud.r-project.org')"

# Install GitHub R packages
RUN Rscript -e "\
    remotes::install_github('twillis209/tomics', upgrade = 'never'); \
    remotes::install_github('ajbass/sffdr', upgrade = 'never')"

# Install Python packages
RUN pip3 install --no-cache-dir --break-system-packages \
    requests \
    pandas \
    "gql[aiohttp]"
```

- [ ] **Step 2: Verify the Dockerfile builds locally (optional but recommended before pushing)**

```bash
docker build -t sffdr_for_dizziness:test .
```

Expected: Build completes without errors. This may take 10-20 minutes on first run due to Bioconductor package compilation.

If `docker` is not available locally (HPC-only workflow), skip to Task 2 and let CI validate the build.

- [ ] **Step 3: Commit**

```bash
git add Dockerfile
git commit -m "feat: add Dockerfile with R/Bioconductor/Python dependencies

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Write the GitHub Actions workflow

**Files:**
- Create: `.github/workflows/docker.yml`

The workflow:
- Triggers on push to `main` when `Dockerfile` changes, and on manual dispatch
- Logs in to GHCR using the built-in `GITHUB_TOKEN`
- Builds and pushes two tags: `latest` and `sha-<short-commit>`

- [ ] **Step 1: Create `.github/workflows/docker.yml`**

```yaml
name: Build and push Docker image

on:
  push:
    branches:
      - main
    paths:
      - Dockerfile
  workflow_dispatch:

jobs:
  build-and-push:
    runs-on: ubuntu-latest
    permissions:
      contents: read
      packages: write

    steps:
      - name: Checkout
        uses: actions/checkout@v4

      - name: Log in to GHCR
        uses: docker/login-action@v3
        with:
          registry: ghcr.io
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - name: Extract metadata
        id: meta
        uses: docker/metadata-action@v5
        with:
          images: ghcr.io/${{ github.repository }}
          tags: |
            type=raw,value=latest,enable={{is_default_branch}}
            type=sha,prefix=sha-,format=short

      - name: Build and push
        uses: docker/build-push-action@v5
        with:
          context: .
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
```

- [ ] **Step 2: Commit**

```bash
git add .github/workflows/docker.yml
git commit -m "ci: add GitHub Actions workflow to build and push Docker image to GHCR

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

- [ ] **Step 3: Push to GitHub and verify the workflow runs**

```bash
git push
```

Navigate to the repository's **Actions** tab on GitHub. The "Build and push Docker image" workflow should trigger. Verify it completes successfully and the image appears at:

```
https://github.com/<your-username>/<your-repo>/pkgs/container/<your-repo>
```

---

### Task 3: Integrate the container into the Snakefile

**Files:**
- Modify: `workflow/Snakefile`

Add a global `container:` directive so every rule uses the GHCR image when `--use-apptainer` is passed. The directive goes near the top of the file, after the `min_version` call and imports.

- [ ] **Step 1: Edit `workflow/Snakefile` — add the `container:` directive**

Insert after the `min_version("9.0")` line and the imports block, before the first rule:

```python
container: "docker://ghcr.io/<your-github-username>/<your-repo-name>:latest"
```

Replace `<your-github-username>/<your-repo-name>` with the actual value, e.g.:
```python
container: "docker://ghcr.io/twillis209/sffdr_for_dizziness:latest"
```

The relevant section of `workflow/Snakefile` should look like:

```python
from snakemake.utils import min_version
min_version("9.0")

import re
import os
import pandas as pd

container: "docker://ghcr.io/twillis209/sffdr_for_dizziness:latest"

def script_path(script):
    return os.path.join(workflow.basedir, "scripts", script)
```

- [ ] **Step 2: Verify the directive is parsed correctly**

```bash
cd /path/to/sffdr_for_dizziness
snakemake --use-apptainer --dry-run
```

Expected: Snakemake prints the dry-run plan without errors. It will note it would use Apptainer for rules that have container support.

- [ ] **Step 3: Add a note about AnnotationHub cache bind-mount**

The `download_ensembl_db` rule uses `AnnotationHub`, which writes a local cache on first run. This cache directory must be bind-mounted into the container so it persists across runs. Add a comment to the rule:

```python
# NB: Needs to create cache when running AnnotationHub for the first time.
# Pass --apptainer-args "--bind $HOME/.cache/R/AnnotationHub:$HOME/.cache/R/AnnotationHub"
# (or your cluster's equivalent) to persist the AnnotationHub cache across runs.
rule download_ensembl_db:
```

- [ ] **Step 4: Commit**

```bash
git add workflow/Snakefile
git commit -m "feat: integrate GHCR container into Snakefile via container: directive

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Running the pipeline

Once the image is pushed to GHCR, run the pipeline with:

```bash
snakemake --use-apptainer \
    --apptainer-args "--bind /path/to/data:/path/to/data --bind $HOME/.cache/R/AnnotationHub:$HOME/.cache/R/AnnotationHub" \
    --cores <N>
```

Apptainer will pull the image from GHCR on first run and cache it as a `.sif` file locally.
