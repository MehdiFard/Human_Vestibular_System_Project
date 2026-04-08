FROM rocker/tidyverse:4.4

# System dependencies for Python venv
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Python virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python packages
RUN pip install --no-cache-dir \
    requests \
    pandas \
    "gql[aiohttp]"

# Install Bioconductor packages
RUN Rscript -e "\
    options(repos = BiocManager::repositories()); \
    BiocManager::install(c('AnnotationHub', 'ensembldb', 'GenomicRanges'), \
        ask = FALSE, update = FALSE)"

# Install CRAN packages not already in rocker/tidyverse
RUN Rscript -e "\
    install.packages(c('data.table', 'magrittr', 'ggtext', 'ggrepel', 'R.utils'), \
        repos = 'https://cloud.r-project.org')"

# Install GitHub R packages
RUN Rscript -e "\
    remotes::install_github('twillis209/tomics', ref = 'master', upgrade = 'never'); \
    remotes::install_github('ajbass/sffdr', ref = 'main', upgrade = 'never')"
