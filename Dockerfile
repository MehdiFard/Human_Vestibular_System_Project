FROM rocker/tidyverse:4.4

# System dependencies for Python venv
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    # libcblas3 \
    && rm -rf /var/lib/apt/lists/*

# Python virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:/usr/local/bin:$PATH"

# Install Python packages
RUN pip install --no-cache-dir \
    requests \
    pandas \
    numpy \
    scipy \
    "gql[aiohttp]" \
    git+https://github.com/twillis209/PRED-LD.git@pip

# Install LDSC from fork with Python 3 compatibility fixes applied.
# Pin to a specific commit so that updating this SHA both triggers a rebuild
# and busts the Docker layer cache for this step.
RUN git clone https://github.com/twillis209/ldsc.git /opt/ldsc \
    && git -C /opt/ldsc checkout 934a94d052fb4c12b713b79f7ffe445b7e303645 \
    && pip install --no-cache-dir bitarray

# Install pak for fast parallel R package installation
RUN Rscript -e "install.packages('pak', repos = 'https://cloud.r-project.org')"

# Bioconductor and CRAN packages (stable — add new ones here)
RUN Rscript -e "\
    pak::pak(c( \
        'AnnotationHub', \
        'ensembldb', \
        'GenomicRanges', \
        'data.table', \
        'magrittr', \
        'ggtext', \
        'ggrepel', \
        'R.utils' \
    ))"

# GitHub R packages (volatile — add new ones here without touching the layer above)
RUN Rscript -e "\
    pak::pak(c( \
        'twillis209/tomics', \
        'ajbass/sffdr' \
    ))"

# plink2 and plink binaries
RUN wget -q -O /tmp/plink2.zip \
        https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20260311.zip \
    && unzip -j /tmp/plink2.zip plink2 -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink2 \
    && rm /tmp/plink2.zip \
    && wget -q -O /tmp/plink.zip \
        https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip \
    && unzip -j /tmp/plink.zip plink -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink \
    && rm /tmp/plink.zip
