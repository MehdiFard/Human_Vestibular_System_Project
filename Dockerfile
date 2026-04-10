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
ENV PATH="/opt/venv/bin:$PATH"

# Install Python packages
RUN pip install --no-cache-dir \
    requests \
    pandas \
    numpy \
    scipy \
    "gql[aiohttp]" \
    git+https://github.com/twillis209/PRED-LD.git@pip

# Install LDSC (Python 3.9+ branch) and its dependencies
RUN git clone -b ldsc39 https://github.com/CBIIT/ldsc.git /opt/ldsc \
    && pip install --no-cache-dir -r /opt/ldsc/requirements.txt bitarray

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

