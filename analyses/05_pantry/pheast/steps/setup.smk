"""This script is loaded near the beginning of the main Pheast Snakefile.
It validates and processes the configuration and sets global variables.
"""

from pathlib import Path
import pandas as pd

# === Require --config dataset=NAME ===
if "dataset" not in config:
    raise ValueError("Missing required parameter: dataset. Use --config dataset=YOUR_DATASET_NAME")

dataset = config["dataset"]

def validate_config(config: dict):
    if "datasets" not in config:
        raise ValueError("Missing 'datasets' block in config file.")
    if dataset not in config["datasets"]:
        raise ValueError(f"Dataset '{dataset}' not found in config['datasets']")
    if "modalities" not in config or "analyses" not in config:
        raise ValueError("Missing 'modalities' or 'analyses' in config.")

def process_config(config: dict):
    """Load dataset-specific paths and sample list."""
    dataset_cfg = config["datasets"][dataset]

    # Required fields from the dataset block
    required_fields = ["phenotype_dir", "geno_prefix", "samples_file", "age_sex_info"]
    for field in required_fields:
        if field not in dataset_cfg:
            raise ValueError(f"Missing field '{field}' in datasets -> {dataset}")

    config["phenotype_dir"] = Path(dataset_cfg["phenotype_dir"]).expanduser()
    config["geno_prefix"] = Path(dataset_cfg["geno_prefix"]).expanduser()
    config["samples_file"] = Path(dataset_cfg["samples_file"]).expanduser()
    config["age_sex_info"] = Path(dataset_cfg["age_sex_info"]).expanduser()

    # Output directories
    config["intermediate_dir"] = Path("intermediate") / dataset
    config["output_dir"] = Path("output") / dataset

    # Read sample IDs
    config["samples"] = pd.read_csv(config["samples_file"], sep='\t', header=None, dtype=str)[0].tolist()

# Run setup
validate_config(config)
process_config(config)

# Assign short aliases
pheno_dir = config["phenotype_dir"]
interm_dir = config["intermediate_dir"]
output_dir = config["output_dir"]
geno_prefix = config["geno_prefix"]
samples = config["samples"]
age_sex_info = config["age_sex_info"]


# Check if all samples exist in .fam
geno_samples = pd.read_csv(str(geno_prefix) + '.fam', sep=r'\s+', header=None, dtype=str)[1].tolist()
missing_samples = [s for s in samples if s not in geno_samples]
if missing_samples:
    raise Exception(f"Samples missing from genotype FAM file: {missing_samples}")

# Add chromosome info to modalities
modalities = config["modalities"]
for modality, params in modalities.items():
    bed_file = pheno_dir / f"{modality}.bed.gz"
    if not bed_file.exists():
        raise FileNotFoundError(f"Missing BED file for modality '{modality}': {bed_file}")
    chroms = pd.read_csv(bed_file, sep='\t', usecols=[0], dtype=str)
    params["chroms"] = chroms.iloc[:, 0].unique().tolist()

# Expand expected output files
analyses = config["analyses"]
outputs = []
for analysis, params in analyses.items():
    for f in params["files"]:
        outputs.append(expand(output_dir / analysis / f, modality=modalities.keys()))

### END