"""This script is loaded near the beginning of the main phenotyping Snakefile.

It validates and processes the configuration and sets global variables.
"""

import pandas as pd
from pathlib import Path

# === Load dataset-specific settings ===
if "dataset" not in config:
    raise ValueError("Missing 'dataset' in config. Use: --config dataset=your_dataset")
if "datasets" not in config:
    raise KeyError("Missing 'datasets' section in config.yml")

dataset = config["dataset"]
if dataset not in config["datasets"]:
    raise KeyError(f"Dataset '{dataset}' not found in config['datasets']")

dataset_config = config["datasets"][dataset]

# === Pull dataset values into top-level config for consistency ===
config["fastq_dir"] = Path(dataset_config["fastq_dir"]).expanduser()
config["fastq_map"] = Path(dataset_config["fastq_map"]).expanduser()
config["samples_file"] = Path(dataset_config["samples_file"]).expanduser()
config["paired_end"] = dataset_config.get("paired_end", True)
config["read_length"] = dataset_config.get("read_length", 101)

if not config["paired_end"]:
    config["fragment_length_mean"] = dataset_config.get("fragment_length_mean")
    config["fragment_length_sd"] = dataset_config.get("fragment_length_sd")

# === Output directories
config["intermediate_dir"] = Path(f"intermediate/{dataset}")
config["intermediate_ref_dir"] = config["intermediate_dir"] / "reference"
output_dir = Path(f"output/{dataset}")

# === Validation function ===
def validate_config(config: dict):
    if len(config.keys()) == 0:
        raise Exception('No config file provided, or it is empty.')

    for field in ["fastq_dir", "fastq_map", "samples_file", "read_length"]:
        if field not in config:
            raise KeyError(f"Required config field missing: {field}")

    for field in ["ref_genome", "ref_anno", "modality_groups"]:
        if field not in config:
            raise KeyError(f"Required global config field missing: {field}")

    if not config["ref_anno"].lower().endswith(".gtf"):
        raise ValueError("ref_anno must be an uncompressed GTF file")

    if not isinstance(config["read_length"], int):
        raise ValueError("read_length must be an integer")

    if "paired_end" in config and not isinstance(config["paired_end"], bool):
        raise ValueError("paired_end must be True or False")

    if config["paired_end"] is False:
        if "fragment_length_mean" not in config:
            raise KeyError("fragment_length_mean required for single-end data")
        if "fragment_length_sd" not in config:
            raise KeyError("fragment_length_sd required for single-end data")

# === Load fastq_map file ===
def load_fastq_map(map_file: Path, fastq_dir: Path, paired_end=None) -> dict:
    paths = {}
    with open(map_file, "r") as f:
        for line in f.read().splitlines():
            fields = line.split("\t")
            if len(fields) == 4:
                r1, r2, sample_id, read_len = fields
                is_paired = True
            elif len(fields) == 3:
                r1, sample_id, read_len = fields
                is_paired = False
            else:
                raise ValueError(f"Unexpected number of fields in line: {line}")

            if sample_id not in paths:
                paths[sample_id] = (([], []), is_paired, int(read_len)) if is_paired else ([], is_paired, int(read_len))
            elif paths[sample_id][1] != is_paired:
                raise ValueError(f"Sample {sample_id} has mixed paired-end and single-end reads")

            if is_paired:
                paths[sample_id][0][0].append(str(Path(r1)))
                paths[sample_id][0][1].append(str(Path(r2)))
            else:
                paths[sample_id][0].append(str(Path(r1)))

    if paired_end is not None:
        for sample_id, (_, is_paired, _) in paths.items():
            if is_paired != paired_end:
                raise ValueError(f"Sample {sample_id} has {'paired' if is_paired else 'single'}-end reads, "
                                f"but config specifies {'paired' if paired_end else 'single'}-end")

    return paths

# === Process the config
def process_config(config: dict):
    for path_key in ["ref_genome", "ref_anno", "samples_file"]:
        config[path_key] = Path(config[path_key]).expanduser()

    samples = pd.read_csv(config["samples_file"], sep="\t", header=None, dtype=str)[0].tolist()
    config["samples"] = samples

    full_fastq_map = load_fastq_map(
        config["fastq_map"],
        config["fastq_dir"],
        config["paired_end"]
    )

    # Save full info separately
    config["fastq_map_full"] = full_fastq_map

    # Return only (paths, is_paired) for compatibility with Snakefile
    config["fastq_map"] = {k: (v[0], v[1]) for k, v in full_fastq_map.items()}

    # Store per-sample read lengths
    config["read_lengths"] = {k: v[2] for k, v in full_fastq_map.items()}

# === Reference validation
def validate_reference(ref_genome: Path, ref_anno: Path):
    if not ref_genome.exists():
        raise FileNotFoundError(f"Reference genome file not found: {ref_genome}")
    if not ref_anno.exists():
        raise FileNotFoundError(f"Reference annotation file not found: {ref_anno}")

    with open(ref_anno) as f:
        lines = [next(f) for _ in range(1000)]

    has_biotype = any("transcript_biotype" in line for line in lines)
    has_ensembl_source = any('transcript_source "ensembl' in line for line in lines)

    if not has_biotype:
        raise ValueError(f"GTF file {ref_anno} missing 'transcript_biotype' field (likely not from Ensembl)")
    if not has_ensembl_source:
        raise ValueError(f"GTF file {ref_anno} missing Ensembl transcript source")

# === Execute
validate_config(config)
process_config(config)

interm_dir = config["intermediate_dir"]
ref_dir = config["intermediate_ref_dir"]
output_dir = Path(f"output/{dataset}")

samples_file = config["samples_file"]
samples = config["samples"]

paired_end = config["paired_end"]
read_length = config["read_length"]
fastq_dir = config["fastq_dir"]
fastq_map = config["fastq_map"]

if not paired_end:
    fragment_length_mean = config["fragment_length_mean"]
    fragment_length_sd = config["fragment_length_sd"]

ref_genome = config["ref_genome"]
ref_anno = config["ref_anno"]
ref_cdna = ref_dir / "cDNA.fa.gz"

validate_reference(ref_genome, ref_anno)

modality_groups = config["modality_groups"]


if 'RNA_editing' in modality_groups:
    edit_sites_bed = config['edit_sites_bed']
    edit_sites_min_coverage = config['edit_sites_min_coverage']
    edit_sites_min_samples = min(config['edit_sites_min_samples'], len(samples))
    edit_sites_bed = Path(edit_sites_bed).expanduser()
    if not edit_sites_bed.exists():
        raise FileNotFoundError(f'Edit sites BED file not found: {edit_sites_bed}')
    if edit_sites_min_samples > len(samples):
        print(f'Note: edit_sites_min_samples ({edit_sites_min_samples}) is greater than the number of samples ({len(samples)}). Setting to {len(samples)}.')
        edit_sites_min_samples = len(samples)


outputs = []
for modality_group, params in modality_groups.items():
    for f in params["files"]:
        outputs.append(output_dir / f)

### END