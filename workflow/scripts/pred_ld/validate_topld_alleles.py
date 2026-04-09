"""
Validate REF2/ALT2 alleles in a TOP-LD LD_info file against the Ensembl hg38 REST API.

Samples up to `n_snps` unique rsID-bearing rows from rsID2/REF2/ALT2, queries
Ensembl in batches, and writes a summary to the log output.
"""

import time
import requests
import pandas as pd

ENSEMBL_URL = "https://rest.ensembl.org/variation/homo_sapiens"
BATCH_SIZE = 200
RETRY_AFTER = 5   # seconds to wait on HTTP 429


def query_ensembl_batch(rsids: list[str]) -> dict:
    payload = {"ids": rsids}
    for attempt in range(3):
        resp = requests.post(
            ENSEMBL_URL,
            json=payload,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        if resp.status_code == 200:
            return resp.json()
        if resp.status_code == 429:
            time.sleep(RETRY_AFTER * (attempt + 1))
        else:
            resp.raise_for_status()
    resp.raise_for_status()


def main(ld_info_path: str, log_path: str, n_snps: int, chrom: str):
    df = pd.read_csv(ld_info_path, sep="\t", usecols=["pos2", "rsID2", "REF2", "ALT2"])

    # Keep only rows with real rsIDs; deduplicate on rsID2
    df = df[df["rsID2"].str.startswith("rs")].drop_duplicates(subset="rsID2")

    if len(df) > n_snps:
        df = df.sample(n=n_snps, random_state=42)

    rows = df.to_dict("records")

    results: dict = {}
    for i in range(0, len(rows), BATCH_SIZE):
        batch_ids = [r["rsID2"] for r in rows[i : i + BATCH_SIZE]]
        results.update(query_ensembl_batch(batch_ids))
        time.sleep(0.3)

    match = pos_mismatch = allele_mismatch = not_found = 0
    mismatch_lines: list[str] = []

    for r in rows:
        info = results.get(r["rsID2"])
        if not info:
            not_found += 1
            continue

        mappings = [m for m in info.get("mappings", []) if m.get("seq_region_name") == chrom]
        if not mappings:
            not_found += 1
            continue

        m = mappings[0]
        ensembl_alleles = set(m["allele_string"].split("/"))
        pos_diff = m["start"] - r["pos2"]

        if pos_diff not in (0, 1):
            pos_mismatch += 1
            mismatch_lines.append(
                f"POS_MISMATCH\t{r['rsID2']}\ttopld={r['pos2']}\tensembl={m['start']}"
            )
            continue

        if r["REF2"] in ensembl_alleles and r["ALT2"] in ensembl_alleles:
            match += 1
        else:
            allele_mismatch += 1
            mismatch_lines.append(
                f"ALLELE_MISMATCH\t{r['rsID2']}\ttopld={r['REF2']}/{r['ALT2']}\tensembl={m['allele_string']}"
            )

    total = len(rows)
    with open(log_path, "w") as fh:
        fh.write(f"TOP-LD allele validation against Ensembl hg38\n")
        fh.write(f"Input file : {ld_info_path}\n")
        fh.write(f"Chromosome : {chrom}\n")
        fh.write(f"SNPs sampled: {total}\n\n")
        fh.write(f"Match          : {match} ({100*match/total:.1f}%)\n")
        fh.write(f"Allele mismatch: {allele_mismatch} ({100*allele_mismatch/total:.1f}%)\n")
        fh.write(f"Pos mismatch   : {pos_mismatch} ({100*pos_mismatch/total:.1f}%) [expected for indels: off by 1]\n")
        fh.write(f"Not found      : {not_found} ({100*not_found/total:.1f}%)\n")
        if mismatch_lines:
            fh.write("\nMismatch details:\n")
            for line in mismatch_lines:
                fh.write(line + "\n")


main(
    ld_info_path=snakemake.input[0],
    log_path=snakemake.output[0],
    n_snps=snakemake.params.n_snps,
    chrom=snakemake.wildcards.chr.replace("chr", ""),
)
