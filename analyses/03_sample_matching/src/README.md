# 03_sample_matching/src (not in git)

- `path_to_fqs.tsv`: one row per RNA-seq library: R1, R2 (paths under data/rna_fq), read_length, library_id, human_id
  (donor), batch, tissue. Participant-level sample sheet: stays on the server.
- `Homo_sapiens.GRCh38.110.chr.collapsed.genes.gtf`: Ensembl 110 GTF collapsed to one transcript per gene (GTEx
  collapse_annotation.py; rule `collapse_annotation`)
- `GRCh38_exon_regions.tsv.gz`: exon regions of the collapsed annotation, for allele counting
