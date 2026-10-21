# phenotyping/input (not in git; Pantry layout)

- `ref/`: links to the reference files: Ensembl 110 genome FASTA (+ .fai) and GTF (`../../../../resources`), `tx2gene.tsv`
  and the RNA editing sites (`../../src`)
- `<dataset>/` (vestibular = utricle + saccule; utricle; saccule): `<dataset>.samples.tsv` (sample IDs),
  `<dataset>.fastq_map.tsv` (FASTQ -> sample) and `fastq/` (links to data/rna_fq). Participant-level.
