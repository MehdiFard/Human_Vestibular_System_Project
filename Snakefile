import pandas as pd

# # Load fastqs list
# fq_list = pd.read_csv("src/path_to_fqs.tsv", sep="\t")
# # Create list of (library_id, human_id) tuples
# bam_info = list(fq_list[['library_id', 'human_id']].itertuples(index=False, name=None))


# Load fastqs list
fq_list = pd.read_csv("src/path_to_fqs.tsv", sep="\t")

# Define libraries to exclude
excluded_libs = ["FTL_500", "FTL_501"]

# Filter out excluded libraries
fq_list = fq_list[~fq_list['library_id'].isin(excluded_libs)]

# Create list of (library_id, human_id) tuples
bam_info = list(fq_list[['library_id', 'human_id']].itertuples(index=False, name=None))


rule all:
    input:
        "data/rna_fq/.symlinks_done",
        expand("data/bam/{lib_id}.Aligned.sortedByCoord.out.bam", lib_id=[lib for lib, hum in bam_info]),
        "data/bam/FTL_500_501.Aligned.sortedByCoord.out.bam"


rule symlink_fastq:
    """Create symbolic links to all files and folders in the top level of the RNA-seq FASTQ directory."""
    output:
        "data/rna_fq/.symlinks_done"
    params:
        src="/data/hps/assoc/private/gdml/user/sesmae/projects_data/vestibular_project_data/human/rnaseq_fastq",
        dest="data/rna_fq/"
    shell:
        r"""
        mkdir -p {params.dest}
        ln -sf {params.src}/* {params.dest}
        touch {output}
        """

# rule make_metadata:
#     """Create metadata file for RNA-seq FASTQ files."""


rule map_star:
    input:
        R1 = lambda wildcards: fq_list.loc[fq_list['library_id'] == wildcards.lib_id, 'R1'].values[0],
        R2 = lambda wildcards: fq_list.loc[fq_list['library_id'] == wildcards.lib_id, 'R2'].values[0]
    output:
        bam = "data/bam/{lib_id}.Aligned.sortedByCoord.out.bam"
    params:
        star = "/data/hps/assoc/private/gdml/user/sesmae/tools/STAR-2.7.11b/source/STAR",
        out_prefix = "data/bam/{lib_id}.",
        star_index = "/data/hps/assoc/private/gdml/user/sesmae/refs_gtfs_indx_vcfs/GRCh38.p14_ensembl/index_star",
        sample_id = lambda wildcards: fq_list.loc[fq_list['library_id'] == wildcards.lib_id, 'human_id'].values[0]
    threads: 16
    log:
        "logs/map_star/{lib_id}.log"
    shell:
        """
        # Making read group info
        IDENTIFIRE=$(head -n 1 <(zcat {input.R1}))
        INSTRO_NAME=$(echo $IDENTIFIRE | awk -F ':' '{{print $1}}' | cut -f 2 -d " " | sed 's/@//')
        RUN_ID=$(echo $IDENTIFIRE | awk -F ':' '{{print $2}}')
        FLOWCELL=$(echo $IDENTIFIRE | awk -F ':' '{{print $3}}')
        LANE=$(echo $IDENTIFIRE | awk -F ':' '{{print $4}}')
        TILE=$(echo $IDENTIFIRE | awk -F ':' '{{print $5}}')

        RGID=$FLOWCELL.$LANE
        RGSM={params.sample_id}
        RGLB={wildcards.lib_id}
        RGPU=$FLOWCELL.$LANE.$TILE
        RGPL="ILLUMINA"

        mkdir -p data/bam/
        mkdir -p logs/map_star/

        # Run STAR alignment
        {params.star}   --runMode alignReads \
                        --runThreadN {threads} \
                        --genomeDir {params.star_index} \
                        --readFilesIn {input.R1} {input.R2} \
                        --readFilesCommand zcat \
                        --twopassMode Basic \
                        --quantMode TranscriptomeSAM GeneCounts \
                        --outSAMstrandField intronMotif \
                        --outSAMattrRGline ID:$RGID SM:$RGSM LB:$RGLB PL:$RGPL PU:$RGPU \
                        --outSAMtype BAM SortedByCoordinate \
                        --outSAMunmapped Within \
                        --outFileNamePrefix {params.out_prefix} \
                        &> {log}
        """


rule map_star_FTL_500_FTL_501:
    input:
        RUN1_R1="data/rna_fq/batch1/FTL_500_S9_L003_R1_001.fastq.gz",
        RUN1_R2="data/rna_fq/batch1/FTL_500_S9_L003_R2_001.fastq.gz",
        RUN2_R1="data/rna_fq/batch1/FTL_501_S10_L003_R1_001.fastq.gz",
        RUN2_R2="data/rna_fq/batch1/FTL_501_S10_L003_R2_001.fastq.gz"
    output:
        bam = "data/bam/FTL_500_501.Aligned.sortedByCoord.out.bam"
    params:
        star = "/data/hps/assoc/private/gdml/user/sesmae/tools/STAR-2.7.11b/source/STAR",
        out_prefix = "data/bam/FTL_500_501.",
        star_index = "/data/hps/assoc/private/gdml/user/sesmae/refs_gtfs_indx_vcfs/GRCh38.p14_ensembl/index_star",
        sample_id = "31661958"
    threads: 16
    log:
        "logs/map_star/FTL_500_FTL_501.log"
    shell:
        r"""
        # Making read group info
        IDENTIFIRE=$(head -n 1 <(zcat {input.RUN1_R1}))
        INSTRO_NAME=$(echo $IDENTIFIRE | awk -F ':' '{{print $1}}' | cut -f 2 -d " " | sed 's/@//')
        RUN_ID=$(echo $IDENTIFIRE | awk -F ':' '{{print $2}}')
        FLOWCELL=$(echo $IDENTIFIRE | awk -F ':' '{{print $3}}')
        LANE=$(echo $IDENTIFIRE | awk -F ':' '{{print $4}}')
        TILE=$(echo $IDENTIFIRE | awk -F ':' '{{print $5}}')

        RGID=$FLOWCELL.$LANE
        RGSM={params.sample_id}
        RGLB="FTL_500_FTL_501"
        RGPU=$FLOWCELL.$LANE.$TILE
        RGPL="ILLUMINA"

        # Run STAR alignment
        {params.star}   --runMode alignReads \
                        --runThreadN {threads} \
                        --genomeDir {params.star_index} \
                        --readFilesIn {input.RUN1_R1},{input.RUN2_R1} {input.RUN1_R2},{input.RUN2_R2} \
                        --readFilesCommand zcat \
                        --twopassMode Basic \
                        --quantMode TranscriptomeSAM GeneCounts \
                        --outSAMstrandField intronMotif \
                        --outSAMattrRGline ID:$RGID SM:$RGSM LB:$RGLB PL:$RGPL PU:$RGPU \
                        --outSAMtype BAM SortedByCoordinate \
                        --outSAMunmapped Within \
                        --outFileNamePrefix {params.out_prefix} \
                        &> {log}
        """

