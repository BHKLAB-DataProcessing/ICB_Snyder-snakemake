from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]
data_source  = "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Snyder-data/main/"

rule get_MultiAssayExp:
    output:
        S3.remote(prefix + filename)
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        # S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "annotation/Gencode.v19.annotation.RData")
    resources:
        mem_mb=4000
    shell:
        """
        Rscript -e \
        '
        load(paste0("{prefix}", "annotation/Gencode.v19.annotation.RData"))
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "Snyder", input_dir = paste0("{prefix}", "processed")), 
            "{prefix}{filename}"
        );
        '
        """

# rule format_snv:
#     output:
#         S3.remote(prefix + "processed/SNV.csv")
#     input:
#         S3.remote(prefix + "download/SNV.txt"),
#         S3.remote(prefix + "processed/cased_sequenced.csv")
#     resources:
#         mem_mb=2000
#     shell:
#         """
#         Rscript scripts/Format_SNV.R \
#         {prefix}download \
#         {prefix}processed \
#         """

rule format_expr:
    output:
        S3.remote(prefix + "processed/EXPR.csv")
    input:
        S3.remote(prefix + "download/EXPR.txt"),
        S3.remote(prefix + "annotation/Gencode.v19.annotation.RData"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_EXPR.R \
        {prefix}download \
        {prefix}processed \
        {prefix}annotation
        """

rule format_clin:
    input:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    output:
        S3.remote(prefix + "processed/CLIN.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_CLIN.R \
        {prefix}download \
        {prefix}processed \
        {prefix}annotation
        """

rule format_cased_sequenced:
    input:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/EXPR.txt")
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_cased_sequenced.R \
        {prefix}download \
        {prefix}processed \
        """

rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v19.annotation.RData"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    shell:
        """
        wget https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/Gencode.v19.annotation.RData?raw=true -O {prefix}annotation/Gencode.v19.annotation.RData
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_drug.csv -O {prefix}annotation/curation_drug.csv
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_tissue.csv -O {prefix}annotation/curation_tissue.csv 
        """

rule format_downloaded_data:
    input:
        S3.remote(prefix + "download/data_clinical.csv"),
        S3.remote(prefix + "download/2850417_Neoantigen_RNA_bams.csv"),
        S3.remote(prefix + "download/kallisto.zip")
    output:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/EXPR.txt")
        # S3.remote(prefix + "download/SNV.txt")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/format_downloaded_data.R {prefix}download
        """ 

rule download_data:
    output:
        S3.remote(prefix + "download/data_clinical.csv"),
        S3.remote(prefix + "download/2850417_Neoantigen_RNA_bams.csv"),
        S3.remote(prefix + "download/kallisto.zip")
    shell:
        '''
        Rscript scripts/download_data.R {prefix}download
        '''