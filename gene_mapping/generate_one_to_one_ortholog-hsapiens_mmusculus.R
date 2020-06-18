############### SYNOPSIS ###################
### AIM: Generate hsapiens to mmusculus unique orthologs gene mapping (via biomaRt)

### DESCRIPTION
# ...

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:
# Ensembl biomaRt guide: http://www.ensembl.org/info/data/biomart/biomart_r_package.html

### AUTHOR:
# Pascal N. Timshel

# ======================================================================= #
# ==============================  SETUP  ============================== #
# ======================================================================= #

library(biomaRt)
library(tidyverse)

# ======================================================================= #
# ==============================  biomaRt  ============================== #
# ======================================================================= #
### Set version variable
#ENSEMBL_VERSION = NULL # for newest version

### GRCh37.ens_v91
# ENSEMBL_VERSION = 91 # Ensembl release 91 = December 2017 | for specific version
# GENOME_VERSION = 37 # GRCh37=hg19; GRCh38=hg38

### GRCh38.ens_v100
ENSEMBL_VERSION = 100 # Ensembl release 91 = December 2017 | for specific version
GENOME_VERSION = 38 # GRCh37=hg19; GRCh38=hg38 # --> "Only 37 can be specified for GRCh version"


### Connect to the BioMart database and dataset hosted by Ensembl
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ENSEMBL_VERSION, GRCh=GENOME_VERSION, verbose=T)


# ======================================================================= #
# ========================= PAGE = GENES FEATURES  ====================== #
# ======================================================================= #

### Get all Ensembl Gene IDs - JUST TO SEE HOW MANY GENES THERE ARE IN TOTAL
df.BM.all.ids <- getBM(attributes = c("ensembl_gene_id"), mart=ensembl)
nrow(df.BM.all.ids) # 66203 for ensembl_v84
str(df.BM.all.ids)

### Features page
df.BM.feature <- getBM(attributes = c("ensembl_gene_id", 
                                  "chromosome_name",
                                  "start_position", # Gene start (bp)
                                  "end_position" # Gene end (bp)
                                  ), 
                          filters = c("with_mmusculus_homolog"), values=TRUE,
                          mart=ensembl)
str(df.BM.feature)
nrow(df.BM.feature) 

length(unique(df.BM.feature$ensembl_gene_id)) == nrow(df.BM.feature) # --> TRUE, all ensembl_gene_id are unique

# ======================================================================= #
# ============================= PAGE = HOMOLOGS ========================= #
# ======================================================================= #

### Homolog page
# We get all genes, but we could also have used a filter.
df.BM.all.homolog <- getBM(attributes = c("ensembl_gene_id",  # human gene ID
                                          "mmusculus_homolog_ensembl_gene", # mouse gene ID
                                          "mmusculus_homolog_orthology_confidence" # gives 0 or 1 scores that can help choose the optimal ortholog for NON-UNIQUE human-mouse orthologs.
                                          ), mart=ensembl)
str(df.BM.all.homolog)

### GRCh37.ens_v91
# df.BM.all.homolog %>% count(mmusculus_homolog_orthology_confidence)
# 0 = 4300
# 1 = 16602
# NA = 45453

### GRCh38.ens_v100
# 0  8747
# 1 17960
# NA 47461

# ======================================================================= #
# =================================  JOIN  ============================== #
# ======================================================================= #

### Merge df.BM.all ###
### OBS: we do LEFT JOIN because df.BM.feature have been filtered to only contain gene IDs with mouse homologs
df.BM <- df.BM.feature %>% left_join(df.BM.all.homolog, by="ensembl_gene_id")

# ======================================================================= #
# ======================== KEEP ONLY 1-1 ORTHOLOGS ====================== #
# ======================================================================= #

# df.BM.high_conf <- df.BM %>% filter(mmusculus_homolog_orthology_confidence==1)

human.non_unique_orthologs <- df.BM %>% group_by(ensembl_gene_id) %>% filter(n()>1) %>% pull(ensembl_gene_id)
mouse.non_unique_orthologs <- df.BM %>% group_by(mmusculus_homolog_ensembl_gene) %>% filter(n()>1) %>% pull(mmusculus_homolog_ensembl_gene)

df.BM.unique_orthologs <- df.BM %>% filter(!ensembl_gene_id %in% human.non_unique_orthologs, 
                                           !mmusculus_homolog_ensembl_gene %in% mouse.non_unique_orthologs)


### Validate that each mouse and human Ensembl ID is only listed ONCE! (Then we have a unique mapping)
length(unique(df.BM.unique_orthologs$ensembl_gene_id)) == nrow(df.BM.unique_orthologs) # --> TRUE, all unique
length(unique(df.BM.unique_orthologs$mmusculus_homolog_ensembl_gene)) == nrow(df.BM.unique_orthologs) # --> TRUE, all unique


# ======================================================================= #
# ===========================  SUMMARY STATS  =========================== #
# ======================================================================= #

str(df.BM.unique_orthologs) 
# GRCh37.ens_v91: N=16228
# GRCh38.ens_v100: N=16766

df.BM.unique_orthologs %>% count(mmusculus_homolog_orthology_confidence)
# 0 = 1268
# 1 = 14960

df.BM.unique_orthologs %>% count(chromosome_name) %>% print(n=50)
# ---> does contain X and Y chromosomes
# ....
# 22               9   637
# 23      GL000191.1     1
# 24      GL000213.1     1
# 25              MT    11
# 26               X   555
# 27               Y     7

# ======================================================================= #
# ===============================  EXPORT  ============================== #
# ======================================================================= #

file.ensmbl_annotation <- "hsapiens_mmusculus_unique_orthologs.GRCh38.ens_v100.txt.gz"
df.BM.unique_orthologs %>% write_tsv(file.ensmbl_annotation)

### Write table - OLD
## file.ensmbl_annotation <- sprintf("biomaRt-annotation-hsapiens_ensembl-v%s.txt.gz", ifelse(is.null(ENSEMBL_VERSION), "NEWEST_RELEASE_CHECK_R_SCRIPT", ENSEMBL_VERSION))
# file.ensmbl_annotation <- "hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
# file.ensmbl_annotation.gz <- gzfile(file.ensmbl_annotation, 'w')
# write.table(df.BM.unique_orthologs, file=file.ensmbl_annotation.gz, col.names=T, row.names=F, quote=F, sep="\t")
# close(file.ensmbl_annotation.gz)




