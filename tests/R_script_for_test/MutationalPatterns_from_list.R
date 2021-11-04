# R script for a list of multiple VCF files
input_data="./test_path_to_vcf_files.list"
genome_version="hg38"
sample_names=NULL # or
#sample_names=readLines("../../tests/input/sample_names_10.list")
cosmic="../../src/lib/cosmic_v2_signatures.tsv"
# import library
library(MutationalPatterns)
if (genome_version=="hg38"){ref_genome <-"BSgenome.Hsapiens.UCSC.hg38"}
if (genome_version=="hg19"){ref_genome <-"BSgenome.Hsapiens.UCSC.hg19"}
library(ref_genome, character.only = TRUE)
# input_data
vcfs=readLines(input_data)
# set sample_names
if (is.null(sample_names)){
  basenames = gsub(pattern = "^.*/",replacement="",vcfs)
  sample_names = gsub(".vcf","",basenames)}

# make grl
grl = read_vcfs_as_granges(vcf_files = vcfs,
                           sample_names = sample_names,
                           genome = ref_genome)
# 96 mutational profile
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
# COSMIC reference signature
cancer_signatures = read.table(cosmic,sep = "\t",header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$MutationType)
# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$MutationType
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:ncol(cancer_signatures)])
# Fit mutation matrix to the COSMIC mutational signatures:
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
out_df=t(fit_res$contribution)
out_df

# write out
# out_f="cancer_samples_10_decomposed.tsv"
#write.table(x=out_df,file=out_f,
#            sep="\t",quote = F,col.names = NA,row.names = T)
