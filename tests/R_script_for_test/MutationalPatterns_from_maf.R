# R script for MAF file
input_data="../input/TCGA_mutect2_random100_1.maf"
genome_version="hg38"
#sample_names=NULL
cosmic="../../src/lib/cosmic_v2_signatures.tsv"
# R script for MAF file
input_data="../input/TCGA_mutect2_random100_1.maf"
genome_version="hg38"
#sample_names=NULL
cosmic="../../src/lib/cosmic_v2_signatures.tsv"
# import library
library(MutationalPatterns)
if (genome_version=="hg38"){ref_genome <-"BSgenome.Hsapiens.UCSC.hg38"}
if (genome_version=="hg19"){ref_genome <-"BSgenome.Hsapiens.UCSC.hg19"}
library(ref_genome, character.only = TRUE)

# From MAF file make grande object
library(GenomicRanges)
maf=read.delim(input_data,comment.char="#")
# remove indels
snp=maf$Variant_Type=="SNP"
maf=maf[snp,]
# adjust status
if ( !grepl(pattern = "chr",x =maf$Chromosome[1])){
  maf$Chromosome=paste0("chr",maf$Chromosome)
}
colnames(maf)[which( colnames(maf)=="Reference_Allele" )]="REF"
colnames(maf)[which( colnames(maf)=="Tumor_Seq_Allele2" )]="ALT"
# make grl
grl=makeGRangesListFromDataFrame(maf,
  keep.extra.columns=TRUE,
  ignore.strand=TRUE,
  seqinfo=NULL,
  seqnames.field=c("Chromosome"),
  start.field="Start_Position",
  end.field="End_Position",
  strand.field="Strand",
  starts.in.df.are.0based=FALSE,
  split.field = "Tumor_Sample_Barcode",
  names.field = "Hugo_Symbol")

chromosomes = paste0('chr', c(1:22,'X'))
seqlevels(grl, pruning.mode = 'tidy') = chromosomes

GenomeInfoDb::genome(grl) = genome_version

# 96 mutational profile
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
# COSMIC mutational signatures
cancer_signatures = read.table( cosmic,sep = "\t",header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$MutationType)
# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$MutationType
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
#Fit mutation matrix to the COSMIC mutational signatures:
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
out_df=t(fit_res$contribution)
out_df

# write out
# out_f="TCGA_mutect2_random100_1_decomposed.tsv"
#write.table(x=out_df,file=out_f,
#            sep="\t",quote = F,col.names = NA,row.names = T)
# import library
library(MutationalPatterns)
if (genome_version=="hg38"){ref_genome <-"BSgenome.Hsapiens.UCSC.hg38"}
if (genome_version=="hg19"){ref_genome <-"BSgenome.Hsapiens.UCSC.hg19"}
library(ref_genome, character.only = TRUE)

# From MAF file make grande object
library(GenomicRanges)
maf=read.delim(input_data,comment.char="#")
# remove indels
snp=maf$Variant_Type=="SNP"
maf=maf[snp,]
# adjust status
if ( !grepl(pattern = "chr",x =maf$Chromosome[1])){
  maf$Chromosome=paste0("chr",maf$Chromosome)
}
colnames(maf)[which( colnames(maf)=="Reference_Allele" )]="REF"
colnames(maf)[which( colnames(maf)=="Tumor_Seq_Allele2" )]="ALT"
# make grl
grl=makeGRangesListFromDataFrame(maf,
  keep.extra.columns=TRUE,
  ignore.strand=TRUE,
  seqinfo=NULL,
  seqnames.field=c("Chromosome"),
  start.field="Start_Position",
  end.field="End_Position",
  strand.field="Strand",
  starts.in.df.are.0based=FALSE,
  split.field = "Tumor_Sample_Barcode",
  names.field = "Hugo_Symbol")

chromosomes = paste0('chr', c(1:22,'X'))
seqlevels(grl, pruning.mode = 'tidy') = chromosomes

GenomeInfoDb::genome(grl) = genome_version

# 96 mutational profile
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
# COSMIC mutational signatures
cancer_signatures = read.table( cosmic,sep = "\t",header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$MutationType)
# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$MutationType
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
#Fit mutation matrix to the COSMIC mutational signatures:
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
out_df=t(fit_res$contribution)
out_df

# write out
# out_f="TCGA_mutect2_random100_1_decomposed.tsv"
#write.table(x=out_df,file=out_f,
#            sep="\t",quote = F,col.names = NA,row.names = T)
