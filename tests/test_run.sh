#!/bin/bash

# first, make classifiers on you env.
# make-clfs

# next, check R environment
# cd R_script_for_test ; Rscript MutationalPatterns_from_single_vcf.R ; cd ..

# now, you can test run
gs-practice -i ./input/cancer_sample01.vcf -o ./output/cancer_sample01

gs-practice -i ./test_path_to_cancer_samples_10.list -o ./output/cancer_samples_10 -sn ./input/sample_names_10.list

gs-practice -i ./input/TCGA_mutect2_random100_1.maf -o ./output/TCGA_mutect2_random100_1
