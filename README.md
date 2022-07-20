# GS-PRACTICE (Genomic Subtyping and Predictive Response Analysis for Cancer Tumor ICi Efficacy)
Tumor genomic subtyping tool based on mutational signatures for cancer samples.  
See [the corresponding paper](https://www.medrxiv.org/content/10.1101/2021.10.03.21264330v2) for details.
___
Overview of the pipeline  
<img src=https://github.com/shirotak/GS-PRACTICE/blob/main/documentation/Pipeline_overview.png width="500">  
___
## Requirements
### Python >= 3.7
- numpy
- numba
- scikit_learn
- joblib
- umap_learn
- pandas
- matplotlib
- rpy2   
  
Please see `requirements.txt` for detailed versions.
### R >= 3.6
- MutationalPatterns
- BSgenome.Hsapiens.UCSC.hg38  
- BSgenome.Hsapiens.UCSC.hg19

Please see `r_requirements.txt` for detailed versions.  
These have been tested on mac OSX and Linux.
___
## Installation
To avoid package dependency problems, I recommend installation in a virtual environment created with anaconda, miniconda, miniforge, pyenv etc.   For example,
```
conda create -n GSP python=3.7
conda activate GSP
```
After downloading this repository,
```
cd GS-PRACTICE
python setup.py install
```
Or, from [PyPI](https://pypi.org/project/GS-PRACTICE),
```
pip install GS-PRACTICE
```
Then,
```
pip install -r requirements.txt
```
After installing required R packages (see `r_requirements.txt`),
```
pip install rpy2
```
___
## Preparations
#### Make classifiers in your environment
After installing the requirements, you should configure classifiers and umap projector in your environment and save them as joblib files.  
By default, joblib files named KNN, SVC, RFC, LRC, and UMAP_projector are generated in `gspractice/data` directory.
```
gs-makeclfs
```
#### Check that R works
Before using rpy2, you should check the behavior of R in your environment.  
For example,
```
cd tests/R_script_for_test
Rscript MutationalPatterns_from_single_vcf.R
Rscript MutationalPatterns_from_list.R
Rscript MutationalPatterns_from_maf.R
```
___
## Usage
```
gs-practice -i {input_file} -o {output_prefix} 
```
### Input files
[VCF file](https://en.wikipedia.org/wiki/Variant_Call_Format)(version >= 4.0) or [MAF file](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) after somatic mutation calling are accepted.  
VCF file needs to contain a header starting with "##fileformat=VCFv**".  
MAF file needs to contain columns named "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2" and "Tumor_Sample_Barcode".  
The default genome version is GRCh38 (hg38), but GRCh37 (hg19) is also accepted if you specify the option `-gv hg19`.
- single VCF file (end with ***.vcf)
- list of paths to multiple VCF files (end with ***.list)
- MAF file (end with ***.maf)
### Example
```
gs-practice -i input.vcf -o output_prefix [-sn samplename ]
gs-practice -i vcffiles.list -o output_prefix [-snl samplenames.list]
gs-practice -i input.maf -o output_prefix
```
### Output files
```
{output_prefix}_decomposed.tsv
{output_prefix}_prediction.txt
{output_prefix}_umap.png
```
___
## Test
```
cd tests
bash test_run.sh
```
This will generates the following files in the "tests/output" directory
-  [cancer_sample01_decomposed.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_sample01_decomposed.tsv)
-  [cancer_sample01_prediction.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_sample01_prediction.tsv)
-  [cancer_sample01_umap.png](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_sample01_umap.png)
<img src=https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_umap.png width="350">

-  [cancer_samples_10_decomposed.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_samples_10_decomposed.tsv)
-  [cancer_samples_10_prediction.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_samples_10_prediction.tsv)
-  [cancer_samples_10_umap.png](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_samples_10_umap.png)
<img src=https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_samples_10_umap.png width="350">  

-  [TCGA_mutect2_random100_1_decomposed.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_decomposed.tsv)
-  [TCGA_mutect2_random100_1_prediction.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_prediction.tsv)
-  [TCGA_mutect2_random100_1_umap.png](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_umap.png)
<img src=https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_umap.png width="350">  

___
## Getting Help  
```
gs-practice -h
```
___
## Notes
- If you want, you can change hyper parameters of classifiers by directly editing the script "gspractice/makeclfs.py". 
- You may use the tool without installation if you satisfy the requirements.
```
python src/gspractice/makeclfs.py
python src/gspractice/run_gspractice.py -h
python src/gspractice/run_gspractice.py -i {input_file} -o {output_prefix}
```
___
## License
Copyright (c) 2021 Shiro Takamatsu (Kyoto University). See [LICENSE](https://github.com/shirotak/GS-PRACTICE/blob/main/LICENSE).  
This software is currently prohibited for commercial use. If you wish to use the Software for commercial or for-profit purposes, please contact [the author](tshiro@kuhp.kyoto-u.ac.jp).
___
## Citation
Currently, the corresonding paper is published as a preprint in MedRxiv.  
[Mutation burden-orthogonal tumor genomic subtypes delineate responses to immune checkpoint therapy](https://www.medrxiv.org/content/10.1101/2021.10.03.21264330)
