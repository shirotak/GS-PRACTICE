# GS-PRACTICE (Genomic Subtyping and Predictive Response Analysis for Cancer Tumor ICi Efficacy)
Tumor genomic subtyping tool based on mutational signatures for cancer samples.  
___
Overview of the pipeline  
<img src=https://github.com/shirotak/GS-PRACTICE/blob/main/documentation/Pipeline_overview.png width="500">  
___
## Requirements
### Python >= 3.7
- Please see requirements.txt
### R >= 4.0
- MutationalPatterns
- BSgenome.Hsapiens.UCSC.hg38  
- BSgenome.Hsapiens.UCSC.hg19  
___
## Installation
To avoid package dependency issues, I recommend installation in a virtual environment created with anaconda, miniconda, miniforge, pyenv etc.  
```
git clone shirotak/GS-PRACTICE.git
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
After install R enviroments,
```
pip install rpy2
```
___
## Preparations
#### Make classifiers in your environment
After installing the requirements, you should configure classifiers and umap projector in your environment and save them as joblib files.
```
cd src/gspractice
python makeclfs.py
```
If you want, you can change hyper parameters of classifiers by directly editing the script `src/gspractice/makeclfs.py`.  
By default, joblib files named KNN, SVC, RFC, LRC, and UMAP_projector are generated in `src/lib` directory.
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
Or,
```
python src/gspractice/run.py -i {input_file} -o {output_prefix} 
```

### Input file
[VCF file](https://en.wikipedia.org/wiki/Variant_Call_Format)(version >= 4.0) or [MAF file](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) after somatic mutation calling are accepted.  
VCF file needs to contain a header starting with "##fileformat=VCFv**".  
MAF file needs to contain columns with named "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2" and "Tumor_Sample_Barcode".  
The default genome version is GRCh38 (hg38), but genomeGRCh37 (hg19) is also accepted if you specify the option `-gv hg19`.
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
cd test
bash test_run.sh
```
This will generates the following files in the "tests/output" directory
-  [cancer_sample01_decomposed.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_sample01_decomposed.tsv)
-  [cancer_sample01_predict.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_sample01_prediction.txt)
-  [cancer_sample01_umap.png](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_sample01_umap.png)
<img src=https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_umap.png width="350">

-  [cancer_samples_10_decomposed.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_samples_10_decomposed.tsv)
-  [cancer_samples_10_table.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_samples_10_prediction.txt)
-  [cancer_samples_10_umap.png](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_samples_10_umap.png)
<img src=https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/cancer_samples_10_umap.png width="350">  

-  [TCGA_mutect2_random100_1_decomposed.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_decomposed.tsv)
-  [TCGA_mutect2_random100_1_table.tsv](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_prediction.txt)
-  [TCGA_mutect2_random100_1_umap.png](https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_umap.png)
<img src=https://github.com/shirotak/GS-PRACTICE/blob/main/tests/output/TCGA_mutect2_random100_1_umap.png width="350">  

___
## Getting Help  
```
gs-practice -h
```
___
## Citation
Currently, the software is being published as a preprint in MedRxiv.  
[Tumor genomic subtypes orthogonal to mutation burden predict the efficacy of immune checkpoint therapy](https://www.medrxiv.org/content/10.1101/2021.10.03.21264330v1)
