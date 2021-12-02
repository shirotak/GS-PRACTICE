#!/usr/bin/env python3
"""
Main CLI interface
This tool implements a pipeline to classify cancer patients into a specific
subgroup based on their mutational profile.
It does not handle raw NGS FASTQ data, but rather expects VCF/MAF input.
"""
import os
import sys
import argparse
import logging
logger = logging.getLogger(__name__)
#logger.setLevel(logging.INFO)
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import joblib
import umap

# specify library directory
script_dir=os.path.dirname(os.path.abspath(__file__))
lib_dir=script_dir+"/data/"


# ----- Begin code for this module. -----

_toolName = "GS-PRACTICE"
_defOutPrefix = "GS-PRACTICE"  # Default output prefix.
_defDecompositionSuffix = "_decomposed.tsv"  # from MutationalPatterns
_defPredictionSuffix = "_prediction.tsv"  # from subtype prediction
_defPlotSuffix = "_umap"  # from umap projection with TCGA
_defFigureFormat="png"
_defFigureSize="5,5"

############################################################
def getRequiredInputsParser():
    """Constructs parser for absolute minimum required inputs."""

    parser = argparse.ArgumentParser(add_help=False)
    inputOpts = parser.add_argument_group("Input options")
    inputOpts.add_argument(
        '-i', "--input_data", type=str, metavar='.vcf, .list, .maf', required=True,
        help='Input filename for single VCF file, single list of paths to multiple vcf files, single MAF file')
    return parser
############################################################

############################################################
def getOptionalInputsParser():
    """Constructs parser for optional arguments."""

    parser = argparse.ArgumentParser(add_help=False)
    extraInOpts = parser.add_argument_group('Extra input options')
    extraInOpts.add_argument(
        "-sn", "--sample_name", type=str,metavar='string', required=False,
        help="""a sample name for single VCF file. This option is ignored for multiple VCFs or MAF file""")
    extraInOpts.add_argument(
        "-snl", "--sample_name_list", type=str,metavar='path to sample name list', required=False,
        help="""path to a list of sample names for multiple VCF files. This is ignored for single VCF or MAF file""")
    extraInOpts.add_argument(
        "-gv", '--genome_version', type=str, metavar='string',
        required=False, choices=["hg19", "hg38"], default="hg38",
        help='Reference genome version, choose from hg19 or hg38, default:hg38')
    extraInOpts.add_argument(
        "-udc","--use_default_clfs", action="store_true",required=False, #metavar='flag',
        help="""Use the default classifiers. It is recommended you make and save
        classifiers in your own environment. Default:False """)
    extraInOpts.add_argument(
        "-udu","--use_default_umap", action="store_true",required=False, #metavar='flag',
        help="""Use the default umap projector. It is recommended you make and save
        umap projector in your own environment. Default:False""")
    extraInOpts.add_argument(
        "-if","--input_format", choices=["maf","vcf","vcfs"],required=False,default=None, 
        help="""You can specify the format of the input file by selecting from 'maf', 'vcf', 'vcfs'.
        If not specified (default), will be inferred from extension name.""")

    return parser
############################################################

############################################################
def getOptionalOutputsParser():
    """Constructs parser for optional arguments."""

    parser = argparse.ArgumentParser(add_help=False)
    outOpts = parser.add_argument_group('Output options')
    outOpts.add_argument(
        "-o", '--out_prefix', type=str, metavar='prefix',
        required=False, default=_defOutPrefix,
        help="""Prefix (Path) to generating output files from analyses.
        Default: '%s'.""" % _defOutPrefix)
    outOpts.add_argument(
        "-nts",'--no_table_sort', action="store_true", #metavar='flag',
        help="""Do not sort the prediction result table.
        Default behavior is false (to sort).""")
    outOpts.add_argument(
        "-np", '--no_plot', action="store_true", #metavar='flag',
        help="""Do not plot mapping figure.
        Default behavior is false (to plot).""")
    outOpts.add_argument(
        "-ff", '--figure_format', type=str, metavar='string',
        required=False, default=_defFigureFormat, 
        dest='figure_format',
        help="""Formats supported by matplotlib are accepted,
        such as 'pdf','png','jpg','tif'... dafault is 'png' """)
    outOpts.add_argument(
        "-fs", '--figure_size', type=str, metavar='string',
        required=False, default=_defFigureSize, 
        dest='figure_size',
        help=""" Size of mapping figure. In two dimentional numbers separated by comma.
        Default: 5,5 """)
    outOpts.add_argument(
        "-ms", '--marker_size', type=float, metavar='float',
        required=False, default=None, 
        dest='marker_size',
        help=""" Size of markers in the figure. This is automatically set by the sample size. Default: None
         """)  

    logOptions = parser.add_argument_group('Logging options')
    logOptions.add_argument(
        '--logfile', type=str, metavar='filename',
        required=False, default=None,
        help="Output messages to given logfile, default is stderr.")
    logOptions.add_argument(
        "-v", "--verbose", action="store_true", 
        help="Increase output verbosity")

    return parser
############################################################

############################################################
def analyzeSignature(input_data, genome_version="hg38", sample_name=None,
                        sample_name_list=None,input_format=None):
    # import data and library
    cosmic=lib_dir+"cosmic_v2_signatures.tsv"
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        robjects.r("genome_version = '%s'" % genome_version)
        robjects.r("""
        suppressMessages(library(MutationalPatterns))
        if (genome_version=="hg38"){ref_genome <-"BSgenome.Hsapiens.UCSC.hg38"}
        if (genome_version=="hg19"){ref_genome <-"BSgenome.Hsapiens.UCSC.hg19"}
        suppressMessages( library(ref_genome, character.only = TRUE) ) 
        """)
        robjects.r("cosmic='%s'" % cosmic)

    # check format of input_data, make GenomeRange
    if (input_format=="maf") | (".maf" in input_data):
        logger.info("input file format = MAF file")
        robjects.r("in_f='%s'" % input_data)
        robjects.r("""
        suppressMessages(library(GenomicRanges))
        maf=read.delim(in_f,comment.char="#")
        # remove indels
        snp=maf$Variant_Type=="SNP"
        maf=maf[snp,]
        # adjust status
        if (!grepl(pattern="chr", x=maf$Chromosome[1])){
        maf$Chromosome=paste0("chr",maf$Chromosome)
        }
        colnames(maf)[which( colnames(maf)=="Reference_Allele" )]="REF"
        colnames(maf)[which( colnames(maf)=="Tumor_Seq_Allele2" )]="ALT"
        # make grl
        grl=makeGRangesListFromDataFrame(
        maf,
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
        # limit chromosomes
        chromosomes = paste0('chr', c(1:22,'X'))
        seqlevels(grl, pruning.mode = 'tidy') = chromosomes
        # set genome version 
        GenomeInfoDb::genome(grl)=genome_version
        """)
    elif (input_format=="vcf") | (".vcf" in input_data):
        logger.info("input file format = single VCF file")
        robjects.r("in_f='%s'" % input_data)
        if sample_name:
            robjects.r("sample_names='%s'" % sample_name)
        else:
            basename=os.path.basename( input_data )
            sample_name=basename.replace(".vcf","")
            robjects.r("sample_names='%s'" % sample_name)
        robjects.r("""
        # make grl
        grl = read_vcfs_as_granges(
        vcf_files = in_f,
        sample_names = sample_names,
        genome = ref_genome)
        """)
    elif (input_format=="vcfs") | (".list" in input_data):
        logger.info("input file format = a list of path to multiple VCF files")
        robjects.r("vcfs=readLines('%s')" % input_data )
        if sample_name_list:
            robjects.r("sample_names='%s'" % sample_name_list)
        else:
            robjects.r("sample_names=NULL")
        robjects.r("""
        if (is.null(sample_names)){
            basenames = gsub(pattern = "^.*/",replacement="",vcfs)
            sample_names = gsub(".vcf","",basenames)} else {
            sample_names = readLines(sample_names)}
        # make grl
        grl = read_vcfs_as_granges(
        vcf_files = vcfs,
        sample_names = sample_names,
        genome = ref_genome)
        """)
    else:
        logger.critical("Check input file name or format. Youcan select format by --input-format option.")
        sys.exit(-1)

    # MutationalPatterns
    robjects.r("""
    # 96 mutational profile
    mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
    # COSMIC reference signature
    cancer_signatures = read.table(cosmic, sep = "\t", header = TRUE)
    # Match the order of the mutation types to MutationalPatterns standard
    new_order = match(row.names(mut_mat), cancer_signatures$MutationType)
    # Reorder cancer signatures dataframe
    cancer_signatures = cancer_signatures[as.vector(new_order),]
    # Add trinucletiode changes names as row.names
    row.names(cancer_signatures) = cancer_signatures$MutationType
    # Keep only 96 contributions of the signatures in matrix
    cancer_signatures = as.matrix(cancer_signatures[,4:ncol(cancer_signatures)])
    # Fit mutation matrix to the COSMIC mutational signatures:
    fit_res = fit_to_signatures(mut_mat, cancer_signatures)
    out_df=t(fit_res$contribution)
    columns=colnames(out_df)
    index=rownames(out_df)
    """) 

    # back to python, make pandas dataframe
    df=pd.DataFrame(robjects.globalenv["out_df"],columns=robjects.globalenv["columns"],index=robjects.globalenv["index"]) 
    logger.info("finished analyzing mutational signature.")

    return df

############################################################

############################################################
def predictSubtype(df_decomp,no_table_sort=False, use_default_clfs=False):
    # import four classifiers
    if use_default_clfs:
        logger.warning("You are using default-version classifiers, we recommend you to make and save them in your environment.")
        svc=joblib.load(lib_dir+"TCGA_7181_svc_c1_g01.joblib")
        rfc=joblib.load(lib_dir+"TCGA_7181_rfc_ne100.joblib")
        knn=joblib.load(lib_dir+"TCGA_7181_knn_n5_dis.joblib")
        lrc=joblib.load(lib_dir+"TCGA_7181_lrc_c1.joblib")
    else:
        logger.info("Import four classifiers.")
        svc=joblib.load(lib_dir+"SVC.joblib")
        rfc=joblib.load(lib_dir+"RFC.joblib")
        knn=joblib.load(lib_dir+"KNN.joblib")
        lrc=joblib.load(lib_dir+"LRC.joblib")
    # import annotations & colors
    dict_cluster_numbers={"SMK":0, "UVL":1, "APB":2, "POL":3, "MRD":4,"HRD":5,
                        "GNS":6,"AGE":7,"UND":8}
    dict_cluster_points={"SMK":1, "UVL":1, "APB":1, "POL":1, "MRD":1,"HRD":0,"GNS":0,"AGE":0}
    
    # prediction            
    input_df=np.log10(df_decomp+1)
    preds=[]
    logger.info("Predicting subtype by the four classifiers")
    for clf in [knn,svc, rfc,lrc]:
        preds.append(clf.predict(input_df))
    df_preds=pd.DataFrame( {"KNN":preds[0],"SVC":preds[1],
    "RF":preds[2], "LR":preds[3]},
                        index=input_df.index)
    details,cons,topcounts=[],[],[]
    # consensus
    for idx in df_preds.index:
        counter=Counter( df_preds.loc[idx,"KNN":"LR"]).most_common()
        topcounts.append( counter[0][1] )
        detail=""
        for x in counter:
            detail+=x[0]+":"+str(x[1])+","
        details.append( detail[:-1] )
        if counter[0][1]>=3:
            cons.append(counter[0][0] )
        else:
            cons.append("UND")
    df_preds["TGS"]=cons
    df_preds["Details"]=details
    df_point=df_preds.iloc[:,0:4].replace( dict_cluster_points )
    df_preds["irGS_count"]=np.sum(df_point,axis=1)
    df_preds["irGS"]=(df_preds["irGS_count"]>=3).astype(int)

    # return
    if no_table_sort:
        logger.warning("Resulting table are not sorted by predicted subtypes.")
        return df_preds
    else:
        logger.info("Resulting table are sorted by predicted subtypes.")
        df_numbers=df_preds.iloc[:,0:5].replace(dict_cluster_numbers)
        df_numbers["TopCounts"]=topcounts
        df_numbers_sort=df_numbers.sort_values( ["TGS","TopCounts","KNN","SVC","RF","LR"],
                                        ascending=[True,False,True,True,True,True])
        df_numbers_sort=df_numbers_sort.astype(float)
        df_preds_sort=df_preds.loc[df_numbers_sort.index,:]

        return df_preds_sort 

############################################################

############################################################
def visualizePrediction(df_decomp,df_preds,out_f,figure_size="5,5",marker_size=None,use_default_umap=False):
    # import TCGA mapping
    if use_default_umap:
        logger.warning("You use default UMAP projection.")
        u=joblib.load(lib_dir+"TCGA_7181_umap_projector.joblib")
    else:
        logger.debug("Load UMAP configuration.")
        u=joblib.load(lib_dir+"UMAP_projector.joblib")
    tcga_mutsig_log10=pd.read_csv(lib_dir+"TCGA_7181_df_mutsig_log10.tsv",sep="\t",index_col=0)
    X0=u.transform(tcga_mutsig_log10)
    with open(lib_dir+"TCGA_7181_cluster_colors.tsv") as f:
        tcga_colors=f.read().splitlines()
    dict_cluster_colors={'SMK': 'red','UVL': 'blue', 'APB': 'green','POL': 'brown', 'MRD': 'purple',
                            'HRD': 'hotpink', 'GNS': 'c', 'AGE': 'y', 'UND': 'grey'}
    ## random shuffle before plotting
    shuffle_numbers=np.arange(len(X0))
    np.random.seed(777)
    np.random.shuffle(shuffle_numbers)
    X0=X0[shuffle_numbers,:]
    tcga_colors=[tcga_colors[x] for x in shuffle_numbers]

    # configure results
    input_df=np.log10(df_decomp+1)
    X=u.transform(input_df)
    sample_order=input_df.index
    df_preds_order=df_preds.loc[sample_order,:]
    colors=[ dict_cluster_colors[cluster] for cluster in df_preds_order["TGS"]]
    ## set marker size
    markeredgewidth = 1
    if marker_size is None:
        if input_df.shape[0] <= 30:
            marker_size=150
        elif input_df.shape[0] <=100:
            marker_size=100
            markeredgewidth=0.75
        elif input_df.shape[0] <=250:
            marker_size=50
            markeredgewidth=0.5
        else:
            marker_size=25
            markeredgewidth=0.2
    else:
        logger.warning("You set marker size = %d and figure size =%s " %(marker_size, figure_size))
    
    # plot and save
    figsize=tuple(float(x) for x in figure_size.split(","))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.scatter(X0[:,0],X0[:,1],c=tcga_colors,s=1)
    ax.scatter(X[:,0],X[:,1],
            edgecolors="k", c=colors,s=marker_size,marker="*",linewidth=markeredgewidth)
    ax.set_xlabel("UMAP1",fontsize=10)
    ax.set_ylabel("UMAP2",fontsize=10)
    ax.tick_params(labelsize=8)
    ax.set_xticks([-5,0,5,10,15])
    ax.set_yticks([-4,0,4,8,12,16]) 
    plt.rcParams["figure.dpi"]=300
    plt.savefig(out_f,dpi=300,bbox_inches="tight")

############################################################




############################################################
def getStandaloneParser():
    """Constructs parser to run this script as standalone application."""

    parser = argparse.ArgumentParser(
        description="""Perform genome variant decomposition by 
            MutationalPatterns package, subsequently classify the sample by
            predictions from models built by scikit-learn, and
            plot samples by umap algorithm projection with TCGA samples (optional)""",
        #epilog="""Results of processing and prediction are output in
        #    multiple files in the local directory.""",
        parents=(
            getRequiredInputsParser(),
            getOptionalInputsParser(),
            getOptionalOutputsParser()
            ))

    return parser
############################################################

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runAsStandalone():
    """Contains script execution flow when run as standalone application.
    """
    
    #### Pre-execution internal configuation check.
    #initLogger = logging.Logger(name=_toolName)
    #initLogger.addHandler(logging.StreamHandler(sys.stderr))
    
    #### Configuration fine, execute as planned.
    useParser = getStandaloneParser()
    args = useParser.parse_args()
    
    # Setup main logger for messages.
    logStream = logging.StreamHandler(
        open(args.logfile, "w") if args.logfile else sys.stderr)
    # Either streams of a logger or the logger itself can have its level set,
    #   but ONLY the logger itself can have the logger level retrieved.
    #logStream.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    logStream.setFormatter(logging.Formatter(
        "%(asctime)-15s %(levelname)s:%(message)s",
        datefmt="%Y-%m-%S %H:%M"))
    logger = logging.Logger(name=_toolName)
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    logger.addHandler(logStream)
    logger.debug("Invocation arguments: %s" % args)

    logging.getLogger("GS-PRACTICE").addHandler(logStream)
    logging.getLogger("GS-PRACTICE").setLevel(logging.INFO)
    
    # Start to analyze mutational signatures.
    logger.info("Start to analyze mutational signatures.")
    df_decomp=analyzeSignature(input_data=args.input_data,
                    genome_version=args.genome_version,
                    sample_name=args.sample_name,
                    sample_name_list=args.sample_name_list,
                    input_format=args.input_format)
    # Write results
    logger.debug("Decomposition result to %s " % args.out_prefix + _defDecompositionSuffix)
    df_decomp.to_csv(args.out_prefix + _defDecompositionSuffix,sep="\t")

    # Proceed to prediction by multi-estimator.
    logger.info("Start to predict subtypes.")
    df_preds=predictSubtype(df_decomp=df_decomp,
                    no_table_sort=args.no_table_sort,
                    use_default_clfs=args.use_default_clfs) 
    # Write out results
    logger.debug("Prediction result to %s" % args.out_prefix + _defPredictionSuffix)
    df_preds.to_csv(args.out_prefix + _defPredictionSuffix,sep="\t")
   
    # Proceed to plot umap
    if args.no_plot==False:
        out_f=args.out_prefix + _defPlotSuffix + "." + args.figure_format
        logger.debug("Plot and save mapping figure to %s" % out_f)
        visualizePrediction(df_decomp=df_decomp,
                df_preds=df_preds,
                out_f=out_f,
                figure_size=args.figure_size,
                marker_size=args.marker_size,
                use_default_umap=args.use_default_umap)
    else:
        logger.warning("Skip plotting figure.")
        pass

    # Clean up and terminate.
    logger.info("Execution completed.")

# end of running script as standalone application.

# What to do if script is a main driver program.
if __name__ == "__main__":

    import doctest
    numFail, numTests = doctest.testmod(verbose=False)
    if numFail:
        logger.critical('Expected functionality in doctests fail. '
              + 'Aborting standalone execution.')
        sys.exit(-1)
    # Otherwise, doctests passed, run standalone application.
    runAsStandalone()

# ----- End code for this module. -----
