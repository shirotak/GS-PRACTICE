#!/usr/bin/env python3
"""
Main CLI interface
This tool implements a pipeline to classify cancer patients into a specific
subgroup based on their mutational profile.
It does not handle raw NGS FASTQ data, but rather expects VCF/MAF input.
"""
import sys
import argparse
import logging
from gspractice import *

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
