use clap::{ArgAction, Parser};

const GROUP: &str = "Column names. NB: case insensitive.";

#[derive(Parser, Debug)]
#[command(
    name = "munge_sumstats",
    version = "0.1",
    author = "Wenjie Wei <weiwenjie@westlake.edu.cn>",
    about = "Munge summary statistics"
)]
struct Args {
    #[arg(long, default_value = None, help = "Input filename.")]
    sumstats: Option<String>,

    #[arg(long="N", default_value = None, help = "Sample size. If this option is not set, will try to infer the sample size from the input file. If the input file contains a sample size column, and this flag is set, the argument to this flag has priority.")]
    n: Option<f64>,

    #[arg(long="N-cas", default_value = None, help = "Number of cases. If this option is not set, will try to infer the number of cases from the input file. If the input file contains a number of cases column, and this flag is set, the argument to this flag has priority.")]
    n_cas: Option<f64>,

    #[arg(long="N-con", default_value = None, help = "Number of controls. If this option is not set, will try to infer the number of controls from the input file. If the input file contains a number of controls column, and this flag is set, the argument to this flag has priority.")]
    n_con: Option<f64>,

    #[arg(long, default_value = None, help = "Output filename prefix.")]
    out: Option<String>,

    #[arg(long, default_value_t = 0.9, help = "Minimum INFO score.")]
    info_min: f64,

    #[arg(long, default_value_t = 0.01, help = "Minimum MAF.")]
    maf_min: f64,

    #[arg(long, action = ArgAction::SetTrue, help = "Use this flag to parse Stephan Ripke's daner* file format.")]
    daner: bool,

    #[arg(long, action = ArgAction::SetTrue, help = "Use this flag to parse more recent daner* formatted files, which include sample size column 'Nca' and 'Nco'.")]
    daner_n: bool,

    #[arg(long, action = ArgAction::SetTrue, help = "Don't require alleles. Useful if only unsigned summary statistics are available and the goal is h2 / partitioned h2 estimation rather than rg estimation.")]
    no_alleles: bool,

    #[arg(long, default_value = None, help = "Same as --merge, except the file should have three columns: SNP, A1, A2, and all alleles will be matched to the --merge-alleles file alleles.")]
    merge_alleles: Option<String>,

    #[arg(long, default_value = None, help = "Minimum N (sample size). Default is (90th percentile N) / 2.")]
    n_min: Option<f64>,

    #[arg(long, default_value_t = 5e6 as i64, help = "Chunksize.")]
    chunksize: i64,

    #[arg(long, default_value = None, help = "Name of SNP column (if not a name that ldsc understands). ")]
    snp: Option<String>,

    #[arg(long="N-col", default_value = None, help = "Name of N column (if not a name that ldsc understands). ", help_heading=Some(GROUP))]
    n_col: Option<String>,

    #[arg(long="N-cas-col", default_value = None, help = "Name of N column (if not a name that ldsc understands). ", help_heading=Some(GROUP))]
    n_cas_col: Option<String>,

    #[arg(long="N-con-col", default_value = None, help = "Name of N column (if not a name that ldsc understands). ", help_heading=Some(GROUP))]
    n_con_col: Option<String>,

    #[arg(long, default_value = None, help = "Name of A1 column (if not a name that ldsc understands). ", help_heading=Some(GROUP))]
    a1: Option<String>,

    #[arg(long, default_value = None, help = "Name of A2 column (if not a name that ldsc understands). ", help_heading=Some(GROUP))]
    a2: Option<String>,

    #[arg(long, default_value = None, help = "Name of p-value column (if not a name that ldsc understands). ", help_heading=Some(GROUP))]
    p: Option<String>,

    #[arg(long, default_value = None, help = "Name of FRQ or MAF column (if not a name that ldsc understands). ", help_heading=Some(GROUP))]
    frq: Option<String>,

    #[arg(long, default_value = None, help = "Name of signed sumstat column, comma null value (e.g., Z,0 or OR,1). ", help_heading=Some(GROUP))]
    signed_sumstats: Option<String>,

    #[arg(long, default_value = None, help = "Name of INFO column (if not a name that ldsc understands). ", help_heading=Some(GROUP))]
    info: Option<String>,

    #[arg(long, default_value = None, help = "Comma-separated list of INFO columns. Will filter on the mean. ", help_heading=Some(GROUP))]
    info_list: Option<String>,

    #[arg(long, default_value = None, help = "Name of NSTUDY column (if not a name that ldsc understands). ", help_heading=Some(GROUP))]
    nstudy: Option<String>,

    #[arg(long, default_value = None, help = "Minimum # of studies. Default is to remove everything below the max, unless there is an N column, in which case do nothing.", help_heading=Some(GROUP))]
    nstudy_min: Option<f64>,

    #[arg(long, default_value = None, help = "Comma-separated list of column names to ignore.", help_heading=Some(GROUP))]
    ignore: Option<String>,

    #[arg(long, action = ArgAction::SetTrue, help = "A1 is the increasing allele.", help_heading=Some(GROUP))]
    a1_inc: bool,

    #[arg(long, action = ArgAction::SetTrue, help = "Keep the MAF column (if one exists).", help_heading=Some(GROUP))]
    keep_maf: bool,
}

fn main() {
    let args = Args::parse();

    // Access the values of the arguments here
    if let Some(sumstats) = args.sumstats {
        println!("Input filename: {}", sumstats);
    }
    if let Some(n) = args.n {
        println!("Sample size: {}", n);
    }
    // Add more argument handling as needed
}
