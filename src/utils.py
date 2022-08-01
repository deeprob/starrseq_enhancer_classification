import os
import json
from argparse import Namespace
import subprocess
import multiprocessing as mp
import pandas as pd
import pybedtools

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath("__file__"))


###############################
# read meta file; create args #
###############################

def create_args(meta_file, lib_name):
    with open(meta_file, "r") as f: 
        meta_dict = json.load(f)
        
    args = Namespace(
        # from metadata file
        library_prefix = meta_dict[lib_name]["prefix"],
        library_reps = meta_dict[lib_name]["replicates"],
        library_pair= meta_dict[lib_name]["read_pairs"],
        library_umi = meta_dict[lib_name]["umi"],
        library_suffix = meta_dict[lib_name]["suffix"],
        library_short = meta_dict[lib_name]["shortform"],
        reference_genome = meta_dict["genome"]["ref_fasta"],
        reference_genome_twobit = meta_dict["genome"]["ref_twobit"],
        roi_file = meta_dict["roi"]["filtered"]
    )

    return args


###################
# filepath parser #
###################

def get_lib_peak_filepath(peak_dir, lib_short, peak_caller="starrpeaker"):
    peak_filename_dict = {"starrpeaker": "peaks.peak.final.bed"}
    peak_filepath = os.path.join(
        peak_dir, lib_short, peak_caller, peak_filename_dict[peak_caller]
        )
    return peak_filepath

def get_lib_diff_activity_peak_filepath(store_dir, lib_short, diff_activity_type):
    peak_filepath = os.path.join(
        store_dir, lib_short, f"{diff_activity_type}.bed"
        )
    return peak_filepath

def get_rep_bam_files(bam_dir, lib_short, lib_prefix, lib_reps):
    bam_files = [os.path.join(bam_dir, lib_short, f"{'_'.join([lib_prefix, rep])}.bam") for rep in lib_reps.split()]
    return bam_files


#########################
# diff activity helpers #
#########################

def save_diff_activity_peak_file(diff_activity_peak_bed, store_dir, lib_short, diff_activity_type):
    diff_activity_peak_file = get_lib_diff_activity_peak_filepath(store_dir, lib_short, diff_activity_type)
    os.makedirs(os.path.dirname(diff_activity_peak_file), exist_ok=True)
    diff_activity_peak_bed.moveto(diff_activity_peak_file)
    return

def get_replicate_norm_cov(peak_bed, bam_file):
    cov_bed = peak_bed.coverage(bam_file)
    cov_df = cov_bed.to_dataframe(disable_auto_names=True, header=None)
    cov_df = cov_df.set_index([0,1,2])
    cov_df.index = cov_df.index.rename(["chrom", "start", "end"])
    cov_reads = cov_df.iloc[:, -4]
    return cov_reads

def get_replicate_wise_cov_df(peak_bed, bam_files, lib_name):
    pool_iter = [(peak_bed, bf) for bf in bam_files]
    rep_cov_sers = [get_replicate_norm_cov(*pi) for pi in pool_iter] # multi_args_pool_job(get_replicate_norm_cov, pool_iter)
    rep_cov_dfs = pd.concat(rep_cov_sers, axis=1)
    rep_cov_dfs.columns = [f"{lib_name}_R{i}" for i in range(1, rep_cov_dfs.shape[1] + 1)]
    return rep_cov_dfs

def save_deseq_in_df(peak_cov_df, store_dir, lib_short):
    save_file = os.path.join(store_dir, lib_short, "deseq_in.csv")
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    # convert the chromosomal coordinates to proper format
    peak_cov_df = peak_cov_df.reset_index().drop_duplicates(keep="first")
    peak_cov_df["unique_id"] = peak_cov_df.chrom + "_" + peak_cov_df.start.astype("str") + "_" + peak_cov_df.end.astype("str")
    peak_cov_df = peak_cov_df.drop(columns= ["chrom", "start", "end"]).set_index("unique_id")
    peak_cov_df.to_csv(save_file, index=True)    
    return save_file

def get_deseq_out_file(store_dir, lib_short):
    save_file = os.path.join(store_dir, lib_short, "deseq_out.csv")
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    return save_file

def da_with_deseq(table_in, table_out):
    cmd = [
        "bash", f"{CURRENT_DIR_PATH}/scripts/run_DESeq2.sh", 
        table_in, table_out
        ]
    subprocess.run(cmd)
    return

def save_deseq_diff_activity_file(deseq_outfile, store_dir, lib_short):
    df = pd.read_csv(deseq_outfile).reset_index().rename(columns={"index": "unique_id"})
    df = pd.read_csv(deseq_outfile).reset_index().rename(columns={"index": "unique_id"})
    df = pd.concat((df, df.unique_id.str.split("_", expand=True).rename(columns={0: "chrom", 1: "start", 2: "end"})), axis=1)
    df["strand"] = "."
    df = df.loc[:, ["chrom", "start", "end", "unique_id", "stat", "strand", "log2FoldChange", "baseMean", "lfcSE", "pvalue", "padj"]]

    df_induced = df.loc[((df.log2FoldChange>0) & (df.pvalue<0.05))]
    df_repressed = df.loc[((df.log2FoldChange<0) & (df.pvalue<0.05))]
    df_constitutive = df.loc[(df.pvalue>=0.05)]
    assert len(df) == sum(list(map(len, [df_induced, df_repressed, df_constitutive])))

    induced_file = os.path.join(store_dir, lib_short, "induced.bed")
    repressed_file = os.path.join(store_dir, lib_short, "repressed.bed")
    constitutive_file = os.path.join(store_dir, lib_short, "constitutive.bed")

    df_induced.to_csv(induced_file, sep="\t", header=False, index=False)
    df_repressed.to_csv(repressed_file, sep="\t", header=False, index=False)
    df_constitutive.to_csv(constitutive_file, sep="\t", header=False, index=False)
    return

###################################
# differentially active enhancers #
###################################

def get_diff_active_enhancers(
    lib_peak_file, 
    cc_peak_file, 
    lib_short, 
    store_dir
    ):
    lib_peak_bed = pybedtools.BedTool(lib_peak_file)
    cc_peak_bed = pybedtools.BedTool(cc_peak_file)
    
    induced_peak_bed = lib_peak_bed - cc_peak_bed
    repressed_peak_bed = cc_peak_bed - lib_peak_bed
    constitutive_peak_bed = lib_peak_bed + cc_peak_bed
    
    save_diff_activity_peak_file(induced_peak_bed, store_dir, lib_short, "induced")
    save_diff_activity_peak_file(repressed_peak_bed, store_dir, lib_short, "repressed")
    save_diff_activity_peak_file(constitutive_peak_bed, store_dir, lib_short, "constitutive")
    return

def get_diff_active_enhancers_deseq(
    lib_peak_file, 
    cc_peak_file,
    lib_bam_files,
    cc_bam_files, 
    lib_short,
    cc_short, 
    store_dir
    ):
    lib_peak_bed = pybedtools.BedTool(lib_peak_file)
    cc_peak_bed = pybedtools.BedTool(cc_peak_file)
    # combine peak beds
    combined_peak_bed = cc_peak_bed.cat(lib_peak_bed, postmerge=True, force_truncate=True).sort()
    # find cov of peaks against the replicate bam files
    combined_peak_cov_df = pd.concat([get_replicate_wise_cov_df(combined_peak_bed, bf, ln) for ln,bf in zip((cc_short, lib_short), (cc_bam_files, lib_bam_files))], axis=1)
    # save the cov df 
    deseq_in = save_deseq_in_df(combined_peak_cov_df, store_dir, lib_short)
    # get deseq out table
    deseq_out = get_deseq_out_file(store_dir, lib_short)
    # run deseq
    da_with_deseq(deseq_in, deseq_out)
    # save deseq results as induced, repressed, constitutive
    save_deseq_diff_activity_file(deseq_out, store_dir, lib_short)
    return


#####################
# multiprocess jobs #
#####################

def multi_args_pool_job(func, pool_iter):
    pool = mp.Pool(len(pool_iter))
    res = pool.starmap(func, pool_iter)
    pool.close()
    pool.join()
    return res
