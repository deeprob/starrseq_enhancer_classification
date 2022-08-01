import argparse
import utils as ut


def main(
    lib_short,
    lib_prefix,
    lib_reps, 
    cc_short,
    cc_prefix,
    cc_reps,
    peak_dir,
    bam_dir,
    store_dir
    ):
    lib_peak_file = ut.get_lib_peak_filepath(peak_dir, lib_short)
    cc_peak_file = ut.get_lib_peak_filepath(peak_dir, cc_short)
    lib_bam_files = ut.get_rep_bam_files(bam_dir, lib_short, lib_prefix, lib_reps)
    cc_bam_files = ut.get_rep_bam_files(bam_dir, cc_short, cc_prefix, cc_reps)
    ut.get_diff_active_enhancers_deseq(lib_peak_file, cc_peak_file, lib_bam_files, cc_bam_files, lib_short, cc_short, store_dir)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq differentially active peaks")
    parser.add_argument("meta_file", type=str, help="The file path to the meta file")
    parser.add_argument("lib", type=str, help="library name as given in the meta file")
    parser.add_argument("control_lib", type=str, default="", help="library name of the control line as given in the meta file, compare induced,repressed and constitutive peaks between lib1 and cc")
    parser.add_argument("peak_dir", type=str, help="Dir where peak files are stored")
    parser.add_argument("bam_dir", type=str, help="Dir where bam files are stored")
    parser.add_argument("store_dir", type=str, help="Output dir where results will be stored")

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)
    cc_args = ut.create_args(cli_args.meta_file, cli_args.control_lib)

    main(
        lib_args.library_short,
        lib_args.library_prefix,
        lib_args.library_reps,
        cc_args.library_short,
        cc_args.library_prefix,
        cc_args.library_reps,
        cli_args.peak_dir,
        cli_args.bam_dir,
        cli_args.store_dir
    )
