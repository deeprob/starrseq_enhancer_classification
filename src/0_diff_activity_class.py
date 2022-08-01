import argparse
import utils as ut

def main(
    lib_short, 
    cc_short,
    peak_dir,
    store_dir
    ):
    lib_peak_file = ut.get_lib_peak_filepath(peak_dir, lib_short)
    cc_peak_file = ut.get_lib_peak_filepath(peak_dir, cc_short)
    ut.get_diff_active_enhancers(lib_peak_file, cc_peak_file, lib_short, store_dir)
    return






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq differentially active peaks")
    parser.add_argument("meta_file", type=str, help="The file path to the meta file")
    parser.add_argument("lib", type=str, help="library name as given in the meta file")
    parser.add_argument("control_lib", type=str, default="", help="library name of the control line as given in the meta file, compare induced,repressed and constitutive peaks between lib1 and cc")
    parser.add_argument("peak_dir", type=str, help="Dir where peak files are stored")
    parser.add_argument("store_dir", type=str, help="Output dir where results will be stored")

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)
    cc_args = ut.create_args(cli_args.meta_file, cli_args.control_lib)

    main(
        lib_args.library_short,
        cc_args.library_short,
        cli_args.peak_dir,
        cli_args.store_dir
    )
