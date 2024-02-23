"""
Script to run full analysis for D resonances
"""

import argparse
import os

def run_full_analysis(config,
                      an_res_file,
                      suffix,
                      doEP,
                      skip_resolution,
                      skip_projection,
                      skip_rawyield
                      ):
    """
    function for full analysis

    Parameters
    ----------

    - config (str): path of directory with config files
    - an_res_file (str): path of directory with analysis results
    - suffix (str): suffix for output files
    - doEP (bool): do EP resolution
    - skip_resolution (bool): skip resolution extraction
    - skip_projection (bool): skip projection extraction
    - skip_rawyield (bool): skip raw yield extraction
    """

    # get all parameters needed
    if suffix != "":
        suffix_withopt = f" -s {suffix}"
    else:
        suffix = an_res_file.split("AnalysisResults")[-1].replace(".root", "")
        suffix_withopt = f" -s {suffix}"
    if doEP:
        outputdir = "~/scalarProd/ep"
    else:
        outputdir = "~/scalarProd/sp"

    if not skip_resolution:
        # resolution extraction
        outputdir_reso = f"-o {outputdir}/resolution"
        command_reso = f"python3 compute_reso.py {config} {an_res_file} {suffix_withopt} {outputdir_reso}"
        if doEP:
            command_reso += " --doEP True"
        print("\n\033[92m Starting resolution extraction\033[0m")
        print(f"\033[92m {command_reso}\033[0m")
        os.system(command_reso)

    if not skip_projection:
        # projection
        outputdir_proj = f"-o {outputdir}/proj"
        if doEP:
            command_proj = f"python3 compute_event_plane.py {config} {an_res_file} {suffix_withopt} {outputdir_proj}"
        else:
            command_proj = f"python3 compute_scalar_product.py {config} {an_res_file} {suffix_withopt} {outputdir_proj}"
        print("\n\033[92m Starting projection\033[0m")
        print(f"\033[92m {command_proj}\033[0m")
        os.system(command_proj)

    if not skip_rawyield:
        # raw yield
        outputdir_rawyield = f"-o {outputdir}/ry"
        proj_file = f"{outputdir}/proj/"
        if doEP:
            proj_file += f"event_plane{suffix}.root"
            #command_rawyield = f"python3 invmass_fitter.py {config} {proj_file} {suffix_withopt} {outputdir_rawyield}"
            command_rawyield_out = f"python3 GetRawYieldsDplusDs.py {config} k3050 {proj_file} {outputdir_rawyield} {suffix_withopt} --isOutOfPlane"
            command_rawyield_in = f"python3 GetRawYieldsDplusDs.py {config} k3050 {proj_file} {outputdir_rawyield} {suffix_withopt} --isInPlane"
        else:
            print("\n\033[92m The raw yield extraction is not implemented for SP\033[0m")
            return
        print("\n\033[92m Starting raw yield extraction\033[0m")
        print(f"\033[92m {command_rawyield_in}\033[0m")
        os.system(command_rawyield_in)
        os.system(command_rawyield_out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument("an_res_file", metavar="text",
                        default="an_res.root", help="input ROOT file with anres")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--doEP", action="store_true", default=False,
                        help="do EP resolution")
    parser.add_argument("--skip_resolution", action="store_true", default=False,
                        help="skip resolution extraction")
    parser.add_argument("--skip_projection", action="store_true", default=False,
                        help="skip projection extraction")
    parser.add_argument("--skip_rawyield", action="store_true", default=False,
                        help="skip raw yield extraction")
    args = parser.parse_args()

    run_full_analysis(
        args.config,
        args.an_res_file,
        args.suffix,
        args.doEP,
        args.skip_resolution,
        args.skip_projection,
        args.skip_rawyield
    )
