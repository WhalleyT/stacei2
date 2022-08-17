import argparse
import sys
import os

def _parse_args():
    """
    Parse the command line arguments
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-F', dest='infile', type=str,
                        help='PDB file of a TCR-pMHC complex', required=True)
    parser.add_argument('--mhc', '-Mh', dest='mhc_class', type=int,
                        help='MHC class. if not supplied, then class will be automatically predicted',
                        required=False)
    parser.add_argument('--chains', '-C', dest='chains', type=str,
                        help='Chains of TCR-pMHC complex, in order of TCR alpha, TCR beta,'
                             'peptide, MHC alpha and MHC beta',
                        required=False)
    parser.add_argument('--ray', '-R', dest='ray_trace', action='store_true',
                        help='Flag. If provided, structures will be ray traced in Pymol. '
                             'This will affect performance at run time but will produce better images.')
    parser.add_argument('--suppress', '-S', dest='suppress', action='store_true',
                        help='Flag. If provided stdout output will be suppressed (inc. CCP4, Pymol and ANARCI)')
    parser.add_argument('--mtz', '-Mt', dest='mtz', required=False, type=str,
                        help='The MTZ file to be analysed, if None (default) is supplied it will skipped',
                        default="None")
    parser.add_argument("--vanderwaals", "-Vdw", dest='van_der_waals_distance', required=False, type=float,
                        default=4.0, help="Distance for an interaction to be considered a contact. If it does not "
                                          "satisfy other criteria it will be marked a van der Waals interaction")
    parser.add_argument("--hydrogenbond", "-Hb", dest="h_bond_distance", required=False, type=float, default=3.5,
                        help="Distance threshold required for a pair of hydrogen bond donors and acceptors to be "
                             "considered to make a hydrogen bond")
    parser.add_argument("--saltbridge", "-Sb", dest='s_bond_distance', required=False, type=float,
                        default=4.0, help="Distance threshold required for a pair of salt bridge donors and "
                                          "acceptors to be considered to make a salt bridge")
    parser.add_argument("--outdir", "-O", dest="outdir", required=False, type=str, default="",
                        help="Output directory, if not specified the output will be written to the name of the PDB")
    args = parser.parse_args()
    return args


def check_parse():
    """
    Sanitize the parse arguments
    """

    args = _parse_args()

    classes = [1, 2]

    if os.path.isfile(args.infile) is False:
        sys.exit("Input PDB could not be found")

    if args.chains is None and args.mhc_class is None:
        return args, True
    elif args.chains is not None and args.mhc_class is not None:
        if args.mhc_class not in classes:
            sys.exit("MHC class must be 1 or 2")

        if len(args.chains) != 5 and args.chains is not None:
            sys.exit("Chains argument must be exactly 5 letters")
        return args, False
    else:
        sys.exit("Chains argument and MHC argument must be supplied together or not all")