import stacei.cli
import stacei.complex

def main():
    args, auto = stacei.cli.check_parse()

    pdb_obj = stacei.complex.TCRpMHCComplex(args.infile, verbose = True)
    pdb_obj.annotate()
    #pdb_obj.contacts()
    #pdb_obj.save_contacts()
    pdb_obj.bsa()


if __name__ == "__main__":
    main()