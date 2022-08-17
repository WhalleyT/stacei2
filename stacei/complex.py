import pathlib
import anarci
import os

from itertools import chain

from .pdb_functions import PDBFile
from .anarci_functions import Anarci
from .data.data import pdb_ranges
from .chain_pairing import Pairing
from .contacts import Contacts

class TCRpMHCComplex(PDBFile, Anarci, Pairing, Contacts):
    """
    Object for defining and describing a TCR pMHC complex
    """
    def __init__(self, *args, **kwargs):
        """
        Args:
            pdb_file (str): Path to PDB file
            verbose (bool, optional): Give verbose statements. Default = False
            complex (string, optional): If given complex detection will be ignored and TCR-pMHC complex will be
                                        determined in the order of TCRa, TCRb, peptide, MHCa, MHCb
            mhc_class (string, optional): If given mhc class determination will be ignored
            outdir (string, optional): Output folder for any analysis, if ignored the default will be the basename
                                       of the PDB file
        """

        if len(args) != 1:
            assert "One positional argument required, PDB file"
        else:
            self.pdb_file = args[0]

        if not pathlib.Path(self.pdb_file).is_file():
            assert self.pdb_file, "is not a file"

        allowed_kwargs = ['verbose', 'complex', 'mhc_class']

        default_args = {
            "verbose" : False,
            "complex" : None,
            "mhc_class" : None,
            "outdir": None
        }

        for key in kwargs.keys():
            if key in allowed_kwargs:
                default_args[key] = kwargs[key]

        for arg in default_args:
            setattr(self, arg, default_args[arg])

        #single string values for complex chains
        self.tcra, self.tcrb, self.peptide, self.mhca, self.mhcb = None, None, None, None, None

        #dictionary of gene usage
        self.gene_usage = {}
        self.pdb_ranges = pdb_ranges

        #stuff for paths
        self.base_name = os.path.basename(self.pdb_file).split(".")[0]

        if self.outdir is None:
            self.outdir = os.path.basename(self.pdb_file).split(".")[0]

        self.imgt_pdb = self.outdir + "/" + self.base_name + "_IMGT.pdb"


    def annotate(self):
        if self.complex is not None:
            #we can take user input, only need to imgt number things
            complex = self.complex.split()

            if len(complex) != 5:
                assert "Complex must contain 5 characters"

            #imgt number
            sequences = self.get_sequences(self.chains)

            if self.verbose:
                print("Chains collected")
                for i in sequences:
                    print(i + ": " + sequences[i])
                print("Calling ANARCI")

            #call anarci so we can find out what TCRs there are
            self.call_anarci(self.create_anarci_obj(sequences))
            anarci_dict = {complex[0] : 0, complex[1]: 0}
            if self.verbose:
                print("Renumbering pdb file. from %s -> %s" % (self.pdb_file, self.imgt_pdb))

            if os.path.isdir(self.outdir) is False:
                os.mkdir(self.outdir)

            self.renumber_tcrs(self.anarci_result[0], anarci_dict, sequences)

            self.complex = {"TCRa": complex[0],
                            "TCRb": complex[1],
                            "peptide": complex[2],
                            "MHCa": complex[3],
                            "MHCb": complex[4]}

        else:
            # we must automatically detect classes
            if self.verbose:
                print("Automatically detecting chains and annotating them")

            self.chains = PDBFile.get_all_chains(self)

            if self.verbose:
                print("There are %i chains" %len(self.chains))
                for i in self.chains:
                    print("\t-%s" %i)

            sequences = self.get_sequences(self.chains)

            if self.verbose:
                print("Chains collected")
                for i in sequences:
                    print(i + ": " + sequences[i])
                print("Calling ANARCI")

            #call anarci so we can find out what TCRs there are
            self.call_anarci(self.create_anarci_obj(sequences))
            #anarci output order = numbering, alignment_details, hit_tables

            possible_tcra, possible_tcrb, anarci_dict = self.find_tcrs()

            # update object of possible ps and mhcs (i.e. those not tcrs)
            possible_pmhcs = set(self.chains) - set(possible_tcra) - set(possible_tcrb)

            # now renumber all the TCR chains
            if self.verbose:
                print("Renumbering pdb file. from %s -> %s" % (self.pdb_file, self.imgt_pdb))

            if os.path.isdir(self.outdir) is False:
                os.mkdir(self.outdir)

            self.renumber_tcrs(self.anarci_result[0], anarci_dict, sequences)
            tcr_pairs = self.pair_tcrs(possible_tcra, possible_tcrb)

            if len(tcr_pairs) == 0:
                assert "No pairs of TCRs found"
            else:
                if self.verbose:
                    print("There are %i TCR pairs" %len(tcr_pairs))

            mhcs = self.find_mhcs(possible_pmhcs, sequences)
            mhc_pairs, mhc_class = self.pair_mhcs(mhcs)
            if self.verbose:
                print("Print there are %i pairs of MHC detected" % len(mhc_pairs))
                print("They are MHC class %i" %mhc_class)

            possible_peptides = possible_pmhcs - set(list(chain(*mhc_pairs)))
            if self.verbose:
                print("Pairing peptide to MHC; there are %i possible peptides" % len(possible_peptides))

            pmhcs = self.pair_p_to_mhc(possible_peptides, mhc_pairs, mhc_class)
            complexes = self.pair_tcr_pmhc(tcr_pairs, pmhcs)

            if self.verbose:
                self.show_tcr_pmhcs(complexes)

            self.complex = complexes[0]
            print(self.complex)


    def contacts(self):
        self.contacts = self.calculate_contacts()
        return self.contacts










