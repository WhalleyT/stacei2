import warnings
from itertools import compress, chain, product
from re import sub

from .data.data import mhc_fastas, mhc_ranges, cdr_ranges

class Pairing:
    """
    Class for functions entailing TCR-pMHC chain pairing
    """
    def __init__(self):
        pass

    def pair_mhcs(self, mhcs):
        print("Pairing MHC")
        classes = []

        class_dict = {}
        alphas, betas = [], []
        mhc_pairs = []

        for mhc in mhcs:
            mhc_type = mhcs[mhc]
            if mhc_type.startswith(("DR", "DQ", "DP")):
                classes.append(2)
                class_dict[mhc] = 2
                if sub("\d", "", mhc_type.split("*")[0]).endswith("A"):
                    alphas.append(mhc)
                else:
                    betas.append(mhc)
            elif mhc_type.startswith(("A", "B", "C", "B2M")):
                classes.append(1)
                class_dict[mhc] = 1
                if mhc_type.startswith("B2M"):
                    alphas.append(mhc)
                else:
                    betas.append(mhc)
            else:
                warnings.warn("MHC class not 1 or 2")

        possible_pairs = product(alphas, betas)

        if len(set(classes)) != 1:
            warnings.warn("Multiple MHCs detected")

        mhc_class = classes[0]

        for pair in possible_pairs:
            alpha, beta = pair[0], pair[1]

            distance_constraints = None
            distances = None

            if class_dict[alpha] == 1:
                a15 = self.extract_coord(["15"], alpha)
                b23 = self.extract_coord(["23"], beta)

                dist_1 = self.mean_euclidian_distance(a15, b23)

                b104 = self.extract_coord(["104"], beta)

                dist_2 = self.mean_euclidian_distance(a15, b104)

                b51 = self.extract_coord(["51"], beta)

                dist_3 = self.mean_euclidian_distance(a15, b51)

                distance_constraints = [32, 32, 37]
                distances = [dist_1, dist_2, dist_3]
            else:
                a29 = self.extract_coord(["19"], alpha)
                b64 = self.extract_coord(["64"], beta)

                dist_1 = self.mean_euclidian_distance(a29, b64)

                a37 = self.extract_coord(["37"], beta)

                dist_2 = self.mean_euclidian_distance(a37, b64)

                b39 = self.extract_coord(["39"], beta)

                dist_3 = self.mean_euclidian_distance(a37, b39)

                distance_constraints = [34, 22, 28]
                distances = [dist_1, dist_2, dist_3]

            paired = True
            for dist, target in zip(distances, distance_constraints):
                if dist > target:
                    paired = False

            if paired:
                mhc_pairs.append(pair)

        return mhc_pairs, mhc_class

    def pair_p_to_mhc(self, possible_peptides, mhc_pairs, mhc_class):
        pmhcs = []
        for peptide in possible_peptides:
            contacts = []

            for mhc in mhc_pairs:
                if mhc_class == 1:
                    residue_range = list(mhc_ranges["class_1_MHCa1"]) + list(mhc_ranges["class_1_MHCa2"])
                    mhc_residues = self.extract_coord(residue_range, mhc[0])
                else:
                    alpha = self.extract_coord(list(mhc_ranges["class_1_MHCa1"]), mhc[0])
                    beta = self.extract_coord(list(mhc_ranges["class_2_MHCb1"]), mhc[1])
                    mhc_residues = alpha + beta

                # todo renumber peptide 1 to n
                peptide_residues = self.extract_coord(list(range(0, 1000)), peptide)

                ncont = 0
                for i in mhc_residues:
                    for j in peptide_residues:
                        dist = self.euclidian_distance(i, j)
                        if dist < 4:
                            ncont += 1

                contacts.append(ncont)

            nonzero = False
            for c in contacts:
                if c > 0:
                    nonzero = True

            if nonzero:
                idx = contacts.index(max(contacts))
                mhc = mhc_pairs[idx]
                pmhcs.append({"MHCa": mhc[0],
                                   "MHCb": mhc[1],
                                   "peptide": peptide})
            else:
                warnings.warn("No cognate MHC found for peptide %s" % peptide)
        return pmhcs


    def find_mhcs(self, possible_pmhcs, sequences):
        class_one = ["A*", "B*", "C*"]

        db = mhc_fastas

        mhcs = {}

        for chain in possible_pmhcs:
            sequence = sequences[chain]

            scores = {}

            if len(sequence) > 50:
                for name in db:
                    db_seq = db[name]
                    blast = self.protein_blast(sequence, db_seq)
                    scores[name] = blast.score

                max_key = max(scores, key=scores.get)
                score = scores[max_key]

                #todo check what a reasonable cut off here is
                if score > 250 and max_key is not None:
                    if max_key == "B2M":
                        hla_type = max_key
                    else:
                        hla_type = max_key.split()[1].split(":")[0]

                    mhcs[chain] = hla_type

        return mhcs

    def annotate_pmhc(self):
        self.mhcs = self.find_mhcs()
        self.pair_mhcs()
        print("Print there are %i pairs of MHC detected" %len(self.mhc_pairs))
        self.possible_peptides = self.possible_pmhcs - set(list(chain(*self.mhc_pairs)))
        print("Pairing peptide to MHC; there are %i possible peptides" %len(self.possible_peptides))
        self.pair_p_to_mhc()

    def pair_tcrs(self, possible_tcra, possible_tcrb):
        """
        takes two sets of strings and returns pairs of TCRs
        based on their distance of Cys104 to each other
        :return: list of pairs [(A, B)...]
        """

        tcr_pairs = []
        possible_pairs = list(product(possible_tcra, possible_tcrb))

        for pair in possible_pairs:
            tcra, tcrb = pair[0], pair[1]

            tcra_104 = self.extract_coord(["104"], tcra)
            tcrb_104 = self.extract_coord(["104"], tcrb)

            ed = self.mean_euclidian_distance(tcra_104, tcrb_104)
            if ed < 23:
                tcr_pairs.append(pair)

        return tcr_pairs

    def pair_tcr_pmhc(self, tcr_pairs, pmhcs):
        tcr_pmhcs = []
        for tcr in tcr_pairs:
            tcra, tcrb = tcr[0], tcr[1]

            alpha = self.extract_coord(list(cdr_ranges["CDR3a"]), tcra)
            beta = self.extract_coord(list(cdr_ranges["CDR3b"]), tcrb)
            tcr_res = alpha + beta

            contacts = list()

            for pmhc in pmhcs:
                peptide = pmhc["peptide"]

                # todo renumber peptide 1 to n
                peptide_residues = self.extract_coord(list(range(0, 1000)), peptide)

                ncont = 0
                for i in tcr_res:
                    for j in peptide_residues:
                        dist = self.euclidian_distance(i, j)
                        if dist < 4:
                            ncont += 1

                contacts.append(ncont)

            nonzero = False
            for c in contacts:
                if c > 0:
                    nonzero = True

            if nonzero:
                idx = contacts.index(max(contacts))
                mhc = pmhcs[idx]
                tcr_pmhcs.append({"TCRa": tcra,
                                       "TCRb": tcrb,
                                       "MHCa": mhc["MHCa"],
                                       "MHCb": mhc["MHCb"],
                                       "peptide": mhc["peptide"]})
            else:
                warnings.warn("No cognate TCRs found for peptide %s" % peptide)

            return tcr_pmhcs

    def show_tcr_pmhcs(self, complexes):
        print("There are %i TCR-pMHC complexes found" % len(complexes))
        print("They are:")
        for idx, complex in enumerate(complexes):
            num = idx + 1
            print("Complex no. %i" % num)
            print("\tTCRa: %s (TRAV = %s; TRAJ = %s)" % (
            complex["TCRa"], self.gene_usage["TRAV"], self.gene_usage["TRAJ"]))
            print("\tTCRb: %s (TRBV = %s; TRBJ = %s)" % (
            complex["TCRb"], self.gene_usage["TRBV"], self.gene_usage["TRBJ"]))
            print("\tpeptide: %s" % complex["peptide"])
            print("\tMHCa: %s" % complex["MHCa"])
            print("\tMHCb: %s" % complex["MHCb"])
            print()