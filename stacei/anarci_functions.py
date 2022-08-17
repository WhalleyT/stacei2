import anarci
import swalign

from itertools import  compress, chain
from warnings import warn

from .generic_functions import Generic

class Anarci(Generic):
    """
    Helper class to run ANARCI
    """
    def __init__(self):
        pass

    def create_anarci_obj(self, sequences):
        """
        creates an ANARCI input list
        :return: list of pairs
        """
        anarci_obj = []

        for elem in sequences:
            pair = (elem, sequences[elem])
            anarci_obj.append(pair)

        return anarci_obj

    def call_anarci(self, sequences):
        self.anarci_result = anarci.anarci(sequences, scheme="imgt", assign_germline=True)

    def find_tcrs(self):
        """
        Function that detects whether a sequence is detected as a TCR or not
        :param alignment: ANARCI alignment object, a list containing a dict
        :return: two lists, returned to self
        """

        possible_tcr_a, possible_tcr_b = set(), set()
        anarci_dict = {}
        for i, j in zip(self.anarci_result[1], self.chains):
            if i is not None:
                # then we have detected something
                # so now detect if there are TCRs
                tcr_bool = [x["chain_type"] in ["A", "B"] for x in i]
                tcrs = list(compress(i, tcr_bool))

                # now we have a list of possible tcr hits; take the lowest e value
                evals = [x["evalue"] for x in tcrs]
                idx = evals.index(min(evals))
                tcr_hit = tcrs[idx]

                if tcr_hit["chain_type"] == "A":
                    possible_tcr_a.add(tcr_hit["query_name"])
                    self.gene_usage["TRAV"] = tcr_hit["germlines"]["v_gene"][0][1]
                    self.gene_usage["TRAJ"] = tcr_hit["germlines"]["j_gene"][0][1]
                    anarci_dict[j] = idx
                elif tcr_hit["chain_type"] == "B":
                    possible_tcr_b.add(tcr_hit["query_name"])
                    anarci_dict[j] = idx
                    self.gene_usage["TRBV"] = tcr_hit["germlines"]["v_gene"][0][1]
                    self.gene_usage["TRBJ"] = tcr_hit["germlines"]["j_gene"][0][1]
                else:
                    warn("Warning: Incorrect parsing of ANARCI output of %s" % j)

        return possible_tcr_a, possible_tcr_b, anarci_dict


    def extract_num_seq(self, anarci_num):
        """
        extract sequence and anarci numbering
        :param anarci_num: anarci numbering object
        :return: a list of tuples containing a residue
        and amino acid and then a string of seq
        """
        sequence_list = []
        num_list = []

        for i in anarci_num:
            num = i[0]
            seq = i[1]

            if str(num[1]) != " ":
                num = str(num[0]) + num[1]
            else:
                num = str(num[0])

            if seq != "-":
                sequence_list.append(seq)
                num_list.append(num)

        sequence_string = "".join(sequence_list)
        return list(zip(num_list, sequence_string)), sequence_string


    def renumber_tcrs(self, anarci_num, anarci_dict, sequences):
        renumbered_tcrs = {}
        for i, j in zip(anarci_num, self.chains):
            if j in anarci_dict:
                idx = anarci_dict[j]

                #collect best number and sequence from original
                best_num = i[idx]
                seq = sequences[j]
                num_seq = best_num[0]
                number_list, anarci_seq = self.extract_num_seq(num_seq)

                blast = self.protein_blast(seq, anarci_seq)

                offset = blast.r_pos

                #todo implement a way of dealing with offset; for now we will proceed naively
                residue_list = []
                residue = None
                for index, original_aa in enumerate(seq):
                    #todo take the original residue and make a pair (orig, new) so we can index
                    if index < len(number_list):
                        residue = str(number_list[index][0])
                    else:
                        residue = str(int(residue) + 1)

                    residue_list.append(residue)

                renumbered_tcrs[j] = residue_list

                self.rewrite_pdb_nums(renumbered_tcrs, file=self.imgt_pdb)