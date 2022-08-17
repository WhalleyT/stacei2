import warnings
import re

import numpy as np


class PDBFile:
    def __init__(self):
        pass

    def three2one(self, three):
        """
        Function to convert the 3 letter amino acid code to the 1 letter type
        :param three: three letter amino acid code
        :return: single letter string of amino acid
        """
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        if len(three) % 3 != 0:
            raise ValueError('Input length should be a multiple of three')

        if three not in d:
            raise ValueError("Input is not a valid three letter amino acid code")

        return d[three]

    def get_column(self, path, key):
        list_of_col = []

        with open(path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    indexes = list(self.pdb_ranges[key])
                    contents = "".join([line[x] for x in indexes]).strip()
                    list_of_col.append(contents)

        return list_of_col


    def get_all_chains(self):
        chains = set(self.get_column(self.pdb_file, "CHAIN"))
        return sorted(list(set(chains)))


    def get_sequence(self, chain):
        sequence = ""

        chains = self.get_column(self.pdb_file, "CHAIN")
        three_letter_codes = self.get_column(self.pdb_file, "RESIDUE_SEQID")
        single_letter_codes = [self.three2one(x) for x in three_letter_codes]
        resnums = self.get_column(self.pdb_file, "RESIDUE_NUM")

        previous_residue = None
        for i, j, k in zip(chains, single_letter_codes, resnums):
            if i == chain:
                if previous_residue != k:
                    sequence += j
                previous_residue = k

        return sequence


    def get_sequences(self, chains):
        d ={}

        for chain in chains:
            d[chain] = self.get_sequence(chain)

        return d


    def pad_number(self, number):
        """
        function to pad out a residue number to fit in the PDB column correctly
        :param number: string of resnum. Should be 4 digits/chars max
        :return: string of number with additional spaces
        """

        colsize = 4
        size = len(number)
        if size > 4:
            warnings.warn("Residue number greater than 4 digits, this will not fit PDB format")

        number =  " " * (colsize - size) + number
        return number


    def rewrite_pdb_nums(self, residue_dictionary, file = "tmp.pdb"):
        """
        rewrites a pdb_file. input file is hardcoded to original
        :param residue_dictionary: key = chain, value = residues
        :param file: outputfile name
        :return: None, side effect writes file
        """
        residue_range = self.pdb_ranges["RESIDUE_NUM"]

        out = open(file, "w")

        for chain in residue_dictionary:
            list_of_residues = residue_dictionary[chain]
            residue_idx = 0
            previous = None
            current_residue = None

            with open(self.pdb_file) as f:
                for line in f:
                    if line.startswith("ATOM"):
                        chain_in_pdb = list(self.pdb_ranges["CHAIN"])
                        chain_in_pdb = "".join([line[x] for x in chain_in_pdb]).strip()
                        if chain_in_pdb == chain:
                            #remember that pdb residue is right justified
                            residue_no = "".join([line[x] for x in residue_range]).strip()
                            if residue_no != previous:
                                current_residue = list_of_residues[residue_idx]
                                residue_idx += 1

                            padded_residue = self.pad_number(current_residue)

                            if len(padded_residue) != len(residue_range):
                                warnings.warn("Residue number padding does not match PDB width")

                            list_line = list(line)
                            padded_residue_list = list(padded_residue)
                            for i, j in zip(residue_range, padded_residue_list):
                                list_line[i] = j
                            line = "".join(list_line)

                            previous = residue_no

                    out.write(line)

    def extract_coord(self, residues, chain):
        out = []
        for residue in residues:
            with open(self.imgt_pdb) as f:
                for line in f:
                    if line.startswith("ATOM"):
                        chain_in_pdb = list(self.pdb_ranges["CHAIN"])
                        chain_in_pdb = "".join([line[x] for x in chain_in_pdb]).strip()

                        if chain == chain_in_pdb:
                            residue_in_pdb = list(self.pdb_ranges["RESIDUE_NUM"])
                            residue_in_pdb = "".join([line[x] for x in residue_in_pdb]).strip()

                            if str(residue) == re.sub("\D", "",residue_in_pdb).strip():
                                x_coord = list(self.pdb_ranges["X_COORD"])
                                x_coord = "".join([line[x] for x in x_coord]).strip()

                                y_coord = list(self.pdb_ranges["Y_COORD"])
                                y_coord = "".join([line[x] for x in y_coord]).strip()

                                z_coord = list(self.pdb_ranges["Z_COORD"])
                                z_coord = "".join([line[x] for x in z_coord]).strip()

                                xyz = tuple([float(x_coord), float(y_coord), float(z_coord)])
                                out.append(xyz)
        return out

    def mean_euclidian_distance(self, x, y):
        """
        computes euclidian distance between two lists of tuples (representing 3D space)
        :param x: list of tuples 1
        :param y: list of tuples 2
        :return: float
        """
        dists = list()

        for i in x:
            for j in y:
                i = np.array(i)
                j = np.array(j)

                dist = np.linalg.norm(i-j)
                dists.append(dist)

        return np.mean(np.array(dists))

    def euclidian_distance(self, x, y):
        i = np.array(x)
        j = np.array(y)

        return np.linalg.norm(i - j)