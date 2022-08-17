import os.path

import pandas as pd
import tqdm
from joblib import Parallel, delayed, cpu_count

from .data.data import pdb_ranges

class Contacts:
    """
    Class for calculating contacts
    """

    def parallel_dist(self, i, tcr_pmhc, pairs, line_tuple, dist_cut):
        collect = []
        chain_to = i[2]
        to_chain_name = list(tcr_pmhc.keys())[list(tcr_pmhc.values()).index(chain_to)]
        for j in line_tuple:
            chain_from = j[2]
            from_chain_name = list(tcr_pmhc.keys())[list(tcr_pmhc.values()).index(chain_from)]
            if (to_chain_name, from_chain_name) in pairs or (from_chain_name, to_chain_name) in pairs:
                dist = self.euclidian_distance(i[-3:], j[-3:])
                if dist <= dist_cut:
                    out = list(i) + list(j) + [dist]
                    collect.append(out)
        return collect

    def calculate_contacts(self, distance = 4.0):
        #todo allow for user flag for ncpu
        n_cpu = cpu_count()

        line_tuple = []
        tcr_pmhc = self.complex

        pairs = {("TCRa", "peptide"), ("TCRb", "peptide"), ("TCRa", "MHCa"), ("TCRb", "MHCa"), ("TCRa", "MHCb"),
                 ("TCRb", "MHCb"), ("peptide", "MHCa"), ("peptide", "MHCb")}

        if os.path.exists(self.imgt_pdb) is False:
            assert "No IMGT numbered PDB file exists, run the annotate() method first"

        with open(self.imgt_pdb) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain = list(pdb_ranges["CHAIN"])
                    chain = "".join([line[x] for x in chain]).strip()

                    x_coord = list(pdb_ranges["X_COORD"])
                    x_coord = "".join([line[x] for x in x_coord]).strip()

                    y_coord = list(pdb_ranges["Y_COORD"])
                    y_coord = "".join([line[x] for x in y_coord]).strip()

                    z_coord = list(pdb_ranges["Z_COORD"])
                    z_coord = "".join([line[x] for x in z_coord]).strip()

                    serial = list(pdb_ranges["SERIAL"])
                    serial = "".join([line[x] for x in serial]).strip()

                    atom = list(pdb_ranges["ATOM_NAME"])
                    atom = "".join([line[x] for x in atom]).strip()

                    big_tuple = (serial, atom, chain,
                                 float(x_coord), float(y_coord), float(z_coord))
                    line_tuple.append(big_tuple)

        out_list = []
        dnames = ["donor_atom_no", "donor_atom", "donor_chain", "donor_x_coord", "donor_y_coord", "donor_z_coord"]
        anames = ["acceptor_atom_no", "acceptor_atom", "acceptor_chain", "acceptor_x_coord", "acceptor_y_coord", "acceptor_z_coord"]
        names = dnames + anames + ["distance"]

        print("Calculating all contacts between TCR-pMHC residues")
        print("Using %i cores" %n_cpu)
        out_list = Parallel(n_jobs=n_cpu, verbose=5)(delayed(self.parallel_dist)(i, tcr_pmhc, pairs, line_tuple, distance)
                                                     for i in line_tuple)

        flat_out = []

        for i in out_list:
            for j in i:
                flat_out.append(j)

        self.contacts = pd.DataFrame(flat_out, columns=names)
        return self.contacts

    def save_contacts(self, file = None):
        #hacky way to get file as a self arg
        if file is None:
            file = self.outdir + "/" + self.base_name + "_contacts.csv"
        self.contacts.to_csv(file)