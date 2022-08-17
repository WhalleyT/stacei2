import swalign

class Generic:
    def __init__(self):
        pass

    def protein_blast(self, seq1, seq2):
        """
        generic caller to swalign, takes two sequences; SW aligns and returns swalign obj
        :param seq1: query
        :param seq2: referemce
        :return: swalign object
        """
        MATCH = 2
        MISMATCH = -1
        SCORE = swalign.NucleotideScoringMatrix(MATCH, MISMATCH)

        sw = swalign.LocalAlignment(SCORE)

        return sw.align(seq1, seq2)