SPACER_MIN = 5  # bp between end of SD and start codon
SPACER_MAX = 10


class InternalRBSChecker:
    """
    A class to check for internal ribosome binding sites (RBS) in a given coding DNA sequence (CDS).
    """

    def __init__(self):
        """
        Initializes the InternalRBSChecker with a predefined RBS sequence to search for.
        """
        self.rbs_like_sequences = None
        self.start_codons = None

    def initiate(self):
        """
        Initializes the list of RBS-like sequences that will be checked against the CDS.
        These sequences are derived from the Shine-Dalgarno sequence and its variants.
        """
        self.rbs_like_sequences = {
            "AGGAGG",  # perfect consensus
            "AAGGAG",
            "AGGAG",
            "GAGG",
            "AAGG",
            "AGGA",
        }
        self.start_codons = {"ATG", "GTG", "TTG"}

    def run(self, sequence: str):
        """
        Checks for internal RBS-like sequences in the provided DNA sequence.
        For each RBS-like sequence found, it looks for a start codon within a defined
        spacer region downstream of the RBS. If such a combination is found, it returns False along with the offending sequence. If no internal RBS-like sequences are found with a nearby start codon, it returns True."""
        for motif in self.rbs_like_sequences:
            idx = sequence.find(motif)
            while idx != -1:
                # look for a start codon in the spacer window
                window_start = idx + len(motif) + SPACER_MIN
                window_end = idx + len(motif) + SPACER_MAX
                window = sequence[window_start:window_end]
                for codon in self.start_codons:
                    if codon in window:
                        return False, sequence[idx : window_end + 3]
                idx = sequence.find(motif, idx + 1)
        return True, None
