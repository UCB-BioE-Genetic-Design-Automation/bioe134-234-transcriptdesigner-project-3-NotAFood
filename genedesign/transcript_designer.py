from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.models.transcript import Transcript
from genedesign.rbs_chooser import RBSChooser
from genedesign.seq_utils.codon_usage import CodonUsage

MAX_ATTEMPTS = 1000


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.codonUsage = CodonUsage()
        self.rbsChooser = None
        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()
        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()
        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

    def initiate(self) -> None:
        """
        Initializes the codon usage table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        self.codonUsage.initiate()

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        Uses a sliding window approach: starts with highest-CAI codons, then selectively
        randomizes 3-codon windows until all validation checks pass.

        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.

        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.

        Raises:
            RuntimeError: If no valid sequence can be generated after MAX_ATTEMPTS attempts.
        """
        # Generate initial sequence using best CAI codons
        codons = self._generate_best_cai_sequence(peptide)

        # Window-based refinement: validate and randomize windows as needed
        window_size = 3
        num_windows = (len(peptide) + window_size - 1) // window_size

        for attempt in range(MAX_ATTEMPTS):
            dna_sequence = "".join(codons)
            if self._all_checkers_pass(dna_sequence):
                # Found valid sequence
                codons.append("TAA")
                cds = "".join(codons)
                selectedRBS = self.rbsChooser.run(cds, ignores)
                return Transcript(selectedRBS, peptide, codons)

            # Randomize a window in a cycling pattern to try different regions
            window_idx = attempt % num_windows
            start_idx = window_idx * window_size
            end_idx = min(start_idx + window_size, len(peptide))

            for i in range(start_idx, end_idx):
                codons[i] = self.codonUsage.weighted_random_codon(peptide[i])

        raise RuntimeError(
            "Could not generate a valid DNA sequence after maximum attempts."
        )

    def _generate_best_cai_sequence(self, peptide: str) -> list:
        """
        Generates a list of codons for the given peptide sequence using the best (highest CAI) codon for each amino acid.
        This provides an optimized starting point for the sliding window refinement.

        Parameters:
            peptide (str): The protein sequence to translate.
        Returns:
            list: A list of codons corresponding to the peptide sequence.
        """
        return [self.codonUsage.best_cai_codon(aa) for aa in peptide]

    def _all_checkers_pass(self, sequence: str) -> bool:
        """
        Placeholder for sequence validation checks. Implement specific checks as needed.

        Parameters:
            sequence (str): The DNA sequence to validate.
        Returns:
            bool: True if all checks pass, False otherwise.
        """
        # Note: The checkers return false if a problem is found, so we can directly return the result of the checks.
        has_internal_promoter, _ = self.promoterChecker.run(sequence)
        has_hairpin, _ = hairpin_checker(sequence)
        has_forbidden, _ = self.forbiddenChecker.run(sequence)

        return has_internal_promoter and has_hairpin and has_forbidden


if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    import time

    peptide = "MYPFIRTARMTV"

    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    start_time = time.time()
    transcript = designer.run(peptide, ignores)
    elapsed_time = time.time() - start_time

    # Print out the transcript information
    print(transcript)
    print(f"\nCompleted in {elapsed_time:.2f} seconds")
