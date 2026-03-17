from itertools import product

from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.models.transcript import Transcript
from genedesign.rbs_chooser import RBSChooser
from genedesign.seq_utils.codon_usage import CodonUsage


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
        Uses a sliding window optimization approach: processes 3-amino-acid windows left to right,
        enumerating all codon combinations for each window while considering upstream (preamble)
        and downstream context. Selects the best combination based on validation and CAI scoring.

        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.

        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.

        Raises:
            RuntimeError: If no valid sequence can be generated.
        """
        # Initialize the chosen codons list (will be built incrementally)
        chosen_codons = []

        # Window size: 3 amino acids at a time
        window_size = 3
        lookahead_size = 6  # Next 6 amino acids for downstream context

        # Process peptide in windows of 3 amino acids
        for window_start in range(0, len(peptide), window_size):
            window_end = min(window_start + window_size, len(peptide))
            current_window = peptide[window_start:window_end]

            # Get downstream context (next 6 amino acids)
            downstream_start = window_end
            downstream_end = min(downstream_start + lookahead_size, len(peptide))
            downstream_aas = peptide[downstream_start:downstream_end]

            # Find the best codon combination for this window
            best_codons = self._optimize_window(
                current_window, chosen_codons, downstream_aas, ignores
            )

            # Add the best codons to our chosen list
            chosen_codons.extend(best_codons)

        # Add stop codon
        chosen_codons.append("TAA")

        # Final assembly: select RBS and create transcript
        cds = "".join(chosen_codons)
        selectedRBS = self.rbsChooser.run(cds, ignores)

        return Transcript(selectedRBS, peptide, chosen_codons)

    def _optimize_window(
        self, window_aas: str, preamble_codons: list, downstream_aas: str, ignores: set
    ) -> list:
        """
        Optimizes a 3-amino-acid window by enumerating all possible codon combinations
        and selecting the best one based on validation checks and CAI.

        Parameters:
            window_aas (str): The amino acids in the current window (1-3 aa).
            preamble_codons (list): Already-chosen codons from previous windows.
            downstream_aas (str): Next 6 amino acids for downstream context.
            ignores (set): RBS options to ignore.

        Returns:
            list: The best codon combination for this window.
        """
        # Get all possible codons for each amino acid in the window
        codon_options = [self._get_all_codons_for_aa(aa) for aa in window_aas]

        # Generate all possible combinations (cartesian product)
        all_combinations = list(product(*codon_options))

        # Generate downstream codons using best CAI (for context only)
        downstream_codons = [
            self.codonUsage.best_cai_codon(aa) for aa in downstream_aas
        ]

        # Score all combinations
        best_score = -float("inf")
        best_combination = None

        for combination in all_combinations:
            score = self._score_combination(
                list(combination), preamble_codons, downstream_codons, ignores
            )

            if score > best_score:
                best_score = score
                best_combination = combination

        # If we found a good combination, return it
        # Otherwise return the best CAI combination as fallback
        if best_combination is not None:
            return list(best_combination)
        else:
            # Fallback: use best CAI codons
            return [self.codonUsage.best_cai_codon(aa) for aa in window_aas]

    def _get_all_codons_for_aa(self, aa: str) -> list:
        """
        Returns all possible codons for a given amino acid.

        Parameters:
            aa (str): Single-letter amino acid code.

        Returns:
            list: All codons that encode this amino acid.
        """
        if aa not in self.codonUsage.amino_acid_to_codons:
            return ["ATG"]  # Default fallback

        # Extract just the codon strings from (codon, weight) tuples
        return [codon for codon, _ in self.codonUsage.amino_acid_to_codons[aa]]

    def _score_combination(
        self,
        current_codons: list,
        preamble_codons: list,
        downstream_codons: list,
        ignores: set,
    ) -> float:
        """
        Scores a codon combination based on validation checks and CAI.
        Higher scores are better.

        Parameters:
            current_codons (list): Codons for the current window being optimized.
            preamble_codons (list): Already-chosen codons from previous windows.
            downstream_codons (list): Placeholder codons for downstream context.
            ignores (set): RBS options to ignore.

        Returns:
            float: Score for this combination (higher is better).
        """
        # Construct the full sequence for validation
        # preamble + current + downstream (no stop codon yet in middle of sequence)
        all_codons = preamble_codons + current_codons + downstream_codons
        cds = "".join(all_codons)

        # Select RBS to construct full transcript
        try:
            selectedRBS = self.rbsChooser.run(cds, ignores)
            full_sequence = selectedRBS.utr.upper() + cds
        except Exception:
            # If RBS selection fails, return very low score
            return -1000.0

        # Check validation for the full context
        has_internal_promoter, _ = self.promoterChecker.run(full_sequence)
        has_hairpin, _ = hairpin_checker(full_sequence)
        has_forbidden, _ = self.forbiddenChecker.run(full_sequence)

        # Calculate CAI just for the current window being optimized
        _, _, _, window_cai = self.codonChecker.run(current_codons)

        # Scoring strategy:
        # - If all checks pass: reward with CAI value (0-1 range)
        # - If checks fail: penalize but still differentiate by CAI
        if has_internal_promoter and has_hairpin and has_forbidden:
            # Valid sequence: high base score + CAI bonus
            return 100.0 + window_cai
        else:
            # Invalid sequence: still use CAI to differentiate, but keep score low
            return window_cai


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
