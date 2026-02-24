import random
from pathlib import Path

from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.models.transcript import Transcript
from genedesign.rbs_chooser import RBSChooser

MAX_ATTEMPTS = 1000


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.aminoAcidToCodons = {}  # Maps amino acid to list of (codon, weight) tuples
        self.rbsChooser = None
        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()
        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()
        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # Load codon usage data for weighted random selection
        self.aminoAcidToCodons = {}
        codon_usage_path = Path(__file__).parent / "data" / "codon_usage.txt"

        with open(codon_usage_path, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    codon = parts[0]
                    amino_acid = parts[1]
                    weight = float(parts[2])

                    if amino_acid not in self.aminoAcidToCodons:
                        self.aminoAcidToCodons[amino_acid] = []

                    self.aminoAcidToCodons[amino_acid].append((codon, weight))

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.

        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.

        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """

        dna_sequence_candidate = ""
        for attempt in range(MAX_ATTEMPTS):
            codons = self._generate_weighted_random_sequence_of_codons(peptide)

            dna_sequence_candidate = "".join(codons)
            if self._all_checkers_pass(dna_sequence_candidate):
                break  # Valid sequence found, exit the loop
        else:
            # fallback if no valid solution found
            raise RuntimeError(
                "Could not generate a valid DNA sequence after maximum attempts."
            )

        # Append the stop codon (TAA in this case)
        codons.append("TAA")

        # Build the CDS from the codons
        cds = "".join(codons)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)

    def _generate_weighted_random_sequence_of_codons(self, peptide: str) -> list:
        """
        Generates a list of codons for the given peptide sequence using weighted random selection.

        Parameters:
            peptide (str): The protein sequence to translate.
        Returns:
            list: A list of codons corresponding to the peptide sequence.
        """
        codons = []
        for amino_acid in peptide:
            codon = self._weighted_random_codon(amino_acid)
            codons.append(codon)
        return codons

    def _weighted_random_codon(self, amino_acid: str) -> str:
        """
        Returns a randomly selected codon for the given amino acid, weighted by codon usage frequency.

        Parameters:
            amino_acid (str): The amino acid to translate.

        Returns:
            str: The corresponding codon, selected based on weighted probabilities.
        """
        if amino_acid not in self.aminoAcidToCodons:
            # Fallback to highest CAI codon if amino acid not found
            return self.aminoAcidToCodon.get(amino_acid, "ATG")

        codons_with_weights = self.aminoAcidToCodons[amino_acid]
        codons = [codon for codon, _ in codons_with_weights]
        weights = [weight for _, weight in codons_with_weights]

        # Use random.choices to select a codon based on weights
        selected_codon = random.choices(codons, weights=weights, k=1)[0]
        return selected_codon

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
    peptide = "MYPFIRTARMTV"

    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)

    # Print out the transcript information
    print(transcript)
