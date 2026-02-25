import random
from pathlib import Path


class CodonUsage:
    """
    Loads codon usage data and provides codon selection for a given amino acid.

    Attributes:
        amino_acid_to_codons (dict): Maps each amino acid to a list of (codon, weight) tuples.
    """

    def initiate(self) -> None:
        """
        Loads codon usage data from the codon_usage.txt file.
        """
        codon_usage_path = Path(__file__).parent.parent / "data" / "codon_usage.txt"
        self.amino_acid_to_codons = {}
        with open(codon_usage_path, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    codon, amino_acid, weight = parts[0], parts[1], float(parts[2])
                    self.amino_acid_to_codons.setdefault(amino_acid, []).append((codon, weight))

    def best_cai_codon(self, amino_acid: str) -> str:
        """
        Returns the codon with the highest CAI weight for the given amino acid.

        Parameters:
            amino_acid (str): Single-letter amino acid code.

        Returns:
            str: The codon with the highest usage weight.
        """
        if amino_acid not in self.amino_acid_to_codons:
            return "ATG"
        return max(self.amino_acid_to_codons[amino_acid], key=lambda x: x[1])[0]

    def weighted_random_codon(self, amino_acid: str) -> str:
        """
        Returns a randomly selected codon for the given amino acid, weighted by usage frequency.

        Parameters:
            amino_acid (str): Single-letter amino acid code.

        Returns:
            str: A codon selected by weighted random sampling.
        """
        if amino_acid not in self.amino_acid_to_codons:
            return "ATG"
        codons_with_weights = self.amino_acid_to_codons[amino_acid]
        codons = [codon for codon, _ in codons_with_weights]
        weights = [weight for _, weight in codons_with_weights]
        return random.choices(codons, weights=weights, k=1)[0]
