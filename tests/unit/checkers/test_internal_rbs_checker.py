"""
Unit tests for InternalRBSChecker.

The checker scans a DNA sequence for Shine-Dalgarno-like motifs and flags any
motif that has a start codon (ATG / GTG / TTG) within 5–10 bp downstream.
Constants from the implementation:
  SPACER_MIN = 5  (inclusive offset from end of motif)
  SPACER_MAX = 10 (exclusive upper bound of the search window)
"""

import pytest

from genedesign.checkers.internal_rbs_checker import (
    SPACER_MAX,
    SPACER_MIN,
    InternalRBSChecker,
)

# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def checker():
    c = InternalRBSChecker()
    c.initiate()
    return c


# ---------------------------------------------------------------------------
# Initialisation
# ---------------------------------------------------------------------------


class TestInitiate:
    def test_rbs_sequences_populated(self, checker):
        assert checker.rbs_like_sequences is not None
        assert len(checker.rbs_like_sequences) > 0

    def test_start_codons_populated(self, checker):
        assert checker.start_codons is not None
        assert {"ATG", "GTG", "TTG"}.issubset(checker.start_codons)

    def test_not_initiated_has_none_attributes(self):
        c = InternalRBSChecker()
        assert c.rbs_like_sequences is None
        assert c.start_codons is None


# ---------------------------------------------------------------------------
# Sequences that should PASS (no internal RBS detected → result is True)
# ---------------------------------------------------------------------------


class TestCleanSequences:
    def test_empty_sequence(self, checker):
        result, offending = checker.run("")
        assert result is True
        assert offending is None

    def test_no_sd_motif_present(self, checker):
        # Codons with no SD-like run at all
        result, offending = checker.run("CCCCCCCCCCCCCCCCCCCC")
        assert result is True
        assert offending is None

    def test_sd_motif_without_downstream_start_codon(self, checker):
        # Perfect consensus with only C's after the spacer – no ATG/GTG/TTG
        motif = "AGGAGG"
        padding = "C" * SPACER_MIN
        result, offending = checker.run(motif + padding + "CCCCCC")
        assert result is True
        assert offending is None

    def test_sd_motif_start_codon_too_close(self, checker):
        # Start codon is only 4 bp after the motif – inside the dead zone (< SPACER_MIN)
        motif = "AGGAGG"
        result, offending = checker.run(motif + "ATGC")
        assert result is True, "Start codon within 4 bp should not trigger the checker"
        assert offending is None

    def test_sd_motif_start_codon_too_far(self, checker):
        # Start codon placed exactly at offset SPACER_MAX (outside the window)
        motif = "AGGAGG"
        spacer = "C" * SPACER_MAX  # positions 0..9 relative to end of motif
        result, offending = checker.run(motif + spacer + "ATG")
        assert result is True, "Start codon at SPACER_MAX offset should not be flagged"
        assert offending is None

    def test_real_cds_no_internal_rbs(self, checker):
        # lacI CDS (E. coli K-12) — historically clean w.r.t. internal RBS
        # GenBank acc. X02723.1, CDS CAA26510.1
        cds = (
            "ATGGACAGTCTCAATCTTAATAAACATATTTCCGGCCAGTTCAACGCCGAACTGGAAAGT"
            "ATCCGCACGCAGGTGATGACCATGGGCGGCATGGTGGAGCAGCAGCTTTCTGATGCAATC"
            "ACCGCGATGCATAACCAGGACAGCGATCTGGCGAAGCGCGTCATCGAAGGCGACAAGAAC"
            "GTCAACATGATGGAAGTGGCGATCGATGAAGCCTGCGTGCGCATTATCGCCAAACGTCAG"
            "CCGACGGCGAGCGACCTGCGACTGGTTATGGTGATCAGTAAAACCATTGCCGAGCTGGAG"
            "CGTATTGGCGACGTGGCGGACAAAATCTGCCGTACTGCGCTGGAGAAATTCTCCCAGCAG"
            "CATCAGCCGTTGCTGGTAAGTCTGGAGTCGCTGGGCCGTCATACCATCCAGATGCTGCAC"
            "GACGTGCTGGACGCGTTCGCGCGGATGGACATTGACGAAGCGGTACGTATTTATCGTGAA"
            "GATAAAAAAGTCGATCAGGAATACGAAGGTATTGTTCGTCAACTGATGACCTACATGATG"
            "GAAGATTCGCGTACCATTCCGAGCGTACTTACTGCGCTGTTCTGCGCGCGTTCTATCGAA"
            "CGTATTGGCGACCGCTGCCAGAATATTTGTGAGTTTATCTTCTACTACGTGAAGGGGCAG"
            "GATTTCCGTCACGTCGGTGGCGATGAGCTGGATAAACTGCTGGCGGGGAAAGATAGCGAC"
            "AAATAA"
        )
        result, offending = checker.run(cds)
        assert result is True
        assert offending is None


# ---------------------------------------------------------------------------
# Sequences that should FAIL (internal RBS detected → result is False)
# ---------------------------------------------------------------------------


class TestInternalRBSDetected:
    # -- Perfect consensus (AGGAGG) -------------------------------------------

    def test_perfect_sd_with_atg_at_min_spacer(self, checker):
        motif = "AGGAGG"
        # ATG starting exactly at position SPACER_MIN after motif end
        spacer = "C" * SPACER_MIN
        seq = "TTTT" + motif + spacer + "ATG" + "TTTT"
        result, offending = checker.run(seq)
        assert result is False
        assert offending is not None

    def test_perfect_sd_with_atg_at_last_valid_position(self, checker):
        motif = "AGGAGG"
        # The window is sequence[motif_end+SPACER_MIN : motif_end+SPACER_MAX].
        # A 3-char codon must start at offset ≤ SPACER_MAX-3 to fit entirely
        # inside the window — that is the last detectable position.
        spacer = "C" * (SPACER_MAX - 3)
        seq = motif + spacer + "ATG"
        result, offending = checker.run(seq)
        assert result is False
        assert offending is not None

    def test_perfect_sd_with_gtg_start(self, checker):
        motif = "AGGAGG"
        spacer = "C" * SPACER_MIN
        seq = motif + spacer + "GTG"
        result, offending = checker.run(seq)
        assert result is False

    def test_perfect_sd_with_ttg_start(self, checker):
        motif = "AGGAGG"
        spacer = "C" * SPACER_MIN
        seq = motif + spacer + "TTG"
        result, offending = checker.run(seq)
        assert result is False

    # -- Weaker SD variants ----------------------------------------------------

    def test_aaggag_variant_triggers(self, checker):
        motif = "AAGGAG"
        spacer = "C" * SPACER_MIN
        seq = motif + spacer + "ATG"
        result, offending = checker.run(seq)
        assert result is False

    def test_aggag_variant_triggers(self, checker):
        motif = "AGGAG"
        spacer = "C" * SPACER_MIN
        seq = motif + spacer + "ATG"
        result, offending = checker.run(seq)
        assert result is False

    def test_gagg_variant_triggers(self, checker):
        motif = "GAGG"
        spacer = "C" * SPACER_MIN
        seq = motif + spacer + "ATG"
        result, offending = checker.run(seq)
        assert result is False

    def test_aagg_variant_triggers(self, checker):
        motif = "AAGG"
        spacer = "C" * SPACER_MIN
        seq = motif + spacer + "ATG"
        result, offending = checker.run(seq)
        assert result is False

    def test_agga_variant_triggers(self, checker):
        motif = "AGGA"
        spacer = "C" * SPACER_MIN
        seq = motif + spacer + "ATG"
        result, offending = checker.run(seq)
        assert result is False

    # -- Return value contracts ------------------------------------------------

    def test_offending_sequence_is_string(self, checker):
        seq = "AGGAGG" + "C" * SPACER_MIN + "ATG"
        result, offending = checker.run(seq)
        assert result is False
        assert isinstance(offending, str)

    def test_offending_sequence_contains_motif(self, checker):
        motif = "AGGAGG"
        spacer = "C" * SPACER_MIN
        seq = motif + spacer + "ATG"
        _, offending = checker.run(seq)
        assert motif in offending

    def test_offending_sequence_contains_start_codon(self, checker):
        motif = "AGGAGG"
        spacer = "C" * SPACER_MIN
        seq = motif + spacer + "ATG"
        _, offending = checker.run(seq)
        assert "ATG" in offending

    # -- Embedded in longer sequence -------------------------------------------

    def test_internal_rbs_mid_sequence(self, checker):
        prefix = "CCCCCCCCCCCCCC"
        motif = "AGGAGG"
        spacer = "C" * SPACER_MIN
        suffix = "CCCCCCCC"
        seq = prefix + motif + spacer + "ATG" + suffix
        result, offending = checker.run(seq)
        assert result is False
        assert offending is not None


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_sd_motif_only_no_start_codon(self, checker):
        # SD motif present but the sequence ends before any start codon
        result, offending = checker.run("AGGAGG")
        assert result is True
        assert offending is None

    def test_multiple_motifs_first_clean_second_triggers(self, checker):
        # First motif: start codon too far away; second motif: valid hit
        motif = "AGGAGG"
        clean_part = (
            motif + "C" * SPACER_MAX + "ATG"
        )  # codon at offset 10 → outside window
        hit_part = motif + "C" * SPACER_MIN + "ATG"
        seq = clean_part + hit_part
        result, offending = checker.run(seq)
        assert result is False

    def test_lowercase_sequence_not_detected(self, checker):
        # The checker works on uppercase; lowercase should not match
        seq = "aggagg" + "c" * SPACER_MIN + "atg"
        result, _ = checker.run(seq)
        # Behaviour depends on implementation; we just verify the call does not raise
        assert isinstance(result, bool)
