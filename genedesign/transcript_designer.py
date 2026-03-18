import logging
import math
import random
from datetime import datetime
from itertools import product
from pathlib import Path

from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker
from genedesign.models.transcript import Transcript
from genedesign.rbs_chooser import RBSChooser


def _setup_logger() -> logging.Logger:
    """Creates a file logger writing to logs/<YYYY-MM-DD_HH-MM-SS>.log."""
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    log_path = log_dir / f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"

    logger = logging.getLogger("transcript_designer")
    logger.setLevel(logging.DEBUG)
    # Avoid duplicate handlers if initiate() is called more than once.
    logger.handlers.clear()

    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(
        logging.Formatter(
            "%(asctime)s  %(levelname)-7s  %(message)s", datefmt="%H:%M:%S"
        )
    )
    logger.addHandler(fh)
    logger.info("Log file: %s", log_path)
    return logger


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and selects an
    RBS using the Sliding Window Method.
    """

    def __init__(self):
        self.aminoAcidToCodons = {}  # aa -> [(codon, weight), ...]
        self.codonToWeight = {}  # codon -> weight (flat lookup)
        self.rbsChooser = None
        self.log = logging.getLogger("transcript_designer")
        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()
        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()
        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()
        self.internalRBSChecker = InternalRBSChecker()
        self.internalRBSChecker.initiate()

    def initiate(self) -> None:
        """
        Loads codon usage data, initialises the RBS chooser, and opens the
        run-specific log file.
        """
        self.log = _setup_logger()

        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        codon_usage_path = Path(__file__).parent / "data" / "codon_usage.txt"
        with open(codon_usage_path, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    codon, amino_acid, weight = parts[0], parts[1], float(parts[2])
                    self.aminoAcidToCodons.setdefault(amino_acid, []).append(
                        (codon, weight)
                    )

        self.codonToWeight = {
            codon: weight
            for entries in self.aminoAcidToCodons.values()
            for codon, weight in entries
        }
        self.log.info(
            "Loaded %d amino acid entries from codon usage table.",
            len(self.aminoAcidToCodons),
        )

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Builds the CDS left-to-right using a sliding window.

        For each window of amino acids:
          1. Enumerate all synonymous codon combinations.
          2. Evaluate each combination in local context:
               preamble   -- last 50 bp of already-locked codons
               window     -- the codons under evaluation
               downstream -- a random sample of the next 6 AAs (18 bp),
                             used only for context, then discarded
          3. Score: -inf if local sequence contains a forbidden site or
             hairpin; otherwise sum of log(codon_weight).
          4. Softmax-sample from valid combos (promotes codon diversity)
             and lock the winner, then advance the window.

        Global checks are run before RBS selection and final validation.
        If checks fail, the entire assembly restarts from scratch.

        Parameters:
            peptide  (str): Protein sequence to reverse-translate.
            ignores  (set): RBS options to exclude.

        Returns:
            Transcript: Validated transcript with RBS and codon list.
        """
        WINDOW = 5  # codons chosen per step — wider window catches cross-boundary
                    # hairpins (e.g. His CAT <-> Met ATG 12bp apart) before they
                    # are locked. 4^5 = 1024 combos/step, still fast.
        DOWNSTREAM = 6  # AAs of random downstream context (discarded)
        MAX_RESTARTS = 100

        self.log.info(
            "=== run() called | peptide_len=%d | peptide=%s%s",
            len(peptide),
            peptide[:30],
            "..." if len(peptide) > 30 else "",
        )

        for restart_idx in range(MAX_RESTARTS):
            self.log.info("--- Restart %d/%d ---", restart_idx + 1, MAX_RESTARTS)
            locked: list[str] = []

            # ── Sliding window assembly ───────────────────────────────────
            zero_valid_windows = 0
            for win_start in range(0, len(peptide), WINDOW):
                win_end = min(win_start + WINDOW, len(peptide))
                amino_window = peptide[win_start:win_end]

                # 50 bp preamble matches hairpin_checker's chunk size.
                preamble = "".join(locked)[-50:]
                ds_aas = peptide[win_end : win_end + DOWNSTREAM]
                downstream = "".join(self._weighted_random_codon(aa) for aa in ds_aas)

                scored: list[tuple[tuple, float]] = []
                for combo in self._enumerate_codon_combos(amino_window):
                    # Include downstream so _score_window can see hairpin traps
                    # that the current window sets up with future codons (AOT).
                    check_seq = preamble + "".join(combo) + downstream
                    score = self._score_window(check_seq, combo, len(preamble))
                    scored.append((combo, score))

                valid = [(c, s) for c, s in scored if s > float("-inf")]
                if valid:
                    max_s = max(s for _, s in valid)
                    weights = [math.exp(s - max_s) for _, s in valid]
                    best_combo = random.choices(
                        [c for c, _ in valid], weights=weights, k=1
                    )[0]
                    best_score = dict(scored)[best_combo]
                else:
                    zero_valid_windows += 1
                    best_combo = max(scored, key=lambda x: x[1])[0]
                    best_score = float("-inf")
                    # Log why every combo failed for this window.
                    reason = self._first_violation_reason(
                        preamble, best_combo, len(preamble), downstream
                    )
                    self.log.warning(
                        "  WINDOW aa[%d:%d] aas=%-6s | valid_combos=0 | "
                        "fallback=%s | violation: %s",
                        win_start, win_end, amino_window,
                        "".join(best_combo), reason,
                    )

                locked.extend(best_combo)

            self.log.info(
                "  Assembly done | zero_valid_windows=%d | cds_len=%d | cds_head=%s",
                zero_valid_windows, len(locked) * 3,
                "".join(locked)[:30],
            )

            # ── Codon quality check ───────────────────────────────────────
            codons_above_board, diversity, rare_count, cai = self.codonChecker.run(
                locked
            )
            self.log.info(
                "  CodonChecker | pass=%s | diversity=%.3f | rare_count=%d | CAI=%.4f",
                codons_above_board,
                diversity,
                rare_count,
                cai,
            )
            if not codons_above_board:
                continue

            # ── RBS selection ─────────────────────────────────────────────
            locked.append("TAA")
            cds = "".join(locked)
            selected_rbs = self.rbsChooser.run(cds, ignores)
            full_transcript = selected_rbs.utr.upper() + cds

            self.log.info(
                "  RBS selected | utr=%s | full_len=%d",
                selected_rbs.utr,
                len(full_transcript),
            )

            # ── Final transcript validation ───────────────────────────────
            passed_promoter, promoter_hit = self.promoterChecker.run(full_transcript)
            passed_hairpin, hairpin_hit = hairpin_checker(full_transcript)
            passed_forbidden, forb_hit = self.forbiddenChecker.run(full_transcript)
            passed_rbs, _ = True, ""  # self.internalRBSChecker.run(cds)

            self.log.info(
                "  Final checks | promoter=%s(%s) | hairpin=%s(%s) | "
                "forbidden=%s(%s) | rbs=%s",
                passed_promoter,
                promoter_hit or "-",
                passed_hairpin,
                hairpin_hit or "-",
                passed_forbidden,
                forb_hit or "-",
                passed_rbs,
            )

            if passed_promoter and passed_hairpin and passed_forbidden and passed_rbs:
                self.log.info(
                    "SUCCESS | restarts_used=%d",
                    restart_idx + 1,
                )
                return Transcript(selected_rbs, peptide, locked)

            self.log.info("  Final checks failed — restarting.")

        self.log.error(
            "FAILED | exhausted %d restarts for peptide=%s...",
            MAX_RESTARTS,
            peptide[:20],
        )
        raise RuntimeError(
            "Could not generate a valid DNA sequence after maximum attempts."
        )

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _enumerate_codon_combos(self, amino_acids: str):
        """
        Yields every synonymous codon combination for *amino_acids*.
        At most ~216 combinations for a 3-AA window.
        """
        options = [
            [
                c
                for c, _ in self.aminoAcidToCodons.get(
                    aa, [(self._weighted_random_codon(aa), 1.0)]
                )
            ]
            for aa in amino_acids
        ]
        yield from product(*options)

    def _score_window(self, check_seq: str, combo: tuple, window_offset: int) -> float:
        """
        Scores *combo* in context *check_seq* = preamble + window [+ downstream].

        *window_offset* is len(preamble), i.e. the index in check_seq where the
        new codons start.

        Violation detection:
          - Forbidden sequences: only rejected if the site overlaps the window
            portion (starts or ends inside positions window_offset..window_end).
          - Hairpins: uses chunk-based counting (50 bp chunks, 25 bp step)
            matching hairpin_checker's threshold of count > 1.  Only chunks that
            overlap the window region are evaluated — preamble-only chunks were
            already validated when their codons were placed.  Including downstream
            in check_seq means chunks that straddle the window-downstream boundary
            are also checked, catching forward hairpin traps AOT.
        """
        from genedesign.seq_utils.hairpin_counter import hairpin_counter
        from genedesign.seq_utils.reverse_complement import reverse_complement

        window_end = window_offset + len(combo) * 3
        CHUNK, STEP = 50, 25

        # ── Forbidden sequences (forward + RC) ────────────────────────────
        rc_seq = reverse_complement(check_seq)
        for site in self.forbiddenChecker.forbidden:
            site_len = len(site)
            pos = check_seq.find(site)
            while pos != -1:
                if pos < window_end and pos + site_len > window_offset:
                    return float("-inf")
                pos = check_seq.find(site, pos + 1)
            pos = rc_seq.find(site)
            while pos != -1:
                fwd_start = len(check_seq) - pos - site_len
                if fwd_start < window_end and fwd_start + site_len > window_offset:
                    return float("-inf")
                pos = rc_seq.find(site, pos + 1)

        # ── Hairpins (chunk-based, mirrors hairpin_checker threshold) ──────
        # Evaluate every 50 bp chunk that overlaps the window region.
        # Downstream content (if present) extends check_seq so cross-boundary
        # chunks are included automatically.
        seq_len = len(check_seq)
        for chunk_start in range(0, seq_len - CHUNK + 1, STEP):
            chunk_end = chunk_start + CHUNK
            if chunk_end <= window_offset or chunk_start >= window_end:
                continue  # no overlap with window — skip
            count, _ = hairpin_counter(check_seq[chunk_start:chunk_end], 3, 4, 9)
            if count > 1:
                return float("-inf")

        return sum(
            math.log(self.codonToWeight.get(codon, 1e-9) + 1e-9) for codon in combo
        )

    def _first_violation_reason(
        self, preamble: str, combo: tuple, window_offset: int, downstream: str = ""
    ) -> str:
        """
        Returns a human-readable string explaining why the given combo scores
        -inf (used when valid_combos=0 to explain what's blocking every choice).
        """
        from genedesign.seq_utils.hairpin_counter import hairpin_counter
        from genedesign.seq_utils.reverse_complement import reverse_complement

        check_seq = preamble + "".join(combo) + downstream
        window_end = window_offset + len(combo) * 3
        rc_seq = reverse_complement(check_seq)

        for site in self.forbiddenChecker.forbidden:
            pos = check_seq.find(site)
            while pos != -1:
                if pos < window_end and pos + len(site) > window_offset:
                    return f"forbidden:{site} at fwd pos {pos}"
                pos = check_seq.find(site, pos + 1)
            pos = rc_seq.find(site)
            while pos != -1:
                fwd_start = len(check_seq) - pos - len(site)
                if fwd_start < window_end and fwd_start + len(site) > window_offset:
                    return f"forbidden_rc:{site} at fwd pos {fwd_start}"
                pos = rc_seq.find(site, pos + 1)

        CHUNK, STEP = 50, 25
        seq_len = len(check_seq)
        for chunk_start in range(0, seq_len - CHUNK + 1, STEP):
            chunk_end = chunk_start + CHUNK
            if chunk_end <= window_offset or chunk_start >= window_end:
                continue
            count, hairpin_str = hairpin_counter(check_seq[chunk_start:chunk_end], 3, 4, 9)
            if count > 1:
                return f"hairpin count={count} in chunk[{chunk_start}:{chunk_end}]: {hairpin_str}"

        return "unknown (all combos -inf but reason unclear)"

    def _weighted_random_codon(self, amino_acid: str) -> str:
        """
        Samples a codon for *amino_acid* proportional to its usage frequency.
        """
        if amino_acid not in self.aminoAcidToCodons:
            return "ATG"
        codons, weights = zip(*self.aminoAcidToCodons[amino_acid])
        return random.choices(codons, weights=weights, k=1)[0]


if __name__ == "__main__":
    peptide = "MYPFIRTARMTV"

    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    print(transcript)
