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
        Builds the CDS left-to-right using a greedy sliding window.

        For each 3-codon (9 bp) window:
          1. Enumerate all synonymous codon combinations (~216 max).
          2. Score each in context: preamble (9 bp already locked) +
             window (9 bp being optimised) + downstream (18 bp random).
             Forbidden sequences are hard-rejected; hairpins and promoters
             are soft-penalised.
          3. Pick the best-scoring combo (greedy) and lock it.

        Restarts vary because downstream context is randomly sampled each time.

        Parameters:
            peptide  (str): Protein sequence to reverse-translate.
            ignores  (set): RBS options to exclude.

        Returns:
            Transcript: Validated transcript with RBS and codon list.
        """
        WINDOW = 3  # 3 codons (9 bp) per step
        DOWNSTREAM = 6  # 6 codons (18 bp) of random downstream context
        MAX_RESTARTS = 500
        # Minimum context length so hairpin chunks (50 bp) always fit.
        MIN_CONTEXT = 50

        random.seed(42)

        self.log.info(
            "=== run() called | peptide_len=%d | peptide=%s%s",
            len(peptide),
            peptide[:30],
            "..." if len(peptide) > 30 else "",
        )

        # Pre-compute the UTR so the first windows see realistic context.
        utr = ""
        for rbsopt in self.rbsChooser.rbsOptions:
            if rbsopt not in ignores:
                utr = rbsopt.utr.upper()
                break

        for restart_idx in range(MAX_RESTARTS):
            self.log.info("--- Restart %d/%d ---", restart_idx + 1, MAX_RESTARTS)
            locked: list[str] = []

            # ── Sliding window assembly with backtracking ────────────────
            # Advance by WINDOW codons.  If every combo in a window is
            # penalised (best_score includes hairpin/promoter penalties),
            # backtrack by removing the previous BACKTRACK codons and
            # re-evaluate the combined region with a wider window.
            BACKTRACK = 2
            NEAR_BEST_MARGIN = 1.0
            zero_valid_windows = 0
            backtrack_used = set()  # positions already backtracked from
            win_start = 0
            while win_start < len(peptide):
                win_end = min(win_start + WINDOW, len(peptide))
                amino_window = peptide[win_start:win_end]

                full_prefix = utr + "".join(locked)
                preamble = (
                    full_prefix[-MIN_CONTEXT:]
                    if len(full_prefix) > MIN_CONTEXT
                    else full_prefix
                )

                ds_aas = peptide[win_end : win_end + DOWNSTREAM]
                downstream = "".join(self._weighted_random_codon(aa) for aa in ds_aas)
                while (
                    len(preamble) + len(amino_window) * 3 + len(downstream)
                    < MIN_CONTEXT + len(amino_window) * 3
                ):
                    downstream += random.choice("ACGT")

                scored: list[tuple[tuple, float]] = []
                for combo in self._enumerate_codon_combos(amino_window):
                    check_seq = preamble + "".join(combo) + downstream
                    score = self._score_window(check_seq, combo, len(preamble))
                    scored.append((combo, score))

                best_score = max(s for _, s in scored)

                # Compute the best possible pure codon score (no penalties).
                # If the actual best is significantly below this, a penalty
                # was applied — backtrack to escape.
                best_codon_only = max(
                    sum(
                        math.log(self.codonToWeight.get(c, 1e-9) + 1e-9)
                        for c in combo
                    )
                    for combo, _ in scored
                )
                has_penalties = best_score < best_codon_only - 10.0

                if has_penalties and win_start > 0 and win_start not in backtrack_used:
                    backtrack_used.add(win_start)
                    bt = min(BACKTRACK, win_start)
                    locked = locked[:-bt]
                    win_start -= bt
                    # Retry with the wider window (previous codons + current).
                    win_end = min(win_start + WINDOW + bt, len(peptide))
                    amino_window = peptide[win_start:win_end]

                    full_prefix = utr + "".join(locked)
                    preamble = (
                        full_prefix[-MIN_CONTEXT:]
                        if len(full_prefix) > MIN_CONTEXT
                        else full_prefix
                    )
                    ds_aas = peptide[win_end : win_end + DOWNSTREAM]
                    downstream = "".join(
                        self._weighted_random_codon(aa) for aa in ds_aas
                    )
                    while (
                        len(preamble) + len(amino_window) * 3 + len(downstream)
                        < MIN_CONTEXT + len(amino_window) * 3
                    ):
                        downstream += random.choice("ACGT")

                    scored = []
                    for combo in self._enumerate_codon_combos(amino_window):
                        check_seq = preamble + "".join(combo) + downstream
                        score = self._score_window(
                            check_seq, combo, len(preamble)
                        )
                        scored.append((combo, score))

                    best_score = max(s for _, s in scored)

                if best_score <= float("-inf"):
                    zero_valid_windows += 1
                    chosen = max(scored, key=lambda x: x[1])[0]
                else:
                    near_best = [
                        c
                        for c, s in scored
                        if s >= best_score - NEAR_BEST_MARGIN
                    ]
                    chosen = random.choice(near_best)

                locked.extend(chosen)
                win_start += len(chosen)

            self.log.info(
                "  Assembly done | zero_valid_windows=%d | cds_len=%d | cds_head=%s",
                zero_valid_windows,
                len(locked) * 3,
                "".join(locked)[:30],
            )

            # ── Codon quality check (logged, not gated) ──────────────────
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
        Scores *combo* in context.  Forbidden sequences → -inf.
        Hairpins and promoters → large penalties.  Otherwise codon weight sum.
        """
        from genedesign.seq_utils.hairpin_counter import hairpin_counter
        from genedesign.seq_utils.reverse_complement import reverse_complement

        window_end = window_offset + len(combo) * 3
        CHUNK, STEP = 50, 25
        HAIRPIN_PENALTY = 200.0
        PROMOTER_PENALTY = 200.0
        PROMOTER_FRAME = 29
        PROMOTER_THRESHOLD = 9.134

        # ── Forbidden sequences (forward + RC) — hard reject ────────────
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

        # ── Hairpins (soft penalty) ──────────────────────────────────────
        # Check ALL chunks that overlap the window — the guaranteed minimum
        # context ensures at least one full 50 bp chunk always exists.
        hairpin_violations = 0
        seq_len = len(check_seq)
        for chunk_start in range(0, seq_len - CHUNK + 1, STEP):
            chunk_end = chunk_start + CHUNK
            if chunk_end <= window_offset or chunk_start >= window_end:
                continue
            count, _ = hairpin_counter(check_seq[chunk_start:chunk_end], 3, 4, 9)
            if count > 1:
                hairpin_violations += count - 1

        # ── Promoter PWM (soft penalty) ──────────────────────────────────
        promoter_violations = 0
        base_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
        combined = check_seq + "x" + rc_seq
        for i in range(len(combined) - PROMOTER_FRAME + 1):
            frame_end = i + PROMOTER_FRAME
            if i < len(check_seq):
                if frame_end <= window_offset or i >= window_end:
                    continue
            else:
                rc_i = i - len(check_seq) - 1
                if rc_i < 0:
                    continue
                fwd_start = len(check_seq) - rc_i - PROMOTER_FRAME
                if (
                    fwd_start + PROMOTER_FRAME <= window_offset
                    or fwd_start >= window_end
                ):
                    continue
            partseq = combined[i:frame_end]
            score = 0.0
            for x, base in enumerate(partseq):
                y = base_idx.get(base, -1)
                if y != -1:
                    score += self.promoterChecker.pwm[y][x]
            if score >= PROMOTER_THRESHOLD:
                promoter_violations += 1

        codon_score = sum(
            math.log(self.codonToWeight.get(codon, 1e-9) + 1e-9) for codon in combo
        )
        return (
            codon_score
            - HAIRPIN_PENALTY * hairpin_violations
            - PROMOTER_PENALTY * promoter_violations
        )

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
