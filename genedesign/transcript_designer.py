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
        Builds the CDS left-to-right using a sliding 3-codon window.

        For each window of 3 amino acids:
          1. Enumerate all synonymous codon combinations (~<=216 for 3 AAs).
          2. Evaluate each combination in local context:
               preamble   -- last 9 bp of already-locked codons
               window     -- the 3 codons under evaluation (9 bp)
               downstream -- a random sample of the next 6 AAs (18 bp),
                             used only for context, then discarded
          3. Score: -inf if local sequence contains a forbidden site or
             hairpin; otherwise sum of log(codon_weight).
          4. Softmax-sample from valid combos (promotes codon diversity)
             and lock the winner, then advance by 3.

        After assembly a repair pass fixes boundary violations, then global
        checks are run before RBS selection and final validation.

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
        MAX_REPAIRS = 2000

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

                # 50 bp preamble matches hairpin_checker / _find_failing_region.
                preamble = "".join(locked)[-50:]
                ds_aas = peptide[win_end : win_end + DOWNSTREAM]
                downstream = "".join(self._weighted_random_codon(aa) for aa in ds_aas)

                scored: list[tuple[tuple, float]] = []
                for combo in self._enumerate_codon_combos(amino_window):
                    check_seq = preamble + "".join(combo)
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
                        preamble, best_combo, len(preamble)
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

            # ── Repair pass ───────────────────────────────────────────────
            # When enumeration produces no change the violation is caused by
            # flanking locked context that the current window can't escape.
            # We widen by 1 codon each side per STUCK_THRESHOLD stuck reps,
            # but cap at MAX_EXTRA to avoid exponential enumeration blowup
            # (6^n combos per repair attempt).  If widening to the cap still
            # produces no change, the assembly is irresolvable — restart now
            # rather than burning the full repair budget.
            STUCK_THRESHOLD = 3  # stuck reps before widening
            MAX_EXTRA = 2  # max extra codons each side (cap at 7 total)

            repair_count = 0
            stuck_at_pos = -1
            stuck_streak = 0
            irresolvable = False

            for repair_idx in range(MAX_REPAIRS):
                cds = "".join(locked)
                fail_result = self._find_failing_region(cds)
                if fail_result is None:
                    break
                fail_pos, fail_type, fail_detail = fail_result

                codon_idx = fail_pos // 3

                if fail_pos == stuck_at_pos:
                    stuck_streak += 1
                else:
                    stuck_at_pos = fail_pos
                    stuck_streak = 0

                extra = min(stuck_streak // STUCK_THRESHOLD, MAX_EXTRA)
                lo = max(0, codon_idx - 1 - extra)
                hi = min(len(locked), codon_idx + 2 + extra)

                before = "".join(locked[lo:hi])
                self._repair_region(locked, peptide, lo, hi)
                after = "".join(locked[lo:hi])
                repair_count += 1

                changed = before != after
                self.log.debug(
                    "  REPAIR #%d | %s:%s | pos=%d codon=%d | "
                    "stuck=%d extra=%d | aa[%d:%d] | %s->%s | %s | "
                    "ctx=...%s[%s]%s...",
                    repair_count, fail_type, fail_detail,
                    fail_pos, codon_idx,
                    stuck_streak, extra, lo, hi,
                    before, after,
                    "CHANGED" if changed else "STUCK",
                    cds[max(0, lo * 3 - 6) : lo * 3],
                    cds[lo * 3 : hi * 3],
                    cds[hi * 3 : hi * 3 + 6],
                )

                if not changed and extra >= MAX_EXTRA:
                    self.log.warning(
                        "  Irresolvable | %s:%s | pos=%d codon=%d | "
                        "after %d repair(s) at max window — restarting.",
                        fail_type, fail_detail, fail_pos, codon_idx, repair_count,
                    )
                    irresolvable = True
                    break
            else:
                self.log.warning(
                    "  Repair budget exhausted after %d attempts — restarting.",
                    MAX_REPAIRS,
                )
                irresolvable = True

            if irresolvable:
                continue

            self.log.info("  Repair pass done. %d repair(s) applied.", repair_count)

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
                    "SUCCESS | restarts_used=%d | repairs_applied=%d",
                    restart_idx + 1,
                    repair_count,
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

        Violation detection is *differential*:
          - Forbidden sequences: only rejected if the site overlaps the window
            portion (starts or ends inside positions window_offset..end).  Sites
            entirely within the preamble/downstream are the responsibility of the
            codons that placed them and will be caught by earlier/later windows.
          - Hairpins: only rejected if adding the window codons *increases* the
            hairpin count vs the preamble alone.  This avoids blaming the current
            combo for hairpins that were already present in locked context.
        """
        from genedesign.seq_utils.hairpin_counter import hairpin_counter
        from genedesign.seq_utils.reverse_complement import reverse_complement

        window_end = window_offset + len(combo) * 3

        # ── Forbidden sequences (forward + RC) ────────────────────────────
        rc_seq = reverse_complement(check_seq)
        for site in self.forbiddenChecker.forbidden:
            site_len = len(site)
            # Forward strand
            pos = check_seq.find(site)
            while pos != -1:
                # Reject only if site overlaps the window region
                if pos < window_end and pos + site_len > window_offset:
                    return float("-inf")
                pos = check_seq.find(site, pos + 1)
            # Reverse complement
            pos = rc_seq.find(site)
            while pos != -1:
                fwd_start = len(check_seq) - pos - site_len
                if fwd_start < window_end and fwd_start + site_len > window_offset:
                    return float("-inf")
                pos = rc_seq.find(site, pos + 1)

        # ── Hairpins (differential) ────────────────────────────────────────
        preamble_count, _ = hairpin_counter(check_seq[:window_offset], 3, 4, 9)
        full_count, _     = hairpin_counter(check_seq, 3, 4, 9)
        if full_count > preamble_count:
            return float("-inf")

        return sum(
            math.log(self.codonToWeight.get(codon, 1e-9) + 1e-9) for codon in combo
        )

    def _repair_region(self, locked: list[str], peptide: str, lo: int, hi: int) -> None:
        """
        Re-enumerates all codon combos for positions lo..hi-1 using the
        already-locked sequence as context on both sides, and replaces those
        positions with the best-scoring combination found.
        """
        amino_window = peptide[lo:hi]
        preamble = "".join(locked[:lo])[-50:]
        downstream = "".join(locked[hi : hi + 18])

        best_combo: tuple | None = None
        best_score = float("-inf")

        for combo in self._enumerate_codon_combos(amino_window):
            check_seq = preamble + "".join(combo) + downstream
            score = self._score_window(check_seq, combo, len(preamble))
            if best_combo is None or score > best_score:
                best_score = score
                best_combo = combo

        for j, codon in enumerate(best_combo):
            locked[lo + j] = codon

    def _find_failing_region(self, dna: str) -> tuple[int, str, str] | None:
        """
        Returns (approx_pos, violation_type, detail) for the first violation
        found in *dna*, or None if the sequence is clean.

        violation_type is one of: "hairpin", "forbidden", "forbidden_rc", "promoter"
        detail is the offending sequence or site name.
        """
        from genedesign.seq_utils.hairpin_counter import hairpin_counter
        from genedesign.seq_utils.reverse_complement import reverse_complement

        chunk, step = 50, 25
        for i in range(0, len(dna) - chunk + 1, step):
            count, hairpin_str = hairpin_counter(dna[i : i + chunk], 3, 4, 9)
            if count > 1:
                return i + chunk // 2, "hairpin", hairpin_str or "?"

        rc = reverse_complement(dna)
        for site in self.forbiddenChecker.forbidden:
            idx = dna.find(site)
            if idx >= 0:
                return idx + len(site) // 2, "forbidden", site
            idx = rc.find(site)
            if idx >= 0:
                return len(dna) - idx - len(site) // 2, "forbidden_rc", site

        passed, seq = self.promoterChecker.run(dna)
        if not passed and seq:
            idx = dna.find(seq)
            if idx >= 0:
                return idx + len(seq) // 2, "promoter", seq[:20]
            idx = rc.find(seq)
            if idx >= 0:
                return len(dna) - idx - len(seq) // 2, "promoter_rc", seq[:20]

        return None

    def _first_violation_reason(self, preamble: str, combo: tuple, window_offset: int) -> str:
        """
        Returns a human-readable string explaining why the given combo scores
        -inf (used when valid_combos=0 to explain what's blocking every choice).
        """
        from genedesign.seq_utils.hairpin_counter import hairpin_counter
        from genedesign.seq_utils.reverse_complement import reverse_complement

        check_seq = preamble + "".join(combo)
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

        preamble_count, _ = hairpin_counter(check_seq[:window_offset], 3, 4, 9)
        full_count, hairpin_str = hairpin_counter(check_seq, 3, 4, 9)
        if full_count > preamble_count:
            return f"hairpin (preamble={preamble_count}→full={full_count}): {hairpin_str}"

        return "unknown (all combos -inf but reason unclear)"

        Returns:
            list: The best codon combination for this window.
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
    start_time = time.time()
    transcript = designer.run(peptide, ignores)
    print(transcript)
    print(f"\nCompleted in {elapsed_time:.2f} seconds")
