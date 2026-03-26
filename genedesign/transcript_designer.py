import logging
import math
import os
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

# Bump this string when algorithm changes would invalidate cached benchmark results.
# Leave it unchanged for tuning-only edits (e.g. max iterations, retry counts).
CACHE_VERSION = "2"


def _setup_logger() -> logging.Logger:
    """Creates a file logger writing to logs/<YYYY-MM-DD_HH-MM-SS>.log."""
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    log_path = (
        log_dir / f"{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}_{os.getpid()}.log"
    )

    logger = logging.getLogger(f"transcript_designer_{os.getpid()}")
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
        MAX_RESTARTS = 50
        # Minimum context so hairpin chunks (50 bp) always fit; extra
        # context lets the scorer catch hairpins between distant positions.
        MIN_CONTEXT = 100

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

        consecutive_same_hairpin = 0
        last_hairpin_sig = None
        EARLY_STOP = 15  # give up if same hairpin pattern repeats this many times

        # Fail-fast: track how many checks pass each restart.  If we haven't
        # improved in PATIENCE restarts, the gene is likely intractable — bail.
        best_checks_passed = 0
        restarts_since_improvement = 0
        PATIENCE = 30

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
                preamble_tail = preamble[-15:] if len(preamble) >= 15 else preamble
                locked_set = set(locked)
                for combo in self._enumerate_codon_combos(amino_window):
                    combo_seq = "".join(combo)
                    # Cheap pre-filter: skip combos that create 2+ hairpins
                    # at the boundary before running the expensive scorer.
                    if self._boundary_hairpin_count(preamble_tail, combo_seq) >= 2:
                        scored.append((combo, float("-inf")))
                        continue
                    check_seq = preamble + combo_seq + downstream
                    score = self._score_window(
                        check_seq, combo, len(preamble), locked_set
                    )
                    scored.append((combo, score))

                best_score = max(s for _, s in scored)

                # Compute the best possible pure codon score (no penalties).
                # If the actual best is significantly below this, a penalty
                # was applied — backtrack to escape.
                best_codon_only = max(
                    sum(math.log(self.codonToWeight.get(c, 1e-9) + 1e-9) for c in combo)
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
                    preamble_tail = preamble[-15:] if len(preamble) >= 15 else preamble
                    locked_set = set(locked)
                    for combo in self._enumerate_codon_combos(amino_window):
                        combo_seq = "".join(combo)
                        if self._boundary_hairpin_count(preamble_tail, combo_seq) >= 2:
                            scored.append((combo, float("-inf")))
                            continue
                        check_seq = preamble + combo_seq + downstream
                        score = self._score_window(
                            check_seq, combo, len(preamble), locked_set
                        )
                        scored.append((combo, score))

                    best_score = max(s for _, s in scored)

                if best_score <= float("-inf"):
                    zero_valid_windows += 1
                    chosen = max(scored, key=lambda x: x[1])[0]
                else:
                    near_best = [
                        c for c, s in scored if s >= best_score - NEAR_BEST_MARGIN
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

            # ── Targeted hairpin repair ────────────────────────────────────
            locked, repair_ok = self._repair_hairpins(locked, utr, peptide)
            if repair_ok:
                self.log.info("  Hairpin repair: passed (or not needed)")
            else:
                self.log.info("  Hairpin repair: could not fix all hairpins")

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
            passed_rbs, _ = self.internalRBSChecker.run(cds)

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

            checks_passed = sum(
                [
                    passed_promoter,
                    passed_hairpin,
                    passed_forbidden,
                    passed_rbs,
                ]
            )

            if checks_passed > best_checks_passed:
                best_checks_passed = checks_passed
                restarts_since_improvement = 0
            else:
                restarts_since_improvement += 1

            if checks_passed == 4:
                self.log.info(
                    "SUCCESS | restarts_used=%d",
                    restart_idx + 1,
                )
                return Transcript(selected_rbs, peptide, locked)

            # Fail-fast: no improvement in PATIENCE restarts → bail.
            if restarts_since_improvement >= PATIENCE:
                self.log.warning(
                    "EARLY STOP after %d restarts — stuck at %d/4 checks for %d restarts",
                    restart_idx + 1,
                    best_checks_passed,
                    PATIENCE,
                )
                break

            # Early stop: if the same hairpin signature repeats, the region
            # is likely inherently unfixable (e.g. repeated ATG-only codons).
            hairpin_sig = hairpin_hit or ""
            if hairpin_sig == last_hairpin_sig and hairpin_sig:
                consecutive_same_hairpin += 1
            else:
                consecutive_same_hairpin = 0
                last_hairpin_sig = hairpin_sig

            if consecutive_same_hairpin >= EARLY_STOP:
                self.log.warning(
                    "EARLY STOP after %d restarts — same hairpin pattern repeated %d times: %s",
                    restart_idx + 1,
                    EARLY_STOP,
                    hairpin_sig.replace("\n", " "),
                )
                break

            self.log.info(
                "  %d/4 checks passed (best=%d) — restarting.",
                checks_passed,
                best_checks_passed,
            )

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

    def _score_window(
        self,
        check_seq: str,
        combo: tuple,
        window_offset: int,
        locked_codons: set | None = None,
    ) -> float:
        """
        Scores *combo* in context.  Forbidden sequences → -inf.
        Hairpins and promoters → soft penalties.  Rare codons penalised.
        Novel codons (not yet in locked_codons) get a small diversity bonus.
        """
        from genedesign.seq_utils.hairpin_counter import hairpin_counter
        from genedesign.seq_utils.reverse_complement import reverse_complement

        window_end = window_offset + len(combo) * 3
        CHUNK, STEP = 50, 25
        HAIRPIN_PENALTY = 30.0
        PROMOTER_PENALTY = 30.0
        RARE_THRESHOLD = 0.1
        RARE_PENALTY = 5.0
        DIVERSITY_BONUS = 0.5
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
        rare_count = sum(
            1 for c in combo if self.codonToWeight.get(c, 0) < RARE_THRESHOLD
        )
        novel = sum(1 for c in combo if c not in locked_codons) if locked_codons else 0
        return (
            codon_score
            - HAIRPIN_PENALTY * hairpin_violations
            - PROMOTER_PENALTY * promoter_violations
            - RARE_PENALTY * rare_count
            + DIVERSITY_BONUS * novel
        )

    def _repair_hairpins(
        self, locked: list[str], utr: str, peptide: str, max_passes: int = 20
    ) -> tuple[list[str], bool]:
        """
        Targeted post-assembly repair.  Finds 50 bp chunks with >1 hairpin,
        identifies the exact stem pairs, and swaps the cheapest single codon
        synonym to break each stem — only if the swap doesn't increase the
        local hairpin count.  Tracks tried swaps to avoid oscillation.
        Returns (locked, success).
        """
        from genedesign.seq_utils.hairpin_counter import hairpin_counter

        _RC = {"A": "T", "T": "A", "C": "G", "G": "C"}
        CHUNK, STEP = 50, 25
        utr_len = len(utr)
        tried: set[tuple[int, str]] = set()  # (codon_idx, codon) already used
        codon_touch_count: dict[
            int, int
        ] = {}  # how many times each codon idx was swapped

        def _local_hairpin_count(seq: str, region_start: int, region_end: int) -> int:
            """Count hairpin violations in all 50bp chunks overlapping a region."""
            total = 0
            for cs in range(
                max(0, region_start - CHUNK + 1),
                min(len(seq) - CHUNK + 1, region_end),
                STEP,
            ):
                count, _ = hairpin_counter(seq[cs : cs + CHUNK], 3, 4, 9)
                if count > 1:
                    total += count - 1
            return total

        for pass_idx in range(max_passes):
            cds = "".join(locked) + "TAA"
            full_seq = utr + cds

            passed, _ = hairpin_checker(full_seq)
            if passed:
                return locked, True

            # Collect all hairpin stem pairs in failing chunks.
            hairpin_stems: list[tuple[int, int]] = []
            for chunk_start in range(0, len(full_seq) - CHUNK + 1, STEP):
                chunk = full_seq[chunk_start : chunk_start + CHUNK]
                count, _ = hairpin_counter(chunk, 3, 4, 9)
                if count <= 1:
                    continue
                for i in range(len(chunk) - 2):
                    s1 = chunk[i : i + 3]
                    for j in range(i + 7, min(i + 13, len(chunk) - 2)):
                        s2 = chunk[j : j + 3]
                        if (
                            s1[0] == _RC.get(s2[2])
                            and s1[1] == _RC.get(s2[1])
                            and s1[2] == _RC.get(s2[0])
                        ):
                            hairpin_stems.append((chunk_start + i, chunk_start + j))

            if not hairpin_stems:
                break

            # Log failing chunk positions for visibility.
            failing_positions = sorted({s for pair in hairpin_stems for s in pair})
            self.log.info(
                "    Repair pass %d: %d hairpin stems at seq positions %s",
                pass_idx + 1,
                len(hairpin_stems),
                failing_positions[:10],
            )

            # Try to break one hairpin per pass (re-check after each).
            fixed_any = False
            for stem_i, stem_j in hairpin_stems:
                # Codon indices whose bases overlap either stem.
                candidates: set[int] = set()
                for pos in range(stem_i, stem_i + 3):
                    cds_pos = pos - utr_len
                    if 0 <= cds_pos < len(peptide) * 3:
                        candidates.add(cds_pos // 3)
                for pos in range(stem_j, stem_j + 3):
                    cds_pos = pos - utr_len
                    if 0 <= cds_pos < len(peptide) * 3:
                        candidates.add(cds_pos // 3)

                # Measure current local hairpin count around this stem.
                region_start = min(stem_i, stem_j)
                region_end = max(stem_i, stem_j) + 3
                current_local = _local_hairpin_count(full_seq, region_start, region_end)

                best_fix: tuple[int, str] | None = None
                best_weight = -1.0
                best_local = current_local

                for idx in candidates:
                    aa = peptide[idx]
                    original = locked[idx]
                    for codon, weight in self.aminoAcidToCodons.get(aa, []):
                        if codon == original:
                            continue
                        if (idx, codon) in tried:
                            continue
                        # Would this swap break the stem RC match?
                        trial = locked[:]
                        trial[idx] = codon
                        trial_cds = "".join(trial) + "TAA"
                        trial_full = utr + trial_cds
                        s1 = trial_full[stem_i : stem_i + 3]
                        s2 = trial_full[stem_j : stem_j + 3]
                        still_hairpin = (
                            len(s1) == 3
                            and len(s2) == 3
                            and s1[0] == _RC.get(s2[2])
                            and s1[1] == _RC.get(s2[1])
                            and s1[2] == _RC.get(s2[0])
                        )
                        if still_hairpin:
                            continue
                        # Check that this swap doesn't increase local hairpins.
                        new_local = _local_hairpin_count(
                            trial_full, region_start, region_end
                        )
                        if new_local > best_local:
                            continue
                        if new_local < best_local or weight > best_weight:
                            best_fix = (idx, codon)
                            best_weight = weight
                            best_local = new_local

                if best_fix:
                    tried.add((best_fix[0], best_fix[1]))
                    tried.add((best_fix[0], locked[best_fix[0]]))  # mark old codon too
                    locked[best_fix[0]] = best_fix[1]
                    codon_touch_count[best_fix[0]] = (
                        codon_touch_count.get(best_fix[0], 0) + 1
                    )
                    touches = codon_touch_count[best_fix[0]]
                    self.log.info(
                        "      swap codon %d %s→%s (w=%.3f, local %d→%d, touches=%d)",
                        best_fix[0],
                        peptide[best_fix[0]],
                        best_fix[1],
                        best_weight,
                        current_local,
                        best_local,
                        touches,
                    )
                    fixed_any = True
                    break  # Re-check full sequence after each swap

            if not fixed_any:
                # ── Fallback: 2-codon swap for overlapping hairpins ───────
                # When single swaps can't help (e.g. GCA…TGC…GCA where both
                # hairpins share bases), try changing two codons at once.
                for chunk_start in range(0, len(full_seq) - CHUNK + 1, STEP):
                    chunk = full_seq[chunk_start : chunk_start + CHUNK]
                    count, _ = hairpin_counter(chunk, 3, 4, 9)
                    if count <= 1:
                        continue
                    # Codon indices overlapping this chunk.
                    codon_idxs: list[int] = []
                    for pos in range(chunk_start, chunk_start + CHUNK):
                        cds_pos = pos - utr_len
                        if 0 <= cds_pos < len(peptide) * 3:
                            ci = cds_pos // 3
                            if ci not in codon_idxs:
                                codon_idxs.append(ci)

                    current_local = _local_hairpin_count(
                        full_seq, chunk_start, chunk_start + CHUNK
                    )
                    best_pair = None
                    best_pair_weight = -1.0
                    best_pair_local = current_local

                    for ai, idx1 in enumerate(codon_idxs):
                        aa1 = peptide[idx1]
                        for c1, w1 in self.aminoAcidToCodons.get(aa1, []):
                            if c1 == locked[idx1]:
                                continue
                            for idx2 in codon_idxs[ai + 1 :]:
                                aa2 = peptide[idx2]
                                for c2, w2 in self.aminoAcidToCodons.get(aa2, []):
                                    if c2 == locked[idx2]:
                                        continue
                                    trial = locked[:]
                                    trial[idx1] = c1
                                    trial[idx2] = c2
                                    trial_cds = "".join(trial) + "TAA"
                                    trial_full = utr + trial_cds
                                    tv = _local_hairpin_count(
                                        trial_full, chunk_start, chunk_start + CHUNK
                                    )
                                    if tv < best_pair_local or (
                                        tv == best_pair_local
                                        and w1 + w2 > best_pair_weight
                                    ):
                                        best_pair = (idx1, c1, idx2, c2)
                                        best_pair_weight = w1 + w2
                                        best_pair_local = tv

                    if best_pair and best_pair_local < current_local:
                        i1, c1, i2, c2 = best_pair
                        locked[i1] = c1
                        locked[i2] = c2
                        codon_touch_count[i1] = codon_touch_count.get(i1, 0) + 1
                        codon_touch_count[i2] = codon_touch_count.get(i2, 0) + 1
                        self.log.info(
                            "      pair swap codons %d,%d → %s,%s (v %d→%d)",
                            i1,
                            i2,
                            c1,
                            c2,
                            current_local,
                            best_pair_local,
                        )
                        fixed_any = True
                        break  # re-check after pair swap

                if not fixed_any:
                    break

        # Summary of repair activity.
        if codon_touch_count:
            top_touched = sorted(
                codon_touch_count.items(), key=lambda x: x[1], reverse=True
            )[:5]
            self.log.info(
                "    Repair summary: %d passes, %d codons touched, top=%s",
                pass_idx + 1 if "pass_idx" in dir() else 0,
                len(codon_touch_count),
                [(idx, peptide[idx], cnt) for idx, cnt in top_touched],
            )

        cds = "".join(locked) + "TAA"
        full_seq = utr + cds
        passed, _ = hairpin_checker(full_seq)
        return locked, passed

    @staticmethod
    def _boundary_hairpin_count(preamble_tail: str, combo_seq: str) -> int:
        """
        Cheap hairpin count at the preamble–combo boundary.

        Checks whether any 3-mer in the combined region has its reverse
        complement 7-12 positions downstream (matching hairpin_counter's
        min_stem=3, min_loop=4, max_loop=9).  Only counts hairpins where
        at least one stem overlaps the combo region.
        """
        _RC = {"A": "T", "T": "A", "C": "G", "G": "C"}
        combined = preamble_tail + combo_seq
        combo_start = len(preamble_tail)
        count = 0
        clen = len(combined)
        for i in range(clen - 2):
            s1 = combined[i : i + 3]
            for j in range(i + 7, min(i + 13, clen - 2)):
                # At least one stem must touch the combo region.
                if j + 3 <= combo_start and i + 3 <= combo_start:
                    continue
                s2 = combined[j : j + 3]
                if (
                    s1[0] == _RC.get(s2[2])
                    and s1[1] == _RC.get(s2[1])
                    and s1[2] == _RC.get(s2[0])
                ):
                    count += 1
        return count

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
