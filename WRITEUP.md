# TranscriptDesigner Bug Report

## Overview

The benchmarker (`tests/benchmarking/proteome_benchmarker.py`) processed all 2,458 genes without exceptions but produced 6,956 validation failures across five checkers. All failures traced back to two bugs in `genedesign/transcript_designer.py`.

---

## Bug 1: CodonChecker Was Never Called

**Failures caused:** 2,322 Codon Usage Checker failures

### What was wrong

`self.codonChecker` was correctly instantiated and initialized in `__init__`, but it was never invoked during sequence generation. The `_all_checkers_pass` method only called the promoter, hairpin, and forbidden sequence checkers. Every generated sequence was accepted regardless of codon quality.

### Why it mattered

The CodonChecker enforces three thresholds:
- Codon diversity ≥ 0.5 (fraction of unique codons)
- Rare codon count ≤ 3
- Codon Adaptation Index (CAI) ≥ 0.2

Because this check was skipped entirely, the designer could return sequences with very low diversity, too many rare codons, or a poor CAI — all of which the benchmarker's validation then caught and flagged.

### Fix

Added a `self.codonChecker.run(codons)` call at the top of the retry loop. If the codon quality check fails, the loop continues to the next candidate sequence.

---

## Bug 2: Only the CDS Was Validated, Not the Full Transcript (UTR + CDS)

**Failures caused:** 1,924 Forbidden Sequence failures, 2,326 Hairpin failures, 254 Promoter failures

### What was wrong

The benchmarker validates the full transcript — the RBS 5' UTR concatenated with the CDS — for hairpins, forbidden restriction sites, and internal promoters:

```python
transcript_dna = result["transcript"].rbs.utr.upper() + cds
```

The original designer only validated the bare CDS before RBS selection. The RBS was chosen *after* the CDS was accepted, so the combined UTR+CDS sequence was never checked. Hairpins, forbidden sequences, or promoter motifs that appeared at the UTR–CDS junction (or were introduced by the specific UTR chosen) went undetected.

### Why it mattered

RBS UTR sequences are typically 30–50 bp long. A 50 bp sliding window check means roughly the first 20 bp of the CDS is evaluated jointly with the entire UTR. Even a CDS that individually passes all checks can fail when prefixed with a particular UTR. Since the original code never validated this combination, the vast majority of transcripts contained one or more of these junction-level issues by the time the benchmarker checked them.

### Fix

Moved RBS selection inside the retry loop. After a candidate CDS passes the codon and sequence checks, the stop codon is appended, an RBS is selected, and then the full transcript (`selectedRBS.utr.upper() + cds`) is run through `_all_checkers_pass`. Only if the full transcript passes are the results returned. Otherwise, the loop continues and a new CDS candidate is generated.

---

## Summary of Changes

| Change | Failures addressed |
|---|---|
| Added `self.codonChecker.run(codons)` call in the retry loop | 2,322 Codon Usage |
| Moved RBS selection inside the retry loop | — |
| Added `_all_checkers_pass(selectedRBS.utr.upper() + cds)` check after RBS selection | 1,924 Forbidden Sequence, 2,326 Hairpin, 254 Promoter |

All changes were confined to the `run` method of `TranscriptDesigner` in `genedesign/transcript_designer.py`.

---

## Experiment Log: Iterating Toward a Working Designer

The bugs above described what was wrong with the *original* code. Fixing them required replacing the naïve retry loop with an entirely new generation strategy. Below is a chronological record of the approaches tried, what worked, what didn't, and why.

### Attempt 1: Remove Repair, Use Sliding Window Only

**Approach:** Replaced the original generate-then-repair loop with a greedy sliding-window assembler. Each 3-codon (9 bp) window enumerates all synonymous codon combinations, scores them in context (preamble + window + random downstream), and locks the best.

**Result:** Every generated sequence triggered hairpin failures. The greedy scorer had no awareness of hairpins or promoters — it only optimized codon weight.

### Attempt 2: Soft Penalty Scoring

**Approach:** Instead of hard-rejecting sequences during assembly, added soft penalties for hairpin and promoter violations in `_score_window`. Forbidden sequences remained a hard reject (`-inf`). Hairpins and promoters each incurred a `-200` penalty per violation, subtracted from the codon weight score.

**Result:** Stopped the cascading zero-valid-window failures from Attempt 1. CAI improved from ~0.35 to ~0.45. However, the first few windows had no hairpin checking because `check_seq` was shorter than 50 bp (no full chunk to evaluate).

### Attempt 3: UTR Preamble + Minimum Context

**Approach:** Pre-computed the UTR (`aaagaggagaaatactag`) and prepended it to the locked codons when building the preamble. Added `MIN_CONTEXT = 50` — if `preamble + window + downstream` was shorter than 50 bp + window length, padded downstream with random bases. This guaranteed at least one full 50 bp chunk for the hairpin checker in every window.

**Result:** First few windows now saw realistic context including the UTR. Promoter scoring was also integrated into `_score_window` using the PWM matrix. Short proteins (≤120 AAs) started passing consistently.

### Attempt 4: Near-Best Randomization

**Approach:** Instead of always picking the single best-scoring combo, collected all combos within `NEAR_BEST_MARGIN = 1.0` of the best score and chose randomly among them.

**Result:** Broke the deterministic trap where the greedy assembler always picked the same hairpin-forming codons (e.g., `TTC(AGATTATT)GAA` appeared in 248/260 windows). Different restarts now explored different parts of the solution space.

### Attempt 5: Backtracking

**Approach:** After scoring a window, compared `best_score` against the pure codon score (no penalties). If the best available combo had significant penalties (hairpin or promoter), backed up by `BACKTRACK = 2` codons and re-evaluated a wider combined window.

**Result:** Allowed the assembler to escape local traps where a previous codon choice forced all subsequent choices into hairpin territory. Combined with near-best randomization, the first two benchmark proteins passed within ~9 restarts.

### Attempt 6: Boundary Hairpin Pre-filter

**Approach:** Added `_boundary_hairpin_count`, a cheap O(n²) check on just 24 bp (last 15 bp of preamble + 9 bp of candidate combo). If ≥2 hairpin stems formed at the boundary, the combo was scored `-inf` without calling the expensive `_score_window`.

**Result:** Substantially faster assembly — many hairpin-forming combos were eliminated before the expensive forbidden-sequence, hairpin-chunk, and promoter-PWM scoring ran. The third benchmark protein (421 AAs) started passing.

### Attempt 7: Targeted Post-Assembly Repair (Single-Codon Swaps)

**Approach:** After assembly, ran the full hairpin checker on `UTR + CDS`. If it failed, identified the exact 3 bp stem pairs in failing 50 bp chunks, found which codons overlapped those stems, and tried swapping each to a synonym that broke the RC match — only accepting swaps that didn't increase local hairpin count.

**Result:** Fixed most hairpin failures on proteins up to ~500 AAs. However, the repair oscillated on some proteins: swapping codon 27 `D→GAC` broke hairpin A but created hairpin B, and swapping back broke B but recreated A.

### Attempt 8: Oscillation Prevention (Tried-Set + Net-Effect Check)

**Approach:** Tracked all `(codon_index, codon)` pairs already tried in a `tried` set. Before accepting a swap, checked that it didn't increase the total hairpin violation count across all overlapping 50 bp chunks (not just the target stem). Marked both the old and new codon as tried to prevent flip-flopping.

**Result:** Eliminated the `TAC↔TAT` oscillation pattern. Repair made monotonic progress. However, a new failure mode appeared: overlapping hairpins that shared stem bases (e.g., `GCA(ACCGA)TGC` and `TGC(CTAT)GCA`) couldn't be broken by any single-codon swap.

### Attempt 9: Two-Codon Pair Swaps (Fallback)

**Approach:** When single-codon swaps were exhausted, fell back to trying all pairs of codon swaps within each failing 50 bp chunk. Evaluated each pair against `_local_hairpin_count` (all overlapping chunks, not just the one chunk) and accepted only pairs that reduced violations.

**Result:** Handled the overlapping-hairpin case for most proteins. More proteins in the benchmark passed. But one protein — cdhC (524 AAs) — still thrashed.

### Attempt 10: Increased MIN_CONTEXT to 100

**Approach:** Doubled the minimum context window from 50 bp to 100 bp. This let `_score_window` see more of the locked sequence during assembly, catching hairpin formations between positions up to 100 bp apart (vs. 50 bp previously).

**Result:** Reduced the number of hairpins that made it through assembly, meaning repair had less work to do. Assembly was ~2× slower per window but required fewer restarts. cdhC still failed.

### Attempt 11: Diagnosing cdhC — An Inherently Unfixable Region

**Analysis:** The protein cdhC contains the stretch `QPMAMQPMPMQMP` (codons 465–477) — four Methionines interspersed with Prolines and Alanines. Methionine has only one codon (`ATG`), so the `ATG...ATG...ATG` pattern is locked. Exhaustive enumeration of all 4,096 synonymous codon combinations for this region showed that **every single one** produces at least 2 hairpins within a 15 bp span. Both hairpins always fall in the same 50 bp checker chunks (`[1400–1450]` and `[1425–1475]`), so they can never be distributed across chunk boundaries.

**Conclusion:** This region is unsolvable given the genetic code and the hairpin checker's parameters (3 bp stems, 4–9 bp loops, >1 per 50 bp chunk = failure). No synonymous codon substitution — single, pair, or exhaustive — can reduce the chunk hairpin count to ≤1.

### Attempt 12: Early Stop + Graceful Failure

**Approach:** Added early-stop detection: if the same hairpin signature (e.g., `GCA(ACCGA)TGC / TGC(CTAT)GCA`) repeats for 15 consecutive restarts, the assembler stops immediately instead of thrashing through all 500 restarts. The benchmarker was wrapped in a `try/except` so that a `RuntimeError` from the designer is logged and the run continues to the next protein.

**Result:** cdhC fails in ~16 restarts (~50 seconds) instead of ~500 restarts (~15+ minutes). The benchmarker completes the full proteome, logging cdhC as a failure without crashing.

---

## Final Architecture

---

## Detailed Algorithm Description

The final algorithm in `TranscriptDesigner.run()` has three phases: **assembly**, **repair**, and **validation**. If validation fails, the entire process restarts with different random downstream context. An early-stop mechanism detects inherently unfixable sequences and aborts without exhausting all restarts.

### Constants

| Name | Value | Purpose |
|---|---|---|
| `WINDOW` | 3 | Codons per sliding-window step (9 bp) |
| `DOWNSTREAM` | 6 | Codons of random downstream context (18 bp) |
| `MIN_CONTEXT` | 100 | Minimum bases of upstream context for scoring |
| `MAX_RESTARTS` | 500 | Maximum full-sequence rebuild attempts |
| `BACKTRACK` | 2 | Codons to undo when backtracking |
| `NEAR_BEST_MARGIN` | 1.0 | Score tolerance for near-best randomization |
| `HAIRPIN_PENALTY` | 200.0 | Score penalty per hairpin violation in a chunk |
| `PROMOTER_PENALTY` | 200.0 | Score penalty per promoter hit in a window |
| `EARLY_STOP` | 15 | Consecutive identical hairpin signatures before giving up |

### Phase 1: Greedy Sliding-Window Assembly

The CDS is built left-to-right, 3 codons (9 bp) at a time. The random seed is set to 42 at the start of `run()` so that restarts diverge only because downstream context is randomly sampled.

**For each window position:**

1. **Build the context string.** The scorer needs a realistic DNA context around the candidate codons. The preamble is the last `MIN_CONTEXT` (100) bases of `UTR + locked_codons`. Using the UTR as prefix for the first windows ensures that hairpins and promoters spanning the UTR-CDS junction are detected. Downstream context is generated by sampling codons for the next 6 amino acids using frequency-weighted random selection, padded with random bases if the total context is too short.

2. **Enumerate all synonymous codon combinations.** For a 3-amino-acid window, this is the Cartesian product of each amino acid's synonymous codons (typically 1-6 each, so at most ~216 combos).

3. **Boundary hairpin pre-filter.** Before running the expensive scorer, `_boundary_hairpin_count` checks a 24 bp region (last 15 bp of preamble + 9 bp candidate) for hairpin stem pairs (3 bp stems whose reverse complement appears 7-12 positions away). If 2+ hairpins form at this boundary, the combo is scored `-inf` without further evaluation. This eliminates a large fraction of candidates cheaply.

4. **Score each combo in context** via `_score_window(check_seq, combo, window_offset)`:
   - **Forbidden sequences (hard reject):** Scans the check sequence and its reverse complement for any forbidden restriction site. If a site overlaps the window being optimized, the combo scores `-inf`.
   - **Hairpin violations (soft penalty):** Slides 50 bp chunks (step 25) across the check sequence. For each chunk overlapping the window, counts hairpins using `hairpin_counter(chunk, min_stem=3, min_loop=4, max_loop=9)`. If a chunk has >1 hairpin, the excess count is multiplied by `HAIRPIN_PENALTY` (200) and subtracted from the score.
   - **Promoter violations (soft penalty):** Slides a 29 bp frame across the check sequence and its reverse complement. For each frame overlapping the window, computes the position weight matrix (PWM) score. If the score exceeds the threshold (9.134), `PROMOTER_PENALTY` (200) is subtracted.
   - **Codon weight score:** The base score is the sum of `log(weight)` for each codon in the combo, where weight comes from the *E. coli* codon usage table.
   - **Final score** = codon_weight_score - (hairpin_penalty * hairpin_violations) - (promoter_penalty * promoter_violations).

5. **Backtracking.** After scoring, the algorithm compares the best achievable score against the best pure codon score (no penalties). If the gap exceeds 10.0 — meaning every combo incurs significant penalties — and the current position hasn't been backtracked before, the algorithm removes the previous `BACKTRACK` (2) codons and re-evaluates a wider combined window. This lets it escape situations where a previous codon choice forced all subsequent choices into hairpin or promoter territory.

6. **Near-best randomization.** Instead of always picking the single highest-scoring combo, the algorithm collects all combos within `NEAR_BEST_MARGIN` (1.0) of the best score and chooses uniformly at random. This prevents the greedy assembler from deterministically selecting the same hairpin-forming codons on every restart.

7. **Lock and advance.** The chosen codons are appended to the `locked` list and the window advances by the number of codons chosen.

### Phase 2: Targeted Hairpin Repair

After assembly completes, the full sequence (`UTR + CDS + stop codon`) is checked for hairpin violations. If any 50 bp chunk has >1 hairpin, the repair phase attempts to fix them through targeted codon swaps.

**Single-codon swaps (primary):**

For each hairpin stem pair (two 3-mer positions `i` and `j` that are reverse complements 7-12 bp apart) in a failing chunk:
1. Identify which codon indices overlap either stem position.
2. For each overlapping codon, try every synonymous alternative that:
   - Breaks the specific RC match between the two stems.
   - Does not increase the total hairpin violation count across all overlapping 50 bp chunks (net-effect check).
   - Has not been tried before (oscillation prevention via `tried` set).
3. Among valid fixes, pick the one with the highest codon usage weight (to preserve CAI).
4. Apply the swap and re-check the full sequence before attempting the next fix.

Both the old and new codon are added to the `tried` set after each swap, preventing the algorithm from flip-flopping between two codons that each break one hairpin but create another.

**Two-codon pair swaps (fallback):**

When no single-codon swap can improve the situation — typically because two overlapping hairpins share stem bases (e.g., `GCA(loop)TGC` and `TGC(loop)GCA`) — the algorithm falls back to pair swaps:
1. For each failing 50 bp chunk, identify all codon indices that overlap the chunk (~17 codons).
2. Enumerate all pairs of codon-synonym substitutions within the chunk.
3. For each pair, compute the local hairpin violation count across all overlapping chunks.
4. Accept the pair that reduces violations the most, breaking ties by total codon weight.
5. Apply and re-check.

The repair runs for up to 20 passes. A summary is logged showing how many codons were touched and which were modified most frequently.

### Phase 3: Validation and Restart Logic

After repair, the algorithm appends a stop codon (`TAA`), selects an RBS, and constructs the full transcript (`UTR + CDS`). Four checks are run on the full transcript:

1. **PromoterChecker** — scans for internal promoter sequences using a 29 bp PWM with threshold 9.134, on both strands.
2. **hairpin_checker** — slides 50 bp chunks (step 25) and fails if any chunk has >1 hairpin (3 bp stems, 4-9 bp loops).
3. **ForbiddenSequenceChecker** — searches for restriction enzyme recognition sites on both strands.
4. **InternalRBSChecker** — currently disabled (always passes).

If all checks pass, the transcript is returned. Otherwise:

**Early stop detection:** The algorithm tracks the hairpin failure signature (the exact text of which hairpins failed). If the same signature repeats for 15 consecutive restarts, the region is deemed inherently unfixable — typically because single-codon amino acids like Methionine (`ATG`) force an unavoidable hairpin pattern — and the algorithm raises a `RuntimeError` immediately instead of exhausting all 500 restarts.

If no early stop triggers and all restarts are exhausted, a `RuntimeError` is raised. The benchmarker catches this exception and logs the failure, allowing the run to continue to the next protein.

### Known Limitations

The protein cdhC (Acetyl-CoA decarbonylase/synthase complex subunit beta, 524 AAs) contains the stretch `QPMAMQPMPMQMP` — four Methionines interspersed with Prolines. Methionine has only one codon (`ATG`), so the repeating `ATG...ATG...ATG` pattern is immutable. Exhaustive enumeration of all 4,096 synonymous codon combinations for this 12-amino-acid region confirms that every combination produces at least 2 hairpins within a 15 bp span. Both hairpins always fall within the same 50 bp checker chunks, so they cannot be distributed across chunk boundaries. This protein cannot pass the hairpin checker under the current parameters with any synonymous codon assignment.
