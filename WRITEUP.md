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
