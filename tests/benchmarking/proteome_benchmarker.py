import csv
import json
import os
import time
import traceback
from multiprocessing import Pool, cpu_count
from statistics import mean

from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.models.rbs_option import RBSOption
from genedesign.models.transcript import Transcript
from genedesign.seq_utils.translate import Translate
from genedesign.transcript_designer import CACHE_VERSION, TranscriptDesigner

# ---------------------------------------------------------------------------
# Checkpoint / resume design rationale
# ---------------------------------------------------------------------------
# A proteome benchmark can take hours.  To avoid losing progress on a crash
# or keyboard-interrupt, results are saved incrementally to a JSON checkpoint
# file after each gene completes.
#
# Cache invalidation key: CACHE_VERSION constant in genedesign/transcript_designer.py
#
# TranscriptDesigner is the only file whose changes affect output correctness.
# Rather than hashing the whole file (which would bust the cache on every
# tuning edit such as changing max iterations or retry counts), we rely on a
# manually maintained CACHE_VERSION string defined there.  Developers must bump
# CACHE_VERSION only when the algorithm changes in a way that would produce
# different outputs for the same input — leave it unchanged for parameter-only
# edits whose results remain valid.
#
# If the version stored in the checkpoint differs from the current value, the
# checkpoint is discarded and the run starts from scratch.
#
# Checkpointing is disabled when a gene_filter is active (single-gene runs
# are fast and should not pollute the full-proteome checkpoint).
# ---------------------------------------------------------------------------

CHECKPOINT_FILE = "benchmark_checkpoint.json"


def _get_cache_key() -> str:
    """Return the current cache invalidation key for benchmark results.

    This is the CACHE_VERSION constant from transcript_designer.py.  Bump that
    constant (not this function) when algorithm changes would make cached
    results invalid.  Tuning-only edits (max iterations, retry counts, etc.)
    do not require a bump.
    """
    return CACHE_VERSION


def _serialize_success(result: dict) -> dict:
    """Convert a success result (containing a Transcript object) to JSON-safe dict."""
    t = result["transcript"]
    return {
        "gene": result["gene"],
        "protein": result["protein"],
        "codons": list(t.codons),
        "rbs_utr": t.rbs.utr,
        "rbs_cds": t.rbs.cds,
        "rbs_gene_name": t.rbs.gene_name,
        "rbs_first_six_aas": t.rbs.first_six_aas,
    }


def _deserialize_success(data: dict) -> dict:
    """Reconstruct a success result dict (with Transcript) from a checkpoint entry."""
    rbs = RBSOption(
        utr=data["rbs_utr"],
        cds=data["rbs_cds"],
        gene_name=data["rbs_gene_name"],
        first_six_aas=data["rbs_first_six_aas"],
    )
    transcript = Transcript(rbs=rbs, peptide=data["protein"], codons=data["codons"])
    return {"gene": data["gene"], "protein": data["protein"], "transcript": transcript}


def _load_checkpoint(path: str, current_hash: str) -> tuple[list, list, set]:
    """
    Load a checkpoint file if it exists and its stored hash matches current_hash.

    Returns (successful_results, error_results, completed_genes).
    If the file is absent or the hash mismatches, returns empty collections and
    prints a message explaining why the run starts fresh.
    """
    if not os.path.exists(path):
        return [], [], set()

    try:
        with open(path, "r") as f:
            data = json.load(f)
    except (json.JSONDecodeError, KeyError):
        print(f"Checkpoint file '{path}' is corrupt — starting fresh.")
        return [], [], set()

    stored_hash = data.get("transcript_designer_hash", "")
    if stored_hash != current_hash:
        print(
            f"CACHE_VERSION has changed ('{stored_hash}' → '{current_hash}') — "
            f"discarding checkpoint and starting fresh."
        )
        return [], [], set()

    successful_results = [_deserialize_success(r) for r in data.get("successes", [])]
    error_results = data.get("errors", [])
    completed_genes = {r["gene"] for r in successful_results} | {
        r["gene"] for r in error_results
    }
    print(
        f"Resuming from checkpoint: {len(completed_genes)} genes already done "
        f"({len(successful_results)} success, {len(error_results)} error)."
    )
    return successful_results, error_results, completed_genes


def _save_checkpoint(
    path: str, designer_hash: str, successes: list, errors: list
) -> None:
    """Atomically write current results to the checkpoint file."""
    tmp_path = path + ".tmp"
    payload = {
        "transcript_designer_hash": designer_hash,
        "successes": [_serialize_success(r) for r in successes],
        "errors": errors,
    }
    with open(tmp_path, "w") as f:
        json.dump(payload, f)
    os.replace(tmp_path, path)


def parse_fasta(fasta_file):
    """
    Parses the FASTA file to extract gene names and protein sequences.
    """
    sequences = {}
    current_gene = None
    current_sequence = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_gene:
                    sequences[current_gene] = "".join(current_sequence)
                gene_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith("GN="):
                        gene_name = part.split("=")[1]
                        break
                if not gene_name:
                    gene_name = line.split("|")[2].split(" ")[0]
                current_gene = gene_name
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_gene:
            sequences[current_gene] = "".join(current_sequence)

    return sequences


def _worker_init():
    """Initialize a TranscriptDesigner in each MP worker (logging suppressed)."""
    import logging
    logging.disable(logging.CRITICAL)
    global _designer
    _designer = TranscriptDesigner()
    _designer.initiate()


def _init_sequential():
    """Initialize a TranscriptDesigner in the main process (logging enabled)."""
    global _designer
    _designer = TranscriptDesigner()
    _designer.initiate()


def _process_gene(args):
    """Process a single gene."""
    gene, protein = args
    ignores = set()
    try:
        transcript = _designer.run(protein, ignores)
        print(f"  SUCCESS: {gene} (len={len(protein)})")
        return ("success", {"gene": gene, "protein": protein, "transcript": transcript})
    except Exception as e:
        print(f"  FAILED: {gene} (len={len(protein)}) — {e}")
        return ("error", {"gene": gene, "error": traceback.format_exc()})


def benchmark_proteome(fasta_file, gene_filter=None, use_mp=False):
    """
    Benchmarks the proteome using TranscriptDesigner.

    When use_mp=True, genes are processed in parallel across all available
    CPUs and logging is suppressed in workers.  When False, genes are
    processed sequentially in the main process with logging enabled.

    When gene_filter is None, results are checkpointed after every completed
    gene so that an interrupted run can resume without reprocessing finished
    genes.  See module-level docstring for the cache-invalidation strategy.
    """
    proteome = parse_fasta(fasta_file)
    if gene_filter:
        proteome = {g: p for g, p in proteome.items() if g == gene_filter}
        if not proteome:
            print(f"Gene '{gene_filter}' not found in FASTA file.")
            return [], []

    use_checkpoint = gene_filter is None

    if use_checkpoint:
        cache_key = _get_cache_key()
        successful_results, error_results, completed_genes = _load_checkpoint(
            CHECKPOINT_FILE, cache_key
        )
        remaining = {g: p for g, p in proteome.items() if g not in completed_genes}
    else:
        cache_key = None
        successful_results, error_results = [], []
        remaining = proteome

    if not remaining:
        print("All genes already completed — nothing left to process.")
        return successful_results, error_results

    if use_mp:
        num_workers = min(cpu_count(), len(remaining))
        print(f"Processing {len(remaining)} genes with {num_workers} workers (MP, logging disabled)...")
        with Pool(num_workers, initializer=_worker_init) as pool:
            for status, result in pool.imap_unordered(_process_gene, remaining.items()):
                if status == "success":
                    successful_results.append(result)
                else:
                    error_results.append(result)
                if use_checkpoint:
                    _save_checkpoint(
                        CHECKPOINT_FILE, cache_key, successful_results, error_results
                    )
    else:
        print(f"Processing {len(remaining)} genes (sequential, logging enabled)...")
        _init_sequential()
        for args in remaining.items():
            status, result = _process_gene(args)
            if status == "success":
                successful_results.append(result)
            else:
                error_results.append(result)
            if use_checkpoint:
                _save_checkpoint(
                    CHECKPOINT_FILE, cache_key, successful_results, error_results
                )

    return successful_results, error_results


def analyze_errors(error_results):
    """
    Write the error analysis to a text file.
    """
    error_summary = {}
    with open("error_summary.txt", "w") as f:
        for error in error_results:
            error_message = error["error"].split("\n")[0]
            error_summary[error_message] = error_summary.get(error_message, 0) + 1
            f.write(f"Gene: {error['gene']}\n{error['error']}\n\n")

    return error_summary


def validate_transcripts(successful_results):
    """
    Validate the successful transcripts using various checkers, now including CodonChecker.
    """
    forbidden_checker = ForbiddenSequenceChecker()
    forbidden_checker.initiate()
    promoter_checker = PromoterChecker()
    promoter_checker.initiate()
    translator = Translate()
    translator.initiate()
    codon_checker = CodonChecker()  # Initialize CodonChecker
    codon_checker.initiate()  # Load the codon usage data

    validation_failures = []
    for result in successful_results:
        cds = "".join(result["transcript"].codons)
        try:
            # Check if CDS length is a multiple of 3
            if len(cds) % 3 != 0:
                raise ValueError("CDS length is not a multiple of 3.")

            # Verify that the translated protein matches the original protein
            original_protein = result["protein"]
            translated_protein = translator.run(cds)
            if original_protein != translated_protein:
                raise ValueError(
                    f"Translation mismatch: Original {original_protein}, Translated {translated_protein}"
                )

            # Ensure CDS starts with valid start codon and ends with stop codon
            if not (
                cds.startswith(("ATG", "GTG", "TTG"))
                and cds.endswith(("TAA", "TGA", "TAG"))
            ):
                raise ValueError(
                    "CDS does not start with a valid start codon or end with a valid stop codon."
                )
        except ValueError as e:
            validation_failures.append(
                {
                    "gene": result["gene"],
                    "protein": result["protein"],
                    "cds": cds,
                    "site": f"Translation or completeness error: {str(e)}",
                }
            )
            continue

        # Validate against hairpins, forbidden sequences, and internal promoters
        transcript_dna = result["transcript"].rbs.utr.upper() + cds
        passed_hairpin, hairpin_string = hairpin_checker(transcript_dna)
        if not passed_hairpin:
            formatted_hairpin = hairpin_string.replace("\n", " ").replace('"', "'")
            validation_failures.append(
                {
                    "gene": result["gene"],
                    "protein": result["protein"],
                    "cds": transcript_dna,
                    "site": f"Hairpin detected: {formatted_hairpin}",
                }
            )

        passed_forbidden, forbidden_site = forbidden_checker.run(transcript_dna)
        if not passed_forbidden:
            validation_failures.append(
                {
                    "gene": result["gene"],
                    "protein": result["protein"],
                    "cds": transcript_dna,
                    "site": f"Forbidden sequence: {forbidden_site}",
                }
            )

        passed_promoter, found_promoter = promoter_checker.run(transcript_dna)
        if not passed_promoter:
            validation_failures.append(
                {
                    "gene": result["gene"],
                    "protein": result["protein"],
                    "cds": transcript_dna,
                    "site": f"Constitutive promoter detected: {found_promoter}"
                    if found_promoter
                    else "Constitutive promoter detected",
                }
            )

        codons_above_board, codon_diversity, rare_codon_count, cai_value = (
            codon_checker.run(result["transcript"].codons)
        )
        if not codons_above_board:
            validation_failures.append(
                {
                    "gene": result["gene"],
                    "protein": result["protein"],
                    "cds": cds,
                    "site": f"Codon usage check failed: Diversity={codon_diversity}, Rare Codons={rare_codon_count}, CAI={cai_value}",
                }
            )

    return validation_failures


def write_validation_report(validation_failures):
    """
    Writes validation results to a TSV file.
    """
    with open("validation_failures.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gene", "protein", "cds", "site"])
        for failure in validation_failures:
            writer.writerow(
                [failure["gene"], failure["protein"], failure["cds"], failure["site"]]
            )


def generate_summary(
    total_genes, total_time, errors_summary, validation_failures
):
    """
    Generates a streamlined summary report categorizing validation failures by checker.
    """
    total_validation_failures = len(validation_failures)

    # Categorize failures by checker type
    checker_failures = {
        "Forbidden Sequence Checker": 0,
        "Hairpin Checker": 0,
        "Codon Usage Checker": 0,
        "Promoter Checker": 0,
        "Translation/Completeness Checker": 0,
    }

    # Increment the appropriate checker category based on the failure site
    for failure in validation_failures:
        site = failure["site"]
        if "Forbidden sequence" in site:
            checker_failures["Forbidden Sequence Checker"] += 1
        elif "Hairpin detected" in site:
            checker_failures["Hairpin Checker"] += 1
        elif "Codon usage check failed" in site:
            checker_failures["Codon Usage Checker"] += 1
        elif "Constitutive promoter detected" in site:
            checker_failures["Promoter Checker"] += 1
        elif "Translation or completeness error" in site:
            checker_failures["Translation/Completeness Checker"] += 1

    # Generate the summary report
    with open("summary_report.txt", "w") as f:
        f.write(f"Total genes processed: {total_genes}\n")
        f.write(f"Total wall-clock runtime: {total_time:.2f} seconds\n")
        f.write(f"Total exceptions: {sum(errors_summary.values())}\n")

        if errors_summary:
            f.write(f"\nTop 3 most common exceptions:\n")
            for error, count in sorted(
                errors_summary.items(), key=lambda x: x[1], reverse=True
            )[:3]:
                f.write(f"- {error}: {count} occurrences\n")
        else:
            f.write("No exceptions encountered.\n")

        f.write(f"\nTotal validation failures: {total_validation_failures}\n")

        # Categorize validation failures by checker
        f.write("\nValidation Failures by Checker:\n")
        for checker, count in checker_failures.items():
            f.write(f"- {checker}: {count} occurrences\n")


def _run_benchmark(fasta_file, gene_filter=None, use_mp=False):
    """
    Runs the complete benchmark process: parsing, running TranscriptDesigner, validating, and generating reports.
    """
    start_time = time.time()

    # Benchmark the proteome
    successful_results, error_results = benchmark_proteome(fasta_file, gene_filter, use_mp)

    # Analyze and log errors
    errors_summary = analyze_errors(error_results)

    # Validate the successful transcripts
    validation_failures = validate_transcripts(successful_results)

    total_time = time.time() - start_time

    # Write validation and error reports
    write_validation_report(validation_failures)

    # Generate the summary report
    total_genes = len(successful_results) + len(error_results)
    generate_summary(
        total_genes, total_time, errors_summary, validation_failures
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Proteome benchmark for TranscriptDesigner")
    parser.add_argument("gene_filter", nargs="?", default=None, help="Run only this gene (by name)")
    parser.add_argument("--mp", action="store_true", help="Enable multiprocessing (logging auto-disabled)")
    args = parser.parse_args()

    fasta_file = "tests/benchmarking/uniprotkb_proteome_UP000054015_2024_09_24.fasta"
    _run_benchmark(fasta_file, args.gene_filter, use_mp=args.mp)
