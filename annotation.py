# Annotate vcfeval output / failure mode types

import argparse
from cyvcf2 import VCF
from pybedtools import BedTool
import pyBigWig
from pathlib import Path
import pandas as pd
 
 
# helper functions
 
# reads a raw VCF and returns a lookup dict: "chrom_pos_ref_alt" -> filter string
def build_raw_filter_lookup(vcf_path: Path) -> dict:
    lookup = {}
    for var in VCF(str(vcf_path)):
        alt = ",".join(var.ALT) if var.ALT else "."
        key = f"{var.CHROM}_{var.POS}_{var.REF}_{alt}"
        filt = var.FILTER if var.FILTER else "PASS"
        lookup[key] = filt
    return lookup
 
 
# returns SNV if all alleles are single base substitutions, MNV if same length, otherwise INDEL
def get_variant_type(ref: str, alt: str) -> str:
    alt_alleles = alt.split(",")
    if all(len(a) == 1 and len(ref) == 1 for a in alt_alleles):
        return "SNV"
    if all(len(a) == len(ref) for a in alt_alleles):
        return "MNV"
    return "INDEL"
 
 
# write a dict of {sheet_name: (headers, rows)} to separate TSV files in output_dir
def write_tsvs(discordance_dir: Path, sheets: dict, headers_map: dict):
    for sheet_name, rows in sheets.items():
        headers = headers_map[sheet_name]
        df = pd.DataFrame(rows, columns=headers)
        filename = sheet_name.lower().replace(" ", "_") + ".tsv"
        df.to_csv(discordance_dir / filename, sep="\t", index=False)
        print(f"[annotate] Saved {discordance_dir / filename}")
 

# main function

def annotate_vcfeval_output(
        vcfeval_dir: Path,
        baseline_vcf: Path,
        query_vcf: Path,
        gc_bw: Path,
        low_complexity_bed: Path,
        segmental_dup_bed: Path,
        discordance_dir: Path,  # directory to write TSVs into
):
    # convert each BED file into a BedTool object to intersect variants with them
    lowcomp = BedTool(str(low_complexity_bed))
    seg_dup = BedTool(str(segmental_dup_bed))
 
    # load raw filter lookups from original caller VCFs
    baseline_filter_lookup = build_raw_filter_lookup(baseline_vcf)
    query_filter_lookup    = build_raw_filter_lookup(query_vcf)
 
    # output.vcf.gz is the combined vcfeval output produced by module 2
    vcf_path = vcfeval_dir / "output.vcf.gz"
 
    # pass 1: collect all variants and build a BED of their positions
    variants  = []
    bed_lines = []
 
    for var in VCF(str(vcf_path)):
        start = var.POS - 1
        end   = start + len(var.REF)
        variants.append(var)
        bed_lines.append(f"{var.CHROM}\t{start}\t{end}")
 
    if not variants:
        print("[annotate] No variants found in vcfeval output — skipping annotation.")
        return
 
    # build one BED for all variants to intersect with the difficult region BED files
    all_variants_bed = BedTool("\n".join(bed_lines), from_string=True)
    bw = pyBigWig.open(str(gc_bw))
 
    # intersect once per annotation track, return hit keys as a set
    def get_hits(annotation_bed) -> set:
        df = all_variants_bed.intersect(annotation_bed, u=True).to_dataframe()
        if df.empty:
            return set()
        return set(df[["chrom", "start", "end"]].itertuples(index=False, name=None))
 
    lowcomp_hits = get_hits(lowcomp)
    seg_dup_hits = get_hits(seg_dup)
 
    # pass 2: build rows for each sheet
 
    sheets = {
        "All Variants":    [],
        "GT Discordant":   [],
        "High GC":         [],
        "High AT":         [],
        "Low Complexity":  [],
        "Segmental Dup":   [],
        "Caller Specific": [],
    }
 
    for var in variants:
        chrom  = var.CHROM
        pos    = var.POS
        ref    = var.REF
        alt    = ",".join(var.ALT) if var.ALT else "."
        start  = pos - 1
        end    = start + len(ref)
        key    = (chrom, start, end)
        var_id = f"{chrom}_{pos}_{ref}_{alt}"
 
        # variant type — SNV, MNV, or INDEL
        variant_type = get_variant_type(ref, alt)
 
        # vcfeval status for baseline and call samples
        base_status = var.INFO.get("BASE", ".")
        call_status = var.INFO.get("CALL", ".")

        # skip TPs — concordant variants don't belong in any failure mode output
        if base_status == "TP" and call_status == "TP":
            continue

        # skip variants ignored by both callers — neither contributed meaningful data
        if base_status == "IGN" and call_status == "IGN":
            continue
 
        # pull genotype strings from gt_bases
        gt_bases        = var.gt_bases
        baseline_gt_str = gt_bases[0] if gt_bases is not None and len(gt_bases) > 0 else "."
        calls_gt_str    = gt_bases[1] if gt_bases is not None and len(gt_bases) > 1 else "."
 
        # GT discordant
        gt_discordant = (
            "YES" if (
                "FP_CA" in (base_status, call_status) and
                "." not in baseline_gt_str and
                "." not in calls_gt_str
            )
            else "NO"
        )
 
        # caller specific
        if base_status in ("IGN", ".") and call_status not in ("IGN", "."):
            caller_specific = "YES"
            missing_in      = "BASELINE"
            raw_filter      = baseline_filter_lookup.get(var_id, "NOT_CALLED")
        elif call_status in ("IGN", ".") and base_status not in ("IGN", "."):
            caller_specific = "YES"
            missing_in      = "QUERY"
            raw_filter      = query_filter_lookup.get(var_id, "NOT_CALLED")
        else:
            caller_specific = "NO"
            missing_in      = None
            raw_filter      = None
 
        # region flags

        gc_val = bw.stats(chrom, pos - 1, pos, type="mean")[0]
        if gc_val is None:
            high_gc_flag = "NO"
            high_at_flag = "NO"
        elif gc_val >= 80.0:
            high_gc_flag = "YES"
            high_at_flag = "NO"
        elif gc_val <= 20.0:
            high_gc_flag = "NO"
            high_at_flag = "YES"
        else:
            high_gc_flag = "NO"
            high_at_flag = "NO"
        lowcomp_flag = "YES" if key in lowcomp_hits else "NO"
        seg_dup_flag = "YES" if key in seg_dup_hits else "NO"
 
        # All Variants
        sheets["All Variants"].append([
            var_id, variant_type, base_status, call_status,
            gt_discordant, high_gc_flag, high_at_flag,
            lowcomp_flag, seg_dup_flag, caller_specific
        ])
 
        # per-failure-mode sheets
        if gt_discordant == "YES":
            sheets["GT Discordant"].append([
                var_id, chrom, pos, ref, alt, variant_type, base_status, call_status,
                baseline_gt_str, calls_gt_str
            ])
 
        if high_gc_flag == "YES":
            sheets["High GC"].append([
                var_id, chrom, pos, ref, alt, variant_type, base_status, call_status,
                gc_val if gc_val is not None else "N/A"
            ])
 
        if high_at_flag == "YES":
            sheets["High AT"].append([
                var_id, chrom, pos, ref, alt, variant_type, base_status, call_status,
                100 - gc_val if gc_val is not None else "N/A"
            ])
 
        if lowcomp_flag == "YES":
            sheets["Low Complexity"].append([
                var_id, chrom, pos, ref, alt, variant_type, base_status, call_status
            ])
 
        if seg_dup_flag == "YES":
            sheets["Segmental Dup"].append([
                var_id, chrom, pos, ref, alt, variant_type, base_status, call_status
            ])
 
        if caller_specific == "YES":
            sheets["Caller Specific"].append([
                var_id, chrom, pos, ref, alt, variant_type, base_status, call_status,
                missing_in, raw_filter
            ])
    bw.close()
    headers_map = {
        "All Variants":    ["VAR_ID", "VARIANT_TYPE", "BASE_STATUS", "CALL_STATUS",
                            "GT_DISCORDANT", "HIGH_GC", "HIGH_AT", "LOW_COMPLEXITY",
                            "SEGMENTAL_DUPLICATION", "CALLER_SPECIFIC"],
        "GT Discordant":   ["VAR_ID", "CHROM", "POS", "REF", "ALT", "VARIANT_TYPE", "BASE_STATUS", "CALL_STATUS", "BASELINE_GT", "CALLS_GT"],
        "High GC":         ["VAR_ID", "CHROM", "POS", "REF", "ALT", "VARIANT_TYPE", "BASE_STATUS", "CALL_STATUS", "GC_PERCENT"],
        "High AT":         ["VAR_ID", "CHROM", "POS", "REF", "ALT", "VARIANT_TYPE", "BASE_STATUS", "CALL_STATUS", "AT_PERCENT"],
        "Low Complexity":  ["VAR_ID", "CHROM", "POS", "REF", "ALT", "VARIANT_TYPE", "BASE_STATUS", "CALL_STATUS"],
        "Segmental Dup":   ["VAR_ID", "CHROM", "POS", "REF", "ALT", "VARIANT_TYPE", "BASE_STATUS", "CALL_STATUS"],
        "Caller Specific": ["VAR_ID", "CHROM", "POS", "REF", "ALT", "VARIANT_TYPE", "BASE_STATUS", "CALL_STATUS", "MISSING_IN", "RAW_FILTER"],
    }
 
    write_tsvs(discordance_dir, sheets, headers_map)
 
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--vcfeval_dir",        required=True)
    parser.add_argument("--baseline",           required=True)
    parser.add_argument("--query",              required=True)
    parser.add_argument("--gc_bw",              required=True)
    parser.add_argument("--low_complexity_bed", required=True)
    parser.add_argument("--segmental_dup_bed",  required=True)
    parser.add_argument("--discordance_dir",    required=True)

    args = parser.parse_args()

    annotate_vcfeval_output(
        vcfeval_dir        = Path(args.vcfeval_dir),
        baseline_vcf       = Path(args.baseline),
        query_vcf          = Path(args.query),
        gc_bw              = Path(args.gc_bw),
        low_complexity_bed = Path(args.low_complexity_bed),
        segmental_dup_bed  = Path(args.segmental_dup_bed),
        discordance_dir    = Path(args.discordance_dir),
    )