# Run vcfeval 

import argparse
import subprocess
from pathlib import Path
import requests

# helper functions 

# check if the vcfs have .tbi indexes, if not index the input .vcf.gz files
def index_vcf_if_needed(vcf: Path) -> None: 
    index = Path(str(vcf) + ".tbi")

    if not index.exists():
        subprocess.run(["bcftools", "index", "-t", str(vcf)], check=True)


# convert the reference build .fa.gz to .sdf for rtg to use
def prepare_reference(ucsc_url: str, fasta: Path, sdf: Path) -> None:

    if not fasta.exists():
        print(f"[reference] Downloading FASTA from {ucsc_url}")
        with requests.get(ucsc_url, stream=True) as r:
            r.raise_for_status() 
            with open(fasta, "wb") as f: 
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

    if not sdf.exists():
        print(f"[reference] Creating SDF at {sdf}")
        subprocess.run(["rtg", "format", "-o", str(sdf), str(fasta)], check=True)


# main function

# run vcfeval germline OR somatic
def run_vcfeval(args, vcfeval_dir: Path) -> None:

    # ensure VCFs are indexed
    index_vcf_if_needed(args.baseline)
    index_vcf_if_needed(args.query)

    # ensure reference FASTA + SDF exist
    prepare_reference(args.ucsc_url, args.ref_fa, args.ref_sdf)

    # build vcfeval command
    cmd = [
        "rtg", "vcfeval",
        "-b", str(args.baseline),
        "-c", str(args.query),
        "-t", str(args.ref_sdf),
        "-o", str(vcfeval_dir),
        "--output-mode=combine",
    ]

    # restrict to target regions ONLY for targeted assays
    if args.assay_type in ["wes", "targeted_panel", "amplicon"] and args.targets_bed is not None:
        cmd.extend(["--bed-regions", str(args.targets_bed)])

    # squash ploidy for somatic mode
    if args.mode == "somatic":
        cmd.append("--squash-ploidy")

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True) 

    # secondary verification that vcfeval ran
    expected = vcfeval_dir / "summary.txt"
    if not expected.exists():
        raise RuntimeError(f"vcfeval may have failed — missing {expected}")
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline",    required=True)
    parser.add_argument("--query",       required=True)
    parser.add_argument("--vcfeval_dir", required=True)
    parser.add_argument("--mode",        required=True, choices=["germline", "somatic"])
    parser.add_argument("--assay_type",  required=True, choices=["wgs", "wes", "targeted_panel", "amplicon"])
    parser.add_argument("--targets_bed", required=False, default=None)
    parser.add_argument("--build",       required=True, choices=["hg19", "hg38"])

    args = parser.parse_args()

    ref_build = {
        "hg19": {
             "ucsc_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
             "fa":  Path("hg19.fa.gz"),
             "sdf": Path("hg19.sdf"),
             },
        "hg38": {
            "ucsc_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
            "fa":  Path("hg38.fa.gz"),
            "sdf": Path("hg38.sdf"),
            },
    }


    args.ucsc_url = ref_build[args.build]["ucsc_url"]
    args.ref_fa   = ref_build[args.build]["fa"]
    args.ref_sdf  = ref_build[args.build]["sdf"]
    args.baseline = Path(args.baseline)
    args.query    = Path(args.query)

    run_vcfeval(args, Path(args.vcfeval_dir))