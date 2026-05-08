# Prepare inputs for report

import argparse
from pathlib import Path
import pandas as pd
 
 
# main function
 
def prepare_report_inputs(
        summary_txt: Path,
        report_inputs_dir: Path,
        sample_name: str,
        assay_type: str,
        build: str,
        mode: str,
        baseline_name: str,   
        query_name: str,      
):
    # parse summary.txt
    # skip the dashes line (row 1) and read whitespace-separated values
    df = pd.read_csv(summary_txt, sep=r"\s+", skiprows=[1], engine="python")
 
    # get the None threshold row which represents unthresholded performance
    row = df[df["Threshold"].isna()].iloc[0]
 
    # overall counts with percentages
    tp    = row["True-pos-call"]
    fp    = row["False-pos"]
    fn    = row["False-neg"]
    total = tp + fp + fn
 
    counts = pd.DataFrame([{
        "Sample":          sample_name,
        "Total_Evaluated": int(total),
        "TP":              int(tp),
        "TP_pct":          f"{round(tp / total * 100, 2)}%",
        "FP":              int(fp),
        "FP_pct":          f"{round(fp / total * 100, 2)}%",
        "FN":              int(fn),
        "FN_pct":          f"{round(fn / total * 100, 2)}%",
    }])
 
    # overall performance
    performance = pd.DataFrame([{
        "Sample":    sample_name,
        "Precision": row["Precision"],
        "Recall":    row["Sensitivity"],
        "F1_score":  row["F-measure"],
    }])
 
    # sample metadata
    metadata = pd.DataFrame([{
        "Sample":     sample_name,
        "Assay_Type": assay_type,
        "Build":      build,
        "Mode":       mode,
        "Baseline":    baseline_name,
        "Query":       query_name,
    }])

    report_inputs_dir = Path(report_inputs_dir)
    
    counts.to_csv(report_inputs_dir / "overall_counts.tsv", sep="\t", index=False)
    performance.to_csv(report_inputs_dir / "overall_performance.tsv", sep="\t", index=False)
    metadata.to_csv(report_inputs_dir / "metadata.tsv", sep="\t", index=False)
    
    print(f"[report inputs] Saved to {report_inputs_dir}")
 

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--summary_txt",       required=True)
    parser.add_argument("--report_inputs_dir", required=True)
    parser.add_argument("--sample_name",       required=True)
    parser.add_argument("--assay_type",        required=True)
    parser.add_argument("--build",             required=True)
    parser.add_argument("--mode",              required=True)
    parser.add_argument("--baseline_name",     required=True)
    parser.add_argument("--query_name",        required=True)

    args = parser.parse_args()

    prepare_report_inputs(
        summary_txt       = Path(args.summary_txt),
        report_inputs_dir = Path(args.report_inputs_dir),
        sample_name       = args.sample_name,
        assay_type        = args.assay_type,
        build             = args.build,
        mode              = args.mode,
        baseline_name     = args.baseline_name,
        query_name        = args.query_name,
    )