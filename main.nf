#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.outdir      = "out"
params.samplesheet = null
params.mode        = null
params.build       = null
params.assay_type  = null
params.targets_bed = null

println "All input parameters: ${params}"


process RUN_VCFEVAL {
    container 'concordance-pipeline:1.0'
    publishDir "${params.outdir}/${sample}_${run_id}_${params.mode}/vcfeval", mode: 'copy'

    input:
    tuple path(baseline), path(query), val(sample), val(run_id)

    output:
    tuple path('vcfeval_out/'), val(sample), val(run_id), val(baseline.name), val(query.name)

    script:
    def targets = params.targets_bed ? "--targets_bed ${params.targets_bed}" : ""
    """
    python3 /pipeline/run_vcfeval.py \
        --baseline ${baseline} \
        --query ${query} \
        --build ${params.build} \
        --mode ${params.mode} \
        --assay_type ${params.assay_type} \
        --vcfeval_dir vcfeval_out \
        ${targets}
    """
}

process ANNOTATE {
    container 'concordance-pipeline:1.0'
    publishDir "${params.outdir}/${sample}_${run_id}_${params.mode}/discordance_modes", mode: 'copy'

    input:
    tuple path(vcfeval_dir), val(sample), val(run_id), val(baseline_name), val(query_name)
    tuple path(baseline), path(query), val(sample_), val(run_id_)

    output:
    tuple path('*.tsv'), val(sample), val(run_id)

    script:
    """
    python3 /pipeline/annotation.py \
        --vcfeval_dir ${vcfeval_dir} \
        --baseline ${baseline} \
        --query ${query} \
        --gc_bw /resources/${params.build}/${params.build}.gc5Base.bw \
        --low_complexity_bed /resources/${params.build}/${params.build}_low_complexity.bed.gz \
        --segmental_dup_bed /resources/${params.build}/${params.build}_segdups.bed.gz \
        --discordance_dir .
    """
}

process PREPARE_REPORT_INPUTS {
    container 'concordance-pipeline:1.0'
    publishDir "${params.outdir}/${sample}_${run_id}_${params.mode}/report_inputs", mode: 'copy'

    input:
    tuple path(vcfeval_dir), val(sample), val(run_id), val(baseline_name), val(query_name)

    output:
    tuple path('*.tsv'), val(sample), val(run_id)

    script:
    """
    python3 /pipeline/prep_report.py \
        --summary_txt ${vcfeval_dir}/summary.txt \
        --report_inputs_dir . \
        --sample_name ${sample} \
        --assay_type ${params.assay_type} \
        --build ${params.build} \
        --mode ${params.mode} \
        --baseline_name ${baseline_name} \
        --query_name ${query_name}
    """
}

process RENDER_REPORT {
    container 'concordance-pipeline:1.0'
    publishDir "${params.outdir}/${sample}_${run_id}_${params.mode}", mode: 'copy'

    input:
    tuple path(discordance_files), val(sample), val(run_id)
    tuple path(report_inputs_files), val(sample_), val(run_id_)

    output:
    path 'Report.html'
    path 'Report_files/'

    script:
    """
    cp /pipeline/Report.qmd .
    mkdir -p report_inputs discordance_modes
    mv metadata.tsv overall_counts.tsv overall_performance.tsv report_inputs/
    mv *.tsv discordance_modes/
    quarto render Report.qmd --to html
    """
}

workflow {

    if ( !params.samplesheet ) { error "Missing --samplesheet" }
    if ( !params.mode )        { error "Missing --mode" }
    if ( !params.build )       { error "Missing --build" }
    if ( !params.assay_type )  { error "Missing --assay_type" }

    samples_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(
            file(row.baseline),
            file(row.query),
            row.sample_name,
            row.run_id
        )}

    vcfeval_ch       = RUN_VCFEVAL(samples_ch)
    discordance_ch   = ANNOTATE(vcfeval_ch, samples_ch)
    report_inputs_ch = PREPARE_REPORT_INPUTS(vcfeval_ch)
    RENDER_REPORT(discordance_ch, report_inputs_ch)
}