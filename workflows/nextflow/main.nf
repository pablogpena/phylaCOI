nextflow.enable.dsl = 2

/*
 * Phase 1 Nextflow wrapper.
 * This workflow orchestrates the existing scripts without changing them.
 */

repo_dir = file(params.repo_dir ?: workflow.launchDir).toRealPath()
output_root = file(params.output_root)

process RAW_PREPROCESS {
    tag 'raw_fasta'

    publishDir "${output_root}/raw", mode: 'copy', pattern: 'seq_headers.txt'
    publishDir "${output_root}/procesed", mode: 'copy', pattern: 'eKOI_metabarcoding_cleaned.fasta'

    input:
    path raw_fasta

    output:
    path 'seq_headers.txt', emit: seq_headers
    path 'eKOI_metabarcoding_cleaned.fasta', emit: cleaned_fasta

    script:
    """
    python ${repo_dir}/scripts/1_raw_data_processing/fasta_preprocess.py \
      -i ${raw_fasta} \
      -t seq_headers.txt \
      -o eKOI_metabarcoding_cleaned.fasta
    """
}

process VSEARCH_TAXONOMY {
    tag 'vsearch_taxonomy'

    publishDir "${output_root}", mode: 'copy', pattern: 'vsearch_results', overwrite: true

    input:
    path cleaned_fasta
    path reference_db

    output:
    path 'vsearch_results', emit: vsearch_results

    script:
    """
    python ${repo_dir}/scripts/1_raw_data_processing/vsearch_taxonomy.py \
      -d ${reference_db} \
      -i ${params.vsearch_identity} \
      -f . \
      -o vsearch_nested

    python ${repo_dir}/workflows/nextflow/src/flatten_vsearch_results.py \
      -i vsearch_nested \
      -o vsearch_results
    """
}

process ABUNDANCE_ESTIMATION {
    tag 'abundance'

    publishDir "${output_root}", mode: 'copy', pattern: 'abundance', overwrite: true

    input:
    path vsearch_results
    path seq_headers
    path sample_metadata

    output:
    path 'abundance', emit: abundance

    script:
    """
    python ${repo_dir}/scripts/2_abundance_estimation/get_abundance.py \
      -i ${vsearch_results} \
      -o abundance \
      -n ${seq_headers} \
      -m ${sample_metadata} \
      -a ${params.run_mafft}
    """
}

process OTU_GENERATION {
    tag 'otus'

    publishDir "${output_root}", mode: 'copy', pattern: 'otus', overwrite: true

    input:
    path abundance

    output:
    path 'otus', emit: otus

    script:
    """
    python ${repo_dir}/scripts/3_otu_generation/generate_otus.py \
      -a ${abundance} \
      -o otus \
      -i ${params.otu_identity}
    """
}

process INFORMATIVE_OTUS {
    tag 'informative_otus'

    publishDir "${output_root}", mode: 'copy', pattern: 'otus', overwrite: true

    input:
    path abundance
    path otus_in, stageAs: 'otus_input'

    output:
    path 'otus', emit: otus

    script:
    """
    cp -Lr otus_input otus

    Rscript ${repo_dir}/scripts/3_otu_generation/get_informative_otus.R \
      -a ${abundance} \
      -o otus
    """
}

process OTU_METRICS {
    tag 'otu_metrics'

    publishDir "${output_root}", mode: 'copy', pattern: 'otus', overwrite: true

    input:
    path abundance
    path otus_in, stageAs: 'otus_input'

    output:
    path 'otus', emit: otus

    script:
    def map_flag = params.write_maps ? '' : '--no-maps'
    """
    cp -Lr otus_input otus

    Rscript ${repo_dir}/scripts/4_otu_metrics/get_div_abun_conn.R \
      -a ${abundance} \
      -o otus \
      ${map_flag}
    """
}

process METRICS_ANALYSIS {
    tag 'metrics_analysis'

    publishDir "${output_root}/analysis", mode: 'copy', pattern: 'otu_metrics_summary', overwrite: true

    input:
    path otus

    output:
    path 'otu_metrics_summary', emit: summary

    script:
    """
    Rscript ${repo_dir}/scripts/5_metrics_analysis/analyze_otu_metrics.R \
      -i ${otus}/div_abun_conn_master.csv \
      -o otu_metrics_summary \
      --min-otus ${params.min_otus}
    """
}

process HAPLOTYPE_CLUSTERING_ALL {
    tag 'haplotype_all'

    publishDir "${output_root}/analysis/haplotype_clustering", mode: 'copy', pattern: 'all_haplotypes', overwrite: true

    input:
    path otus
    path ocean_metadata

    output:
    path 'all_haplotypes', emit: all_haplotypes

    script:
    def map_flag = params.write_maps ? '' : '--no-maps'
    """
    Rscript ${repo_dir}/scripts/6_haplotype_clustering/run_all_haplotypes_clustering.R \
      -i ${otus} \
      -o all_haplotypes \
      --metadata ${ocean_metadata} \
      ${map_flag}
    """
}

process HAPLOTYPE_CLUSTERING_SAME_DIFF {
    tag 'haplotype_same_diff'

    publishDir "${output_root}/analysis/haplotype_clustering", mode: 'copy', pattern: 'same_vs_diff_currents', overwrite: true

    input:
    path otus
    path all_haplotypes

    output:
    path 'same_vs_diff_currents', emit: same_vs_diff_currents

    script:
    """
    Rscript ${repo_dir}/scripts/6_haplotype_clustering/run_same_vs_diff_currents.R \
      -i ${otus} \
      -o same_vs_diff_currents \
      --all-output ${all_haplotypes}
    """
}

workflow {
    raw_fasta_ch = Channel.fromPath(params.raw_fasta, checkIfExists: true)
    reference_db_ch = Channel.fromPath(params.reference_db, checkIfExists: true)
    sample_metadata_ch = Channel.fromPath(params.sample_metadata, checkIfExists: true)
    ocean_metadata_ch = Channel.fromPath(params.ocean_metadata, checkIfExists: true)

    RAW_PREPROCESS(raw_fasta_ch)
    VSEARCH_TAXONOMY(RAW_PREPROCESS.out.cleaned_fasta, reference_db_ch)
    ABUNDANCE_ESTIMATION(VSEARCH_TAXONOMY.out.vsearch_results, RAW_PREPROCESS.out.seq_headers, sample_metadata_ch)
    OTU_GENERATION(ABUNDANCE_ESTIMATION.out.abundance)
    INFORMATIVE_OTUS(ABUNDANCE_ESTIMATION.out.abundance, OTU_GENERATION.out.otus)
    OTU_METRICS(ABUNDANCE_ESTIMATION.out.abundance, INFORMATIVE_OTUS.out.otus)

    METRICS_ANALYSIS(OTU_METRICS.out.otus)
    HAPLOTYPE_CLUSTERING_ALL(OTU_METRICS.out.otus, ocean_metadata_ch)
    HAPLOTYPE_CLUSTERING_SAME_DIFF(OTU_METRICS.out.otus, HAPLOTYPE_CLUSTERING_ALL.out.all_haplotypes)
}
