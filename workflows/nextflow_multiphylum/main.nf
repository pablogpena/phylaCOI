nextflow.enable.dsl = 2

/*
 * Multiphylum Nextflow workflow.
 * Per-phylum tasks are launched after VSEARCH taxonomy assignment.
 */

params.raw_fasta = params.raw_fasta ?: 'data/raw/eKOI_metabarcoding.fasta'
params.reference_db = params.reference_db ?: 'data/raw/eKOI_database.fasta'
params.sample_metadata = params.sample_metadata ?: 'data/raw/KOI_metadata.csv'
params.ocean_metadata = params.ocean_metadata ?: 'data/raw/ocean_metadata.csv'
params.output_root = params.output_root ?: 'data'

params.vsearch_identity = params.vsearch_identity ?: 0.84
params.otu_identity = params.otu_identity ?: 0.97
params.run_mafft = params.run_mafft ?: 1
params.min_otus = params.min_otus ?: 10
params.max_parallel_phyla = params.max_parallel_phyla ?: 2
params.max_otu_seqs = params.max_otu_seqs ?: 500
params.cluster_radius_km = params.cluster_radius_km ?: 5
params.write_maps = params.write_maps ?: true

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

    python ${repo_dir}/workflows/nextflow_multiphylum/src/flatten_vsearch_results.py \
      -i vsearch_nested \
      -o vsearch_results
    """
}

process ABUNDANCE_ONE {
    tag "${phylum}"
    maxForks params.max_parallel_phyla

    publishDir "${output_root}/abundance", mode: 'copy', overwrite: true

    input:
    tuple val(phylum), path(phylum_dir), path(seq_headers), path(sample_metadata)

    output:
    tuple val(phylum), path("abundance_out/${phylum}"), emit: abundance_dir

    script:
    """
    python ${repo_dir}/workflows/nextflow_multiphylum/src/run_abundance_one_phylum.py \
      --phylum-dir ${phylum_dir} \
      --output-dir abundance_out/${phylum} \
      --names ${seq_headers} \
      --metadata ${sample_metadata} \
      --mafft ${params.run_mafft} \
      --repo-dir ${repo_dir}
    """
}

process OTU_GENERATION_ONE {
    tag "${phylum}"
    maxForks params.max_parallel_phyla

    publishDir "${output_root}/otus", mode: 'copy', overwrite: true

    input:
    tuple val(phylum), path(abundance_dir)

    output:
    tuple val(phylum), path("otus_out/${phylum}"), emit: otus_dir

    script:
    """
    python ${repo_dir}/workflows/nextflow_multiphylum/src/run_otu_generation_one_phylum.py \
      --abundance-dir ${abundance_dir} \
      --output-dir otus_out/${phylum} \
      --identity ${params.otu_identity} \
      --repo-dir ${repo_dir}
    """
}

process INFORMATIVE_OTUS_ONE {
    tag "${phylum}"
    maxForks params.max_parallel_phyla

    publishDir "${output_root}/otus", mode: 'copy', overwrite: true

    input:
    tuple val(phylum), path(abundance_dir, stageAs: 'abundance_input'), path(otus_dir, stageAs: 'otus_input')

    output:
    tuple val(phylum), path("otus_out/${phylum}"), emit: otus_dir

    script:
    """
    mkdir -p abundance_root otus_root otus_out
    cp -Lr abundance_input abundance_root/${phylum}
    cp -Lr otus_input otus_root/${phylum}

    Rscript ${repo_dir}/workflows/nextflow_multiphylum/src/run_informative_otus_one_phylum.R \
      --phylum ${phylum} \
      --abundance-root abundance_root \
      --otus-root otus_root

    cp -Lr otus_root/${phylum} otus_out/${phylum}
    """
}

process OTU_METRICS_ONE {
    tag "${phylum}"
    maxForks params.max_parallel_phyla

    publishDir "${output_root}/otus", mode: 'copy', overwrite: true

    input:
    tuple val(phylum), path(abundance_dir, stageAs: 'abundance_input'), path(otus_dir, stageAs: 'otus_input')

    output:
    tuple val(phylum), path("otus_out/${phylum}"), emit: otus_dir

    script:
    def map_flag = params.write_maps ? '' : '--no-maps'
    """
    mkdir -p abundance_root otus_root otus_out
    cp -Lr abundance_input abundance_root/${phylum}
    cp -Lr otus_input otus_root/${phylum}

    Rscript ${repo_dir}/workflows/nextflow_multiphylum/src/run_otu_metrics_one_phylum.R \
      --phylum ${phylum} \
      --abundance-root abundance_root \
      --otus-root otus_root \
      --max-otu-seqs ${params.max_otu_seqs} \
      --cluster-radius-km ${params.cluster_radius_km} \
      ${map_flag}

    cp -Lr otus_root/${phylum} otus_out/${phylum}
    """
}

process COMBINE_OTU_METRICS {
    tag 'combine_otu_metrics'

    publishDir "${output_root}", mode: 'copy', pattern: 'otus', overwrite: true

    input:
    path phylum_otus_dirs

    output:
    path 'otus', emit: otus_root

    script:
    """
    mkdir -p otus
    cp -Lr ${phylum_otus_dirs} otus/

    Rscript ${repo_dir}/workflows/nextflow_multiphylum/src/combine_otu_metrics.R \
      --otus-root otus
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
    """
    Rscript ${repo_dir}/scripts/6_haplotype_clustering/run_all_haplotypes_clustering.R \
      -i ${otus} \
      -o all_haplotypes \
      --metadata ${ocean_metadata}
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

    VSEARCH_TAXONOMY.out.vsearch_results
        .combine(RAW_PREPROCESS.out.seq_headers)
        .combine(sample_metadata_ch)
        .flatMap { vsearch_root, seq_headers, sample_metadata ->
            def phylum_dirs = vsearch_root.toFile().listFiles() ?: []
            phylum_dirs
                .findAll { it.isDirectory() }
                .sort { it.name }
                .collect { phylum_dir -> tuple(phylum_dir.name, phylum_dir.toPath(), seq_headers, sample_metadata) }
        }
        .set { abundance_inputs_ch }

    ABUNDANCE_ONE(abundance_inputs_ch)
    OTU_GENERATION_ONE(ABUNDANCE_ONE.out.abundance_dir)

    ABUNDANCE_ONE.out.abundance_dir
        .join(OTU_GENERATION_ONE.out.otus_dir)
        .set { informative_inputs_ch }

    INFORMATIVE_OTUS_ONE(informative_inputs_ch)

    ABUNDANCE_ONE.out.abundance_dir
        .join(INFORMATIVE_OTUS_ONE.out.otus_dir)
        .set { otu_metrics_inputs_ch }

    OTU_METRICS_ONE(otu_metrics_inputs_ch)

    OTU_METRICS_ONE.out.otus_dir
        .map { phylum, otus_dir -> otus_dir }
        .collect()
        .set { otu_metric_dirs_ch }

    COMBINE_OTU_METRICS(otu_metric_dirs_ch)

    METRICS_ANALYSIS(COMBINE_OTU_METRICS.out.otus_root)
    HAPLOTYPE_CLUSTERING_ALL(COMBINE_OTU_METRICS.out.otus_root, ocean_metadata_ch)
    HAPLOTYPE_CLUSTERING_SAME_DIFF(COMBINE_OTU_METRICS.out.otus_root, HAPLOTYPE_CLUSTERING_ALL.out.all_haplotypes)
}
