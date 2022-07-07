//
// TEAL measures telomere content, and estimates telomeric length based on WGS read data.
//

include { PICARD_COLLECTWGSMETRICS as COLLECTWGSMETRICS } from '../../modules/nf-core/modules/picard/collectwgsmetrics/main'

include { TEAL as TEAL_PROCESS } from '../../modules/scwatts/nextflow_modules/teal/main'

workflow TEAL {
  take:
    ch_inputs           // channel: [val(meta), tumor_bam, normal_bam, cobalt_dir, purple_dir]
    ref_data_genome_dir //    file: /path/to/genome_dir/
    ref_data_genome_fn  //     val: genome name

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Obtain WGS metrics
    // NOTE(SW): here I remove duplicate files so that we only process each input once
    // channel: [val(meta_teal), bam]
    ch_metrics_inputs = WorkflowTeal.get_metrics_inputs(ch_inputs)
    // channel: [val(meta_teal), bam]
    ch_metrics_inputs = WorkflowTeal.get_unique_input_files(ch_metrics_inputs)
    COLLECTWGSMETRICS(
      ch_metrics_inputs,
      "${ref_data_genome_dir}/${ref_data_genome_fn}",
    )
    ch_versions = ch_versions.mix(COLLECTWGSMETRICS.out.versions)

    // Run TEAL
    // NOTE(SW): performed here to avoid WorkflowHmftools import in WorkflowTeal
    // channel: [val(meta_teal), wgs_metrics, bam]
    ch_bams_and_metrics_unsorted = WorkflowHmftools.group_by_meta(
      ch_metrics_inputs,
      COLLECTWGSMETRICS.out.metrics,
    )
    // channel: [val(meta), tumor_bam, normal_bam, tumor_wgs_metrics, normal_wgs_metrics]
    ch_bam_and_metrics = WorkflowTeal.sort_bams_and_metrics(ch_bams_and_metrics_unsorted)
    // channel: [val(meta), cobalt_dir, purple_dir]
    ch_cobalt_purple = ch_inputs
      .map { meta, tbam, nbam, cobalt_dir, purple_dir ->
        [meta, cobalt_dir, purple_dir]
      }
    // channel: [val(meta), tumor_bam, normal_bam, tumor_wgs_metrics, normal_wgs_metrics, cobalt_dir, purple_dir]
    ch_teal_inputs = WorkflowHmftools.group_by_meta(
      ch_bam_and_metrics,
      ch_cobalt_purple,
    )
    TEAL_PROCESS(
      ch_teal_inputs,
    )
    ch_versions = ch_versions.mix(TEAL_PROCESS.out.versions)

  emit:
    versions = ch_versions // channel: [versions.yml]
}
