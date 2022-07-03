//
// TEAL measures telomere content, and estimates telomeric length based on WGS read data.
//

include { PICARD_COLLECTWGSMETRICS as COLLECTWGSMETRICS } from '../modules/nf-core/modules/picard/collectwgsmetrics/main'

include { TEAL as TEAL_PROCESS } from '../modules/scwatts/nextflow_modules/teal/main'

workflow TEAL {
  take:
    ch_inputs           // channel: [meta, tumor_bam, normal_bam, cobalt_dir, purple_dir]
    ref_data_genome_dir
    ref_data_genome_fn

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Obtain WGS metrics
    ch_metrics_inputs = WorkflowTeal.get_metrics_inputs(ch_inputs)
    ch_metrics_inputs = WorkflowTeal.get_unique_input_files(ch_metrics_inputs)
    COLLECTWGSMETRICS(
      ch_metrics_inputs,
      "${ref_data_genome_dir}/${ref_data_genome_fn}",
    )
    ch_versions = ch_versions.mix(COLLECTWGSMETRICS.out.versions)

    // Run TEAL
    // NOTE(SW): performed here to avoid WorkflowHmftools import in WorkflowTeal
    ch_bams_and_metrics_unsorted = WorkflowHmftools.group_by_meta(
      ch_metrics_inputs,
      COLLECTWGSMETRICS.out.metrics,
    )
    ch_bam_and_metrics = WorkflowTeal.sort_bams_and_metrics(ch_bams_and_metrics_unsorted)
    ch_cobalt_purple = ch_inputs
      .map { meta, tbam, nbam, cobalt_dir, purple_dir ->
        [meta, cobalt_dir, purple_dir]
      }
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
