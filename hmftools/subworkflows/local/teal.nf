//
// TEAL measures telomere content, and estimates telomeric length based on WGS read data.
//

include { PICARD_COLLECTWGSMETRICS as COLLECTWGSMETRICS } from '../../modules/nf-core/modules/picard/collectwgsmetrics/main'

include { TEAL as TEAL_PROCESS } from '../../modules/scwatts/nextflow_modules/teal/main'

workflow TEAL {
  take:
    ch_inputs_bams          // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
    ch_inputs_other         // channel: [val(meta), cobalt_dir, purple_dir]
    ref_data_genome_dir     //    file: /path/to/genome_dir/
    ref_data_genome_fn      //     val: genome name

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Obtain WGS metrics
    // NOTE(SW): here I remove duplicate files so that we only process each input once
    // channel: [val(meta_teal), bam, bai]
    ch_input_bams_split = WorkflowTeal.split_input_bams(ch_inputs_bams)
    // channel: [val(meta_teal), bam, bai]
    ch_input_bams_split_unique = WorkflowTeal.get_unique_input_files(ch_input_bams_split)
    // channel: [val(meta_teal), bam]
    ch_metrics_inputs = ch_input_bams_split_unique.map { meta_teal, bam, bai -> [meta_teal, bam] }
    COLLECTWGSMETRICS(
      ch_metrics_inputs,
      "${ref_data_genome_dir}/${ref_data_genome_fn}",
    )
    ch_versions = ch_versions.mix(COLLECTWGSMETRICS.out.versions)
    // Run TEAL
    // NOTE(SW): performed here to avoid WorkflowHmftools import in WorkflowTeal
    // channel: [val(meta_teal), wgs_metrics, bam, bai]
    ch_bams_and_metrics_unsorted = WorkflowHmftools.group_by_meta(
      ch_input_bams_split_unique,
      COLLECTWGSMETRICS.out.metrics,
    )
    // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, tumor_wgs_metrics, normal_wgs_metrics]
    ch_bam_and_metrics = WorkflowTeal.sort_bams_and_metrics(ch_bams_and_metrics_unsorted)
    // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, tumor_wgs_metrics, normal_wgs_metrics, cobalt_dir, purple_dir]
    ch_teal_inputs = WorkflowHmftools.group_by_meta(
      ch_bam_and_metrics,
      ch_inputs_other,
    )
    TEAL_PROCESS(
      ch_teal_inputs,
    )
    ch_versions = ch_versions.mix(TEAL_PROCESS.out.versions)

  emit:
    versions = ch_versions // channel: [versions.yml]
}
