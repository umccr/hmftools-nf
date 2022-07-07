//
// GRIPSS performs SV filtering.
//

include { GRIPSS_GERMLINE } from '../../modules/scwatts/nextflow_modules/gripss/germline/main'
include { GRIPSS_SOMATIC  } from '../../modules/scwatts/nextflow_modules/gripss/somatic/main'

workflow GRIPSS {
  take:
    ch_inputs           // channel: [meta, gridss_vcf]
    ref_data_genome_dir
    ref_data_genome_fn
    breakend_pon
    breakpoint_pon
    known_fusions

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Germline
    GRIPSS_GERMLINE(
      ch_inputs,
      ref_data_genome_dir,
      ref_data_genome_fn,
      breakend_pon,
      breakpoint_pon,
      known_fusions,
    )
    ch_versions = ch_versions.mix(GRIPSS_GERMLINE.out.versions)

    // Somatic
    GRIPSS_SOMATIC(
      ch_inputs,
      ref_data_genome_dir,
      ref_data_genome_fn,
      breakend_pon,
      breakpoint_pon,
      known_fusions,
    )
    ch_versions = ch_versions.mix(GRIPSS_SOMATIC.out.versions)

  // Pack output
  ch_germline_out = WorkflowHmftools.group_by_meta(
    GRIPSS_GERMLINE.out.vcf_hard,
    GRIPSS_GERMLINE.out.vcf_soft,
  )
  ch_somatic_out = WorkflowHmftools.group_by_meta(
    GRIPSS_SOMATIC.out.vcf_hard,
    GRIPSS_SOMATIC.out.vcf_soft,
  )

  emit:
    germline = ch_germline_out
    somatic  = ch_somatic_out

    versions = ch_versions // channel: [versions.yml]
}
