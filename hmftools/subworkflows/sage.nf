//
// SAGE is a precise and highly sensitive somatic SNV, MNV and small INDEL caller.
//

include { SAGE_GERMLINE } from '../modules/scwatts/nextflow_modules/sage/germline/main'
include { SAGE_SOMATIC  } from '../modules/scwatts/nextflow_modules/sage/somatic/main'

workflow SAGE {
  take:
    ch_inputs                   // channel: [meta, tumor_bam, normal_bam, tumor_bai, normal_bai]
    ref_data_genome_dir
    ref_data_genome_fn
    ref_data_sage_known_hotspots_germline
    ref_data_sage_known_hotspots_somatic
    ref_data_sage_coding_panel_germline
    ref_data_sage_coding_panel_somatic
    ref_data_sage_high_confidence
    ref_data_sage_pon_file
    ref_data_mappability_bed
    ref_data_driver_gene_panel
    ref_data_ensembl_data_dir

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Germline
    SAGE_GERMLINE(
      ch_inputs,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_sage_known_hotspots_germline,
      ref_data_sage_coding_panel_germline,
      ref_data_sage_high_confidence,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(SAGE_GERMLINE.out.versions)

    // Somatic
    SAGE_SOMATIC(
      ch_inputs,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_sage_known_hotspots_somatic,
      ref_data_sage_coding_panel_somatic,
      ref_data_sage_high_confidence,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(SAGE_SOMATIC.out.versions)

  emit:
    germline = SAGE_GERMLINE.out.vcf
    somatic = SAGE_SOMATIC.out.vcf

    versions = ch_versions // channel: [versions.yml]
}
