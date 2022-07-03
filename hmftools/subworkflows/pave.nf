//
// PAVE annotates somatic and germline variant VCFs with gene and transcript coding and protein effects.
//

include { PAVE_GERMLINE } from '../modules/scwatts/nextflow_modules/pave/germline/main'
include { PAVE_SOMATIC  } from '../modules/scwatts/nextflow_modules/pave/somatic/main'

workflow PAVE {
  take:
    ch_inputs_germline                   // channel: [meta, sage_germline_vcf]
    ch_inputs_somatic                    // channel: [meta, sage_somatic_vcf]
    ref_data_genome_dir
    ref_data_genome_fn
    ref_data_sage_pon_file
    ref_data_sage_blacklist_bed
    ref_data_sage_blacklist_vcf
    ref_data_clinvar_vcf
    ref_data_mappability_bed
    ref_data_driver_gene_panel
    ref_data_ensembl_data_dir

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Germline
    PAVE_GERMLINE(
      ch_inputs_germline,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_sage_blacklist_bed,
      ref_data_sage_blacklist_vcf,
      ref_data_clinvar_vcf,
      ref_data_mappability_bed,
      ref_data_driver_gene_panel,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(PAVE_GERMLINE.out.versions)

    // Somatic
    PAVE_SOMATIC(
      ch_inputs_somatic,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_sage_pon_file,
      ref_data_mappability_bed,
      ref_data_driver_gene_panel,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(PAVE_SOMATIC.out.versions)

  emit:
    germline = PAVE_GERMLINE.out.vcf
    somatic = PAVE_SOMATIC.out.vcf

    versions = ch_versions // channel: [versions.yml]
}
