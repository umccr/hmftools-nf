//
// Linx is an annotation, interpretation and visualisation tool for structural variants.
//

include { LINX_GERMLINE } from '../../modules/scwatts/nextflow_modules/linx/germline/main'
include { LINX_SOMATIC  } from '../../modules/scwatts/nextflow_modules/linx/somatic/main'
include { VISUALISER    } from '../../modules/scwatts/nextflow_modules/linx/visualiser/main'

workflow LINX {
  take:
    ch_linx_germline_inputs     // channel: [meta, purple]
    ch_linx_somatic_inputs      // channel: [meta, gripss_sv]
    ref_data_linx_fragile_sites // file: /path/to/linx/fragile/sites.csv
    ref_data_linx_line_elements // file: /path/to/linx/lines.csv
    ref_data_ensembl_data_dir   // file: /path/to/hmf/ensembl/data/dir/
    ref_data_known_fusion_data  // file: /path/to/hmf/known/fusion/data.csv
    ref_data_driver_gene_panel  // file: /path/to/hmf/driver/gene/panel.tsv

  main:
    // Channel for versions.yml files
    ch_versions = Channel.empty()

    // Germline
    LINX_GERMLINE(
      ch_linx_germline_inputs,
      ref_data_linx_fragile_sites,
      ref_data_linx_line_elements,
      ref_data_ensembl_data_dir,
      ref_data_driver_gene_panel,
    )
    ch_versions = ch_versions.mix(LINX_GERMLINE.out.versions)

    // Somatic
    LINX_SOMATIC(
      ch_linx_somatic_inputs,
      ref_data_linx_fragile_sites,
      ref_data_linx_line_elements,
      ref_data_ensembl_data_dir,
      ref_data_known_fusion_data,
      ref_data_driver_gene_panel,
    )
    ch_versions = ch_versions.mix(LINX_SOMATIC.out.versions)

    VISUALISER(
      LINX_SOMATIC.out.annotation_dir,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(VISUALISER.out.versions)

    ch_linx_somatic_out = WorkflowHmftools.group_by_meta(
      LINX_SOMATIC.out.annotation_dir,
      VISUALISER.out.visualiser_dir,
    )

  emit:
    somatic  = ch_linx_somatic_out              // channel: [meta, linx_annotation, linx_visualiser]
    germline = LINX_GERMLINE.out.annotation_dir // channel: [meta, linx_annotation]

    versions = ch_versions                      // channel: [versions.yml]
}
