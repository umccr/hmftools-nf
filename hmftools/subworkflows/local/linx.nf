//
// Linx is an annotation, interpretation and visualisation tool for structural variants.
//

include { ANNOTATION } from '../../modules/scwatts/nextflow_modules/linx/annotation/main'
include { VISUALISER } from '../../modules/scwatts/nextflow_modules/linx/visualiser/main'

workflow LINX {
  take:
    ch_purple_dir               // channel: [meta, purple]
    ref_data_linx_fragile_sites // file: /path/to/linx/fragile/sites.csv
    ref_data_linx_line_elements // file: /path/to/linx/lines.csv
    ref_data_ensembl_data_dir   // file: /path/to/hmf/ensembl/data/dir/
    ref_data_known_fusion_data  // file: /path/to/hmf/known/fusion/data.csv
    ref_data_driver_gene_panel  // file: /path/to/hmf/driver/gene/panel.tsv

  main:
    // Channel for versions.yml files
    ch_versions = Channel.empty()

    ANNOTATION(
      ch_purple_dir,
      ref_data_linx_fragile_sites,
      ref_data_linx_line_elements,
      ref_data_ensembl_data_dir,
      ref_data_known_fusion_data,
      ref_data_driver_gene_panel,
    )
    ch_versions = ch_versions.mix(ANNOTATION.out.versions)

    VISUALISER(
      ANNOTATION.out.annotation_dir,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(VISUALISER.out.versions)

  ch_linx_out = WorkflowHmftools.group_by_meta(
    ANNOTATION.out.annotation_dir,
    VISUALISER.out.visualiser_dir,
  )

  emit:
    results = ch_linx_out  // channel: [meta, linx_annotation, linx_visualiser]

    versions = ch_versions // channel: [versions.yml]
}
