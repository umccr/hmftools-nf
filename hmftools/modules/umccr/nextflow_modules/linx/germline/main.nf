process LINX_GERMLINE {
  //conda (params.enable_conda ? "bioconda::hmftools-linx=1.21" : null)
  container 'quay.io/biocontainers/hmftools-linx:1.21--hdfd78af_0'

  input:
  tuple val(meta), path(gripss_sv)
  val ref_data_genome_ver
  path fragile_sites
  path lines
  path ensembl_data_dir
  path driver_gene_panel

  output:
  tuple val(meta), path('linx_germline/'), emit: annotation_dir
  path 'versions.yml'                    , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  java \
    -Xmx${task.memory.giga}g \
    -jar "${task.ext.jarPath}" \
      ${args} \
      -sample "${meta.get(['sample_name', 'normal'])}" \
      -germline \
      -ref_genome_version "${ref_data_genome_ver}" \
      -sv_vcf "${gripss_sv}" \
      -fragile_site_file "${fragile_sites}" \
      -line_element_file "${lines}" \
      -ensembl_data_dir "${ensembl_data_dir}" \
      -check_drivers \
      -driver_gene_panel "${driver_gene_panel}" \
      -output_dir linx_germline/

  # NOTE(SW): hard coded since there is no reliable way to obtain version information
  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      linx: 1.19
  END_VERSIONS
  """

  stub:
  """
  mkdir linx_germline/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
