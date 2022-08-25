process SAMBAMBA_SLICE {
  //conda (params.enable_conda ? "bioconda::sambamba=0.8.2" : null)
  container 'quay.io/biocontainers/sambamba:0.8.2--h98b6b92_2'

  input:
  tuple val(meta), path(bam), path(bai), path(bed), val(regions)

  output:
  tuple val(meta), path("*sliced.bam"), path("*sliced.bam.bai"), emit: bam
  path 'versions.yml'                                          , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def bed_arg = bed ? "--regions ${bed}" : ''
  def regions_arg = regions ? regions.join(' ') : ''

  """
  sambamba slice \
    ${bed_arg} \
    --output-filename "${bam.simpleName}.sliced.bam" \
    ${bam} \
    ${regions_arg} \

  sambamba index "${bam.simpleName}.sliced.bam"

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      sambamba: \$(sambamba --version 2>&1 | sed -n '/sambamba/ s/^sambamba \\(.\\+\\)/\\1/p' | head -n1)
  END_VERSIONS
  """

  stub:
  """
  touch ${bam.simpleName}.sliced.bam ${bam.simpleName}.sliced.bam.bai
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
