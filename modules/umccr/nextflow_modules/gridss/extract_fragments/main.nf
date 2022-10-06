process EXTRACT_FRAGMENTS {
  //conda (params.enable_conda ? "bioconda::gridss=2.13.2" : null)
  container 'docker.io/scwatts/gridss:2.13.2--3'

  input:
  tuple val(meta), path(bam), path(bai), path(sv_vcf)

  output:
  tuple val(meta), path("${output_fp}"), emit: bam
  path 'versions.yml'                  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  // NOTE(SW): use of global scoping for output_fp to allow access in the output directive above
  output_fp = "gridss_extract_fragments/${bam.getSimpleName()}.targeted.bam"

  """
  # Run
  gridss_extract_overlapping_fragments \
    ${args} \
    --jar "${task.ext.jarPath}" \
    --targetvcf "${sv_vcf}" \
    --workingdir gridss_extract_fragments/work/ \
    --output "${output_fp}" \
    --threads "${task.cpus}" \
    "${bam}"
  # This script can exit silently, check that we have some reads in the output file before proceeding
  if [[ "\$(samtools view "${output_fp}" | head | wc -l)" -eq 0 ]]; then
    exit 1;
  fi;

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      gridss: \$(java -cp "${task.ext.jarPath}" gridss.CallVariants --version 2>&1 | sed 's/-gridss//')
  END_VERSIONS
  """

  stub:
  output_fp = "gridss_extract_fragments/${bam.getSimpleName()}.targeted.bam"

  """
  mkdir -p gridss_extract_fragments/
  touch "${output_fp}"
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
