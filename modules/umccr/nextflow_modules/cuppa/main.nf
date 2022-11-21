process CUPPA {
  container 'docker.io/scwatts/cuppa:1.7--0'

  input:
  tuple val(meta), path(isofox_dir), path(purple_dir), path(linx_dir), path(virus_interpreter_dir)
  path reference_data

  output:
  tuple val(meta), path('cuppa/'), emit: cuppa_dir
  path 'versions.yml'            , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  # Symlink input files into a single directory
  mkdir sample_data/
  input_dirs="${isofox_dir} ${purple_dir} ${linx_dir} ${virus_interpreter_dir}"
  find ${input_dirs} -type f -exec ln -s ../{} sample_data/ \\;

  java \\
    -Xmx${task.memory.giga}g \\
    -jar ${task.ext.jarPath} \\
      -categories ALL \\
      -ref_data_dir ${reference_data} \\
      -sample_data ${meta.id} \\
      -sample_data_dir sample_data/ \\
      -output_dir cuppa/

  # NOTE(SW): hard coded since there is no reliable way to obtain version information.
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      cuppa: 1.7
  END_VERSIONS
  """

  stub:
  """
  mkdir -p cuppa/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}

