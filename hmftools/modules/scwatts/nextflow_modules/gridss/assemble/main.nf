process ASSEMBLE {
  //conda (params.enable_conda ? "bioconda::gridss=2.13.2" : null)
  container 'docker.io/scwatts/gridss:2.13.2'

  input:
  tuple val(meta), path(bams), path(preprocess_dirs), val(labels)
  path gridss_config
  path ref_data_genome_dir
  val ref_data_genome_fn
  path blacklist

  output:
  tuple val(meta), path('gridss_assemble/'), emit: assemble_dir
  path 'versions.yml'                      , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def config_arg = gridss_config ? "--configuration ${gridss_config}" : ''
  def output_dirname = 'gridss_assemble'
  def labels_arg = labels.join(',')
  // NOTE(SW): Nextflow implicitly casts List<TaskPath> to an atomic TaskPath, hence the required check below
  def bams_list = bams instanceof List ? bams : [bams]
  def bams_arg = bams_list.join(' ')

  """
  # Create shadow directory with file symlinks of GRIDSS 'working' dir to prevent cache invalidation
  # NOTE: for reasons that elude me, NF doesn't always stage in the workingdir; remove if it is present
  shadow_input_directory() {
    src=\${1}
    dst="${output_dirname}/work/\${src##*/}"
    for filepath_src in \$(find -L \${src} ! -type d); do
      # Get destination location for symlink
      filepath_src_rel=\$(sed 's#^'\${src}'/*##' <<< \${filepath_src})
      filepath_dst=\${dst%/}/\${filepath_src_rel}
      # Create directory for symlink
      mkdir -p \${filepath_dst%/*};
      # Get path for symlink source file, then create it
      # NOTE(SW): ideally we would get the relative path using the --relative-to but this is only
      # supported for GNU realpath and fails for others such as BusyBox, which is used in Biocontainers
      symlinkpath=\$(realpath \${filepath_src})
      ln -s \${symlinkpath} \${filepath_dst};
    done
    if [[ -L "\${src##*/}" ]]; then
      rm "\${src}"
    fi
  }
  for preprocess_dir in ${preprocess_dirs}; do
    shadow_input_directory \${preprocess_dir};
  done

  # Run
  gridss \
    ${args} \
    --jvmheap "${task.memory.giga}g" \
    --jar "${task.ext.jarPath}" \
    --steps assemble \
    --labels "${labels_arg}" \
    --reference "${ref_data_genome_dir}/${ref_data_genome_fn}" \
    --blacklist "${blacklist}" \
    --workingdir "${output_dirname}/work" \
    --assembly "${output_dirname}/sv_assemblies.bam" \
    --threads "${task.cpus}" \
    ${config_arg} \
    ${bams_arg}

  # NOTE(SW): hard coded since there is no reliable way to obtain version information, see GH issue
  # https://github.com/PapenfussLab/gridss/issues/586
  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      gridss: 2.13.2
  END_VERSIONS
  """

  stub:
  """
  mkdir -p gridss_assemble/
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}
