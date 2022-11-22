process CUPPA_VISUALISER {
  container 'docker.io/scwatts/cuppa:1.7--0'

  input:
  tuple val(meta), path(cuppa_csv)

  output:
  path 'cuppa_chart/*png'
  path 'cuppa_chart/*txt'
  path 'versions.yml'    , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  """
  mkdir -p cuppa_chart/

  python software/cuppa_chart/cuppa-chart.py \\
    -sample ${meta.get(['sample_name', 'tumor'])} \\
    -sample_data ${cuppa_csv} \\
    -output_dir cuppa_chart/

  # NOTE(SW): hard coded since there is no reliable way to obtain version information.
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      cuppa: 1.7
  END_VERSIONS
  """

  stub:
  """
  mkdir -p cuppa_chart/
  touch cuppa_chart/${meta.get(['sample_name', 'tumor'])}.cuppa.chart.png
  touch cuppa_chart/${meta.get(['sample_name', 'tumor'])}.cuppa.conclusion.txt
  echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
  """
}

