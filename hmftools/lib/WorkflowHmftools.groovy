//
// This file holds several functions specific to the workflow/hmftools.nf in the nf-core/hmftools pipeline
//

import nextflow.Channel

class WorkflowHmftools {

  //
  // Check and validate parameters
  //
  public static void initialise(params, log) {
  }

  public static prepare_inputs(samplesheet_fp, stub_run, log) {
    return Channel.fromPath(samplesheet_fp)
      .splitCsv(header: true, strip: true, sep: '\t')
      .map { [it.subject_name, it] }
      .groupTuple()
      .map { subject_name, inputs ->
        def meta = [
          'subject_name': subject_name,
        ]
        inputs.each {
          assert meta.subject_name == it.subject_name
          // Add sample name
          def key_sample_name = ['sample_name', it.sample_type]
          if (meta.containsKey(key_sample_name)) {
            assert meta.get(key_sample_name) == it.sample_name
          } else {
            meta[key_sample_name] = it.sample_name
          }
          // Add file
          def key_file = [it.sample_type, it.filetype]
          assert ! meta.containsKey(key_file)
          meta[key_file] = it.filepath
          // For BAM file inputs, required that it has a co-located index; ignore for stub runs
          if (it.filepath.endsWith('.bam') && ! stub_run) {
            def bam_index_fp = new File("${it.filepath}.bai")
            if (! bam_index_fp.exists()) {
              log.error "\nERROR: No index found for ${it.filepath}"
              System.exit(1)
            }
          }
        }
        return meta
      }
  }

  public static group_by_meta(Map named_args, ... channels) {
    def r = channels
    if (! named_args.get('interleave', false)) {
      r = r.collect { ch -> ch.map { [it[0], it[1..-1]] }}
    }
    r = Channel.empty()
      .concat(
        *r,
      )
    def tuple_size = named_args.get('tuple_size')
    if (tuple_size) {
      r = r.groupTuple(size: tuple_size)
    } else {
      r = r.groupTuple()
    }
    if (named_args.get('flatten', true)) {
      r = r.map { it.flatten() }
    }
    return r
  }

  // NOTE(SW): function signature required to catch where no named arguments are passed
  public static group_by_meta(... channels) {
    return group_by_meta([:], *channels)
  }

  public static Boolean has_records_vcf(filepath) {
    def command = "bcftools view -H ${filepath} | head | wc -l"
    def results = execute_command(command)
    def stdout = results[1]
    return stdout.toInteger() > 0
  }

  public static <T> List<T> execute_command(command) {
    def command_fq = ['/bin/bash', '-c', command]
    def stdout = new StringBuilder()
    def stderr = new StringBuilder()
    def process = command_fq.execute()
    process.waitForProcessOutput(stdout, stderr)
    return [process.exitValue(), stdout, stderr]
  }

  public static get_gridss_reads_input(ch) {
    return get_gridss_input(ch, 'bam').map { meta, sample_name, filepath ->
      [meta, sample_name, filepath, filepath + '.bai']
    }
  }

  public static get_gridss_svs_input(ch) {
    return get_gridss_input(ch, 'sv')
  }

  public static get_gridss_input(ch, filetype) {
    // Gets filepath if it exists otherwise Map.get returns null, and channel items with null
    // filepaths are subsequently filtered.
    def sample_types = ['tumor', 'normal']
    return ch.flatMap { meta ->
      sample_types.collect { sample_type ->
        [
          meta,
          meta.get(['sample_name', sample_type]),
          meta.get([sample_type, filetype]),
        ]
      }
    }.filter { meta, sample_name, filepath -> filepath }
  }
}
