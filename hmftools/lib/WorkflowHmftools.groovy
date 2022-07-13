//
// This file holds several functions specific to the workflows/hmftools.nf in the nf-core/hmftools pipeline
//

import static groovy.io.FileType.FILES

import nextflow.Channel
import nextflow.Nextflow

class WorkflowHmftools {

  //
  // Check and validate parameters
  //
  public static void initialise(params, workflow, log) {
    // Download region file for collectwgsmetrics if in test mode
    if (workflow.profile.contains('test')) {
        def work_dir = new File(workflow.workDir.toString())
        def interval_file = File.createTempFile('collectwgsmetrics.interval_list', null, work_dir)
        interval_file.deleteOnExit()
        interval_file << new URL (params.ref_data_wgsmetrics_intervals_url).getText()
        params.ref_data_wgsmetrics_intervals_local = interval_file
    }
  }

  public static List set_stages(mode, log) {
    def stages = []
    switch(mode) {
      case 'full':
        stages = Stage.values()
        break
      case 'gridss-purple-linx':
        stages = [
          Stage.AMBER,
          Stage.COBALT,
          Stage.GRIDSS,
          Stage.GRIPSS,
          Stage.PURPLE,
          Stage.LINX,
        ]
        break
      case 'gridss':
        stages = [
          Stage.GRIDSS,
          Stage.GRIPSS,
        ]
        break
      case 'purple':
        stages = [
          Stage.PURPLE,
        ]
        break
      case 'linx':
        stages = [
          Stage.LINX,
        ]
        break
      case 'lilac':
        stages = [
          Stage.LILAC
        ]
        break
      case 'teal':
        stages = [
          Stage.TEAL
        ]
        break
      default:
        log.error "\nERROR: recieved invalid mode '${mode}'"
        System.exit(1)
    }
    return stages
  }

  enum Stage {
    AMBER,
    COBALT,
    GRIDSS,
    GRIPSS,
    LILAC,
    LINX,
    PAVE,
    PURPLE,
    SAGE,
    TEAL,
  }

  public static prepare_inputs(ch, stub_run, log) {
    return ch
      .splitCsv(header: true, strip: true, sep: '\t')
      .map { [it.id, it] }
      .groupTuple()
      .map { id, inputs ->
        def meta = [
          'id': id,
        ]
        inputs.each {
          assert meta.id == it.id
          // Add subject name if not already present
          if (meta.containsKey('subject_name')) {
              assert meta.subject_name == it.subject_name
          } else {
              meta.subject_name = it.subject_name
          }

          // Set sample name
          def key = []
          if (it.sample_type == 'tumor') {
            key = ['sample_name', 'tumor']
          } else if (it.sample_type == 'normal') {
            key = ['sample_name', 'normal']
          } else {
            assert false
          }
          if (meta.containsKey(key)) {
            assert meta[key] == it.sample_name
          } else {
            meta[key] = it.sample_name
          }

          // Add file
          def key_file = [it.filetype, it.sample_type]
          assert ! meta.containsKey(key_file)
          meta[key_file] = it.filepath

          if (! stub_run) {
            // For BAM file inputs, require co-located index
            if (it.filepath.endsWith('.bam')) {
              def bam_index_fp_str = "${it.filepath}.bai".toString()
              def bam_index_fp = Nextflow.file(bam_index_fp_str)
              if (! bam_index_fp.exists()) {
                log.error "\nERROR: No index found for ${it.filepath}"
                System.exit(1)
              }
            }

            // For GRIPSS SV VCFs, require co-located index
            if (it.filetype.startsWith('gripss')) {
              def vcf_index_fp = Nextflow.file("${it.filepath}.tbi".toString())
              if (! vcf_index_fp.exists()) {
                log.error "\nERROR: No index found for ${it.filepath}"
                System.exit(1)
              }
            }
          }

          // NOTE(SW): CHECK_SAMPLESHEET curently enforces T/N; this may be relevant in the future
          //// Set sample type: tumor_normal, tumor_only, normal_only
          //if (meta.containsKey(['sample_name', 'tumor']) && meta.containsKey(['sample_name', 'normal'])) {
          //  meta.sample_type = 'tumor_normal'
          //} else if (meta.containsKey(['sample_name', 'tumor'])) {
          //  meta.sample_type = 'tumor_only'
          //} else if (meta.containsKey(['sample_name', 'normal'])) {
          //  meta.sample_type = 'normal_only'
          //} else {
          //  assert false
          //}
        }

        // For PURPLE only runs, we must get normal sample name from inputs since there is no way to provide this
        // in the samplesheet
        if (meta.containsKey(['cobalt_dir', 'tumor']) && ! meta.containsKey(['sample_name', 'normal'])) {
          // Discover files
          def normal_ratio_fps = []
          new File(meta[['cobalt_dir', 'tumor']])
            .eachFileMatch(groovy.io.FileType.FILES, ~/.+\.cobalt\.gc\.median\.tsv/, { normal_ratio_fps << it })
          // Select normal sample file
          def normal_ratio_fp = normal_ratio_fps
            .findAll { ! it.getName().contains(meta[['sample_name', 'tumor']]) }
          assert normal_ratio_fp.size() == 1
          // Set normal sample name
          def m = (normal_ratio_fp =~ /.+\/(.+)\.cobalt\.gc\.median\.tsv/)
          meta[['sample_name', 'normal']] = m[0][1]
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
}
