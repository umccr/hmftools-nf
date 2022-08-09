//
// GRIDSS is a software suite containing tools useful for the detection of genomic rearrangements.
//

include { ASSEMBLE          } from '../../modules/scwatts/nextflow_modules/gridss/assemble/main'
include { CALL              } from '../../modules/scwatts/nextflow_modules/gridss/call/main'
include { EXTRACT_FRAGMENTS } from '../../modules/scwatts/nextflow_modules/gridss/extract_fragments/main'
include { PREPROCESS        } from '../../modules/scwatts/nextflow_modules/gridss/preprocess/main'

workflow GRIDSS {
  take:
    ch_meta                   // channel: [val(meta)]
    gridss_config             //    file: /path/to/gridss_config (optional)
    ref_data_genome_dir       //    file: /path/to/genome_dir/
    ref_data_genome_fn        //     val: genome name
    ref_data_gridss_blacklist //     val: /path/to/gridss_blacklist

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Prepare inputs
    ch_inputs = WorkflowGridss.get_inputs(ch_meta)
      // NOTE(SW): calling .branch within WorkflowGridss.get_inputs_from_meta triggers unavoidable
      // errors. Hence, the .branch is called at this scope.
      .branch {
        // channel: [val(meta_gridss), bam, bai, vcf]
        extract_fragments: it[0]
          return it[1]
        // channel: [val(meta_gridss), bam]
        preprocess: ! it[0]
          return it[1]
      }
    // In instances where the same input BAM is provided (e.g. multiple t/n pairs sharing the same normal), select
    // only one file to process
    // channel: [val(meta_gridss), bam, bai, vcf]
    ch_extract_fragments_inputs = WorkflowGridss.get_unique_input_files(ch_inputs.extract_fragments)
    // channel: [val(meta_gridss), bam]
    ch_preprocess_inputs_raw = WorkflowGridss.get_unique_input_files(ch_inputs.preprocess)

    // Read selection using prior SV calls
    EXTRACT_FRAGMENTS(
      ch_extract_fragments_inputs,
    )
    ch_versions = ch_versions.mix(EXTRACT_FRAGMENTS.out.versions)

    // Preprocess reads
    // channel: [val(meta_gridss), bam]
    ch_preprocess_inputs = Channel
      .empty()
      .concat(
        EXTRACT_FRAGMENTS.out.bam,
        ch_preprocess_inputs_raw,
      )
    PREPROCESS(
      ch_preprocess_inputs,
      gridss_config,
      ref_data_genome_dir,
      ref_data_genome_fn,
    )
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)

    // Joint assembly
    // The data flow in this following section is complex. For joint assembly, we must collect
    // all subject BAMs and preprocess output directories but to do this in a non-blocking way,
    // the number expected files/entries for each subject must be provided.

    // First, we obtain the number of entries of each subject directly from the input meta. This
    // evaluates immediately.
    // channel: [subject_name, count]
    ch_subject_files_expected = Channel.empty()
      .mix(
        ch_extract_fragments_inputs.map { it[0].subject_name },
        ch_preprocess_inputs_raw.map { it[0].subject_name },
      )
      .collect()
      .flatMap { names ->
        names
          .countBy { it }
          .collect { k, v -> [k, v] }
      }

    // Next, we prepare a channel containing the required files to join with the above counts.
    // NOTE(SW): performed here to avoid WorkflowHmftools import in WorkflowGridss
    // channel: [subject_name, [[val(meta_gridss), bam, preprocess_dir]]
    ch_assembly_inputs_base = WorkflowHmftools.group_by_meta(
      ch_preprocess_inputs,
      PREPROCESS.out.preprocess_dir,
    )
      .map { [it[0].subject_name, it] }

    // Finally, we create the basis for the assembly input channel. Each element in ch_subject_files_expected
    // is associated with the corresponding element in ch_subject_files_expected using the cross
    // operator, creating a channel with the format
    //   [[subject_name, count], [subject_name, [val(meta_gridss), bam, preprocess_dir]]]
    // From this data, a groupKey with the expected size emit groupTuple size is created.
    // channel: [subject_name, [[val(meta_gridss), bam, preprocess_dir], ...]]
    ch_bams_and_preprocess = ch_subject_files_expected
      .cross(ch_assembly_inputs_base)
      .map { bams_expected, data ->
        def (subject_name, count) = bams_expected
        def values = data[1]
        return [groupKey(subject_name, count), values]
      }
      .groupTuple()

    // channel: [val(meta_gridss), [bams], [preprocess_dirs], [labels]]
    ch_assemble_inputs = WorkflowGridss.get_assemble_inputs(ch_bams_and_preprocess)
    ASSEMBLE(
      ch_assemble_inputs,
      gridss_config,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_gridss_blacklist,
    )
    ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

    // Joint calling
    // channel: [val(meta_gridss), [bams], assemble_dir, [labels]]
    ch_call_inputs = WorkflowHmftools.group_by_meta(
      ch_assemble_inputs,
      ASSEMBLE.out.assemble_dir,
      flatten: false,
    )
      .map { data ->
        def meta = data[0]
        def (bams, preprocess_dirs, labels) = data[1]
        def (assemble_dir) = data[2]
        return [meta, bams, assemble_dir, labels]
      }
    CALL(
      ch_call_inputs,
      gridss_config,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_gridss_blacklist,
    )
    ch_versions = ch_versions.mix(CALL.out.versions)

    // Pair output VCF with input meta
    // channel: [val(meta), vcf]
    ch_out = Channel.empty()
      .concat(
        ch_meta.map { meta -> [meta.id, meta] },
        CALL.out.vcf.flatMap { meta_gridss, sv ->
          meta_gridss.id_list.collect { id -> [id, sv] }
        }
      )
      .groupTuple(size: 2)
      .map { id, other -> other.flatten() }

  emit:
    results  = ch_out      // channel: [val(meta), vcf]

    versions = ch_versions // channel: [versions.yml]
}
