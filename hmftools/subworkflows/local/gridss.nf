//
// GRIDSS is a software suite containing tools useful for the detection of genomic rearrangements.
//

include { ANNOTATE          } from '../../modules/scwatts/nextflow_modules/gridss/annotate/main'
include { ASSEMBLE          } from '../../modules/scwatts/nextflow_modules/gridss/assemble/main'
include { CALL              } from '../../modules/scwatts/nextflow_modules/gridss/call/main'
include { EXTRACT_FRAGMENTS } from '../../modules/scwatts/nextflow_modules/gridss/extract_fragments/main'
include { PREPROCESS        } from '../../modules/scwatts/nextflow_modules/gridss/preprocess/main'

workflow GRIDSS {
  take:
    ch_meta                   // channel: [val(meta)]
    gridss_config
    ref_data_genome_dir       // file: /path/to/genome/dir/
    ref_data_genome_fn        // val: genome name
    ref_data_gridss_blacklist // val: /path/to/gridss/blacklist.bed

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Prepare inputs
    ch_inputs = WorkflowGridss.get_inputs(ch_meta)
      // NOTE(SW): calling .branch within WorkflowGridss.get_inputs_from_meta triggers unavoidable
      // errors. Hence, the .branch is called at this scope.
      .branch {
        extract_fragments: it[0]
          return it[1]
        preprocess: ! it[0]
          return it[1]
      }
    // In instances where the same input BAM is provided (e.g. multiple t/n pairs sharing the same normal), select
    // only one file to process
    ch_extract_fragments_inputs = WorkflowGridss.get_unique_input_files(ch_inputs.extract_fragments)
    ch_preprocess_inputs_raw = WorkflowGridss.get_unique_input_files(ch_inputs.preprocess)

    // Read selection using prior SV calls
    EXTRACT_FRAGMENTS(
      ch_extract_fragments_inputs,
    )
    ch_versions = ch_versions.mix(EXTRACT_FRAGMENTS.out.versions)

    // Preprocess reads
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
    // NOTE(SW): performed here to avoid WorkflowHmftools import in WorkflowGridss
    ch_bams_and_preprocess = WorkflowHmftools.group_by_meta(
      ch_preprocess_inputs,
      PREPROCESS.out.preprocess_dir,
    )
      .map { [it[0].subject_name, it] }
      .groupTuple()
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
    ch_call_inputs = WorkflowHmftools.group_by_meta(
      ch_assemble_inputs,
      ASSEMBLE.out.assemble_dir,
      flatten: false
    )
      .map { meta, other ->
        def (bams, preprocess_dirs, labels) = other[0]
        def (assemble_dir) = other[1]
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

    // Annotate with RepeatMasker, required for LINX
    ch_annotate_inputs = CALL.out.vcf.filter { meta, vcf ->
        return WorkflowHmftools.has_records_vcf(vcf)
      }
    ANNOTATE(ch_annotate_inputs)
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)

    // Pair annotated output with input meta
    ch_out = Channel.empty()
      .concat(
        ch_meta.map { meta -> [meta.id, meta] },
        ANNOTATE.out.vcf.flatMap { meta_gridss, sv ->
          meta_gridss.id_list.collect { id -> [id, sv] }
        }
      )
      .groupTuple()
      .map { id, other -> other.flatten() }

  emit:
    results  = ch_out      // channel: [val(meta), vcf]

    versions = ch_versions // channel: [versions.yml]
}
