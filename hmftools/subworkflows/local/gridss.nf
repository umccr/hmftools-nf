//
// GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements.
//

include { ANNOTATE          } from '../../modules/scwatts/nextflow_modules/gridss/annotate/main'
include { ASSEMBLE          } from '../../modules/scwatts/nextflow_modules/gridss/assemble/main'
include { CALL              } from '../../modules/scwatts/nextflow_modules/gridss/call/main'
include { EXTRACT_FRAGMENTS } from '../../modules/scwatts/nextflow_modules/gridss/extract_fragments/main'
include { PREPROCESS        } from '../../modules/scwatts/nextflow_modules/gridss/preprocess/main'

workflow GRIDSS {
  take:
    ch_bams_and_index         // channel: [val(meta), val(sample_name), bam, bai]
    ch_sv_vcfs                // channel: [val(meta), val(sample_name), vcf]
    ref_data_genome_dir       // file: /path/to/genome/dir/
    ref_data_genome_fn        // val: genome name
    ref_data_gridss_blacklist // val: /path/to/gridss/blacklist.bed

  main:
    // NOTE(SW): BAM indices are required by the optional process EXTRACT_FRAGMENTS. To avoid
    // introducing further forks/branchs, we just required BAM indices for all inputs
    // irrespective of EXTRACT_FRAGMENTS execution.

    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Match SV VCFs and BAMs, differentiate samples with SV VCFs
    // NOTE(SW): The .branch operator does not appear to support named unpacking in this instance
    // so we access by index; branch on presence of SV VCF path and return meta and sample name.
    ch_matched_inputs = ch_bams_and_index.join(ch_sv_vcfs, by: [0, 1], remainder: true)
    ch_sample_has_sv_vcf = ch_matched_inputs
      .branch {
        yes: it[4]
          return it[0..1]
        no: ! it[4]
          return it[0..1]
      }
    // Get channel of BAMs + SV VCFs to run EXTRACT_FRAGMENTS; join on meta and sample name
    ch_extract_fragments_input = ch_matched_inputs.join(ch_sample_has_sv_vcf.yes, by: [0, 1])

    EXTRACT_FRAGMENTS(
      ch_extract_fragments_input,
    )
    ch_versions = ch_versions.mix(EXTRACT_FRAGMENTS.out.versions)

    // Select inputs that have no SV VCFs and remove BAI and SV VCF channel slot
    ch_full_run_input = ch_matched_inputs
      .join(ch_sample_has_sv_vcf.no, by: [0, 1])
      .map { it[0..-3] }

    ch_bams = Channel.empty()
      .concat(
        EXTRACT_FRAGMENTS.out.bam,
        ch_full_run_input,
      )

    PREPROCESS(
      ch_bams,
      ref_data_genome_dir,
      ref_data_genome_fn,
    )
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)

    ch_bams_split = WorkflowGridss.split_by_sample_type(ch_bams)
    ch_preprocess_output_split = WorkflowGridss.split_by_sample_type(PREPROCESS.out.preprocess_dir)

    ch_assemble_input = WorkflowHmftools.group_by_meta(
        ch_bams_split.get('tumor'),
        ch_bams_split.get('normal'),
        ch_preprocess_output_split.get('tumor'),
        ch_preprocess_output_split.get('normal'),
    )

    ASSEMBLE(
      ch_assemble_input,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_gridss_blacklist,
    )
    ch_versions = ch_versions.mix(ASSEMBLE.out.versions)

    ch_call_input = WorkflowHmftools.group_by_meta(
      ch_bams_split.get('tumor'),
      ch_bams_split.get('normal'),
      ASSEMBLE.out.assembly_dir,
    )

    CALL(
      ch_call_input,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_gridss_blacklist,
    )
    ch_versions = ch_versions.mix(CALL.out.versions)

    ch_annotate_input = CALL.out.vcf.filter { meta, vcf ->
        return WorkflowHmftools.has_records_vcf(vcf)
      }

    // Annotate with RepeatMasker, required for LINX
    ANNOTATE(ch_annotate_input)
    ch_versions = ch_versions.mix(ANNOTATE.out.versions)

  emit:
    results  = ANNOTATE.out.vcf // channel: [val(meta), vcf]

    versions = ch_versions      // channel: [versions.yml]
}
