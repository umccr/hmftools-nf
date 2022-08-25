//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

include { EXTRACT_AND_INDEX_CONTIG              } from '../../modules/umccr/nextflow_modules/custom/extract_and_index_contig/main'
include { LILAC as LILAC_PROCESS                } from '../../modules/umccr/nextflow_modules/lilac/main'
include { REALIGN_READS_LILAC as REALIGN_READS  } from '../../modules/umccr/nextflow_modules/custom/realign_reads_lilac/main'
include { SAMBAMBA_SLICE                        } from '../../modules/umccr/nextflow_modules/sambamba/slice/main'

workflow LILAC {
  take:
    ch_inputs_bams              // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
    ch_inputs_purple_dir        // channel: [val(meta), purple_dir]
    ref_data_genome_dir         //    file: /path/to/genome_dir/
    ref_data_genome_fn          //     val: genome name
    ref_data_genome_version     //     val: genome version
    ref_data_lilac_resource_dir //    file: /path/to/lilac_resource_dir/

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Slice HLA region
    // NOTE(SW): here I remove duplicate files so that we only process each input once
    // NOTE(SW): orphaned reads are sometimes obtained, this is the slicing procedure used
    // in Pipeline5, see LilacBamSlicer.java#L115
    // channel: [val(meta_lilac), bam, bai, bed, regions_list]
    ch_slice_inputs = WorkflowLilac.get_slice_inputs(ch_inputs_bams, "${ref_data_lilac_resource_dir}/hla.${ref_data_genome_version}.bed")
    // channel: [val(meta_lilac), bam, bai, bed, regions_list]
    ch_slice_inputs = WorkflowLilac.get_unique_input_files(ch_slice_inputs)
    SAMBAMBA_SLICE(
      ch_slice_inputs,
    )
    ch_versions = ch_versions.mix(SAMBAMBA_SLICE.out.versions)

    if (ref_data_genome_version == '38') {
      // Align reads with chr6
      // NOTE(SW): the aim of this process is to take reads mapping to ALT contigs and align them
      // to the three relevant HLA genes on chr6. All reads including those previously mapped to chr6
      // are realigned for consistency.
      EXTRACT_AND_INDEX_CONTIG(
        'chr6',
        ref_data_genome_dir,
        ref_data_genome_fn,
      )
      REALIGN_READS(
        SAMBAMBA_SLICE.out.bam,
        EXTRACT_AND_INDEX_CONTIG.out.contig,
        EXTRACT_AND_INDEX_CONTIG.out.bwa_indices,
      )

      // Create input channel for LILAC
      // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, purple_dir]
      ch_lilac_inputs = WorkflowHmftools.group_by_meta(
        WorkflowLilac.sort_slices(REALIGN_READS.out.bam),
        ch_inputs_purple_dir,
      )
    } else {
      // Create input channel for LILAC
      // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, purple_dir]
      ch_lilac_inputs = WorkflowHmftools.group_by_meta(
        WorkflowLilac.sort_slices(SAMBAMBA_SLICE.out.bam),
        ch_inputs_purple_dir,
      )
    }

    // Run LILAC
    LILAC_PROCESS(
      ch_lilac_inputs,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_genome_version,
      ref_data_lilac_resource_dir,
    )
    ch_versions = ch_versions.mix(LILAC_PROCESS.out.versions)

  emit:
    results = LILAC_PROCESS.out.lilac_dir // channel: [val(meta), lilac_dir]

    versions = ch_versions                // channel: [versions.yml]
}
