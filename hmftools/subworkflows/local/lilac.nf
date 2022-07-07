//
// LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
//

include { SAMBAMBA_SLICE } from '../../modules/scwatts/nextflow_modules/sambamba/slice/main'

include { LILAC as LILAC_PROCESS } from '../../modules/scwatts/nextflow_modules/lilac/main'

workflow LILAC {
  take:
    ch_inputs                   // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, purple_dir]
    ref_data_genome_dir         //    file: /path/to/genome_dir/
    ref_data_genome_fn          //     val: genome name
    ref_data_lilac_resource_dir //    file: /path/to/lilac_resource_dir/

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    // Slice HLA region
    // NOTE(SW): here I remove duplicate files so that we only process each input once
    // channel: [val(meta_lilac), bam, bai, bed, regions_list]
    ch_slice_inputs = WorkflowLilac.get_slice_inputs(ch_inputs, "${ref_data_lilac_resource_dir}/hla.38.bed")
    // channel: [val(meta_lilac), bam, bai, bed, regions_list]
    ch_slice_inputs = WorkflowLilac.get_unique_input_files(ch_slice_inputs)
    SAMBAMBA_SLICE(
      ch_slice_inputs,
    )
    ch_versions = ch_versions.mix(SAMBAMBA_SLICE.out.versions)

    // Run LILAC
    // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, purple_dir]
    ch_lilac_inputs = WorkflowHmftools.group_by_meta(
      WorkflowLilac.sort_slices(SAMBAMBA_SLICE.out.bam),
      ch_inputs.map { meta, tbam, nbam, tbai, nbai, purple_dir -> [meta, purple_dir] }
    )
    LILAC_PROCESS(
      ch_lilac_inputs,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_lilac_resource_dir,
    )
    ch_versions = ch_versions.mix(LILAC_PROCESS.out.versions)

  emit:
    results = LILAC_PROCESS.out.lilac_dir // channel: [val(meta), lilac_dir]

    versions = ch_versions                // channel: [versions.yml]
}
