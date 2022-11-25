//
// STAR is short RNA read aligner.
//

include { PICARD_MARKDUPLICATES as MARKDUPLICATES } from '../../modules/local/picard/markduplicates/main'

include { STAR_ALIGN as ALIGN                     } from '../../modules/nf-core/star/align/main'
include { SAMTOOLS_SORT as SORT                   } from '../../modules/nf-core/samtools/sort/main'

workflow STAR {
  take:
    ch_inputs                             // channel: [meta, [fastq_wts_fwd, fastq_wts_rev]]
    ref_data_genome_star_index            //    file: /path/to/genome_star_index
    ref_data_genome_genes_gtf             //    file: /path/to/genome_genes_gtf

  main:
    // Channel for version.yml files
    ch_versions = Channel.empty()

    ALIGN(
      ch_inputs,
      ref_data_genome_star_index,
      ref_data_genome_genes_gtf,
      false, // star_ignore_sjdbgtf
      '',    // seq_platform
      '',    // seq_center
    )
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    SORT(
      ALIGN.out.bam,
    )
    ch_versions = ch_versions.mix(SORT.out.versions)

    MARKDUPLICATES(
      SORT.out.bam,
    )
    ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions)

  emit:
    bam      = MARKDUPLICATES.out.bam // channel: [val(meta), bam]

    versions = ch_versions            // channel: [versions.yml]
}

