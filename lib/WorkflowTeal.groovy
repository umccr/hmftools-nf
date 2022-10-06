//
// This file holds several functions specific to the subworkflows/teal.nf in the umccr/hmftools pipeline
//

class WorkflowTeal {

  public static split_input_bams(ch) {
    // channel: [val(meta_teal), bam]
    def d = ch
      .flatMap { meta, tbam, nbam, tbai, nbai ->
        def sample_types = ['tumor': [tbam, tbai], 'normal': [nbam, nbai]]
        sample_types
          .keySet()
          .collect { sample_type ->
            def fps = sample_types[sample_type]
            def sample_name = meta.get(['sample_name', sample_type])
            def meta_teal = [
              id: sample_name,
              sample_type: sample_type,
              meta_full: meta,
            ]
            return [meta_teal, *fps]
          }
      }
    return d
  }

  public static get_unique_input_files(ch) {
    // channel: [val(meta_teal), bam, bai]
    def d = ch
      .map { [it[1..-1], it[0]] }
      .groupTuple()
      .map { filepaths, meta_teals ->
        def (meta_fulls, sample_types) = meta_teals
          .collect {
            [it.meta_full, it.sample_type]
          }
          .transpose()

        def sample_type = sample_types.unique(false)
        assert sample_type.size() == 1

        def id = meta_fulls.collect { it.id }.join('__')
        def meta_teal_new = [
          id: "${id}_${sample_type[0]}",
          id_simple: id,
          metas_full: meta_fulls,
          sample_type: sample_type[0],
        ]
        return [meta_teal_new, *filepaths]
      }
    return d
  }

  public static sort_bams_and_metrics(ch) {
    // Collect T/N pairs into single channel element
    // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, tumor_wgs_metrics, normal_wgs_metrics]
    def d = ch
      .flatMap{ data ->
        def meta_teal = data[0]
        def files = data[1..-1]
        meta_teal.metas_full.collect { meta -> [meta.id, meta, [meta_teal.sample_type, *files]] }
      }
      .groupTuple(size: 2)
      .map { id, meta, other ->
        def data = [:]
        other.each { sample_type, bam, bai, metrics ->
          data[[sample_type, 'bam']] = bam
          data[[sample_type, 'bai']] = bai
          data[[sample_type, 'metrics']] = metrics
        }
        [
          meta[0],
          data.get(['tumor', 'bam']),
          data.get(['normal', 'bam']),
          data.get(['tumor', 'bai']),
          data.get(['normal', 'bai']),
          data.get(['tumor', 'metrics']),
          data.get(['normal', 'metrics']),
        ]
      }
    return d
  }
}
