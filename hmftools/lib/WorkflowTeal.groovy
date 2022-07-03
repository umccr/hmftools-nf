//
// This file holds several functions specific to the subworkflows/teal.nf in the nf-core/hmftools pipeline
//

class WorkflowTeal {

  public static get_metrics_inputs(ch) {
    def d = ch
      .flatMap { meta, tbam, nbam, cobalt_dir, purple_dir ->
        def sample_types = ['tumor': tbam, 'normal': nbam]
        sample_types
          .keySet()
          .collect { sample_type ->
            def fp = sample_types[sample_type]
            def sample_name = meta.get(['sample_name', sample_type])
            def meta_teal = [
              id: sample_name,
              sample_type: sample_type,
              meta_full: meta,
            ]
            return [meta_teal, fp]
          }
      }
    return d
  }

  public static get_unique_input_files(ch) {
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
          metas_full: meta_fulls,
          sample_type: sample_type[0],
        ]
        return [meta_teal_new, *filepaths]
      }
    return d
  }

  public static sort_bams_and_metrics(ch) {
    def d = ch
      .flatMap{ data ->
        def meta_teal = data[0]
        def files = data[1..-1]
        meta_teal.metas_full.collect { meta -> [meta.id, meta, [meta_teal.sample_type, *files]] }
      }
      .groupTuple()
      .map { id, meta, other ->
        def data = [:]
        other.each { sample_type, bam, metrics ->
          data[[sample_type, 'bam']] = bam
          data[[sample_type, 'metrics']] = metrics
        }
        [
          meta[0],
          data.get(['tumor', 'bam']),
          data.get(['normal', 'bam']),
          data.get(['tumor', 'metrics']),
          data.get(['normal', 'metrics']),
        ]
      }
    return d
  }
}
