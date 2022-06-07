//
// This file holds several functions specific to the subworkflow/gridss.nf in the nf-core/hmftools pipeline
//

import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel

class WorkflowGridss {

  public static Map<String, DataflowWriteChannel> split_by_sample_type(channel) {
    def sample_types = ['tumor', 'normal']
    def d = channel.map { meta, sample_name, filepath ->
      def sample_type = sample_types
        .collectEntries { [meta.get(['sample_name', it]), it] }
        .get(sample_name)
      return [sample_type, meta, filepath]
    }
    return sample_types
      .collectEntries { sample_type ->
        [sample_type, d.filter { it[0] == sample_type }.map { it[1..-1] }]
      }
  }
}
