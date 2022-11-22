import Stages
import Constants


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowHmftools.initialise(params, workflow, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
  params.input,
  // Reference genome
  params.ref_data_genome_fa,
  params.ref_data_genome_fai,
  params.ref_data_genome_dict,
  params.ref_data_genome_bwa_index,
  params.ref_data_genome_bwa_index_image,
  params.ref_data_genome_gridss_index,
  // AMBER and COBALT
  params.ref_data_amber_loci,
  params.ref_data_cobalt_gc_profile,
  // CUPPA
  params.ref_data_cuppa,
  // SVPREP, GRIDSS, GRIPSS
  params.gridss_config,
  params.ref_data_sv_prep_blacklist,
  params.ref_data_gridss_blacklist,
  params.ref_data_gridss_breakend_pon,
  params.ref_data_gridss_breakpoint_pon,
  params.ref_data_repeat_masker_file,
  // Isofox
  params.ref_data_rna_exp_counts,
  params.ref_data_rna_exp_gc_ratios,
  // LINX
  params.ref_data_linx_fragile_sites,
  params.ref_data_linx_lines,
  // SAGE, PAVE
  params.ref_data_sage_blacklist_bed,
  params.ref_data_sage_blacklist_vcf,
  params.ref_data_sage_coding_panel,
  params.ref_data_sage_high_confidence,
  params.ref_data_sage_known_hotspots_germline,
  params.ref_data_sage_known_hotspots_somatic,
  params.ref_data_sage_pon_file,
  // LILAC
  params.ref_data_lilac_resource_dir,
  // VIRUSBreakend, Virus Interpreter
  params.ref_data_virusbreakenddb,
  params.ref_data_virus_taxonomy,
  params.ref_data_virus_reporting,
  // Other
  params.ref_data_purple_germline_del,
  params.ref_data_clinvar_vcf,
  params.ref_data_driver_gene_panel,
  params.ref_data_ensembl_data_dir,
  params.ref_data_known_fusion_data,
  params.ref_data_known_fusions,
  params.ref_data_mappability_bed,
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES
//
include { CHECK_SAMPLESHEET } from '../modules/local/check_samplesheet/main'

include { AMBER            } from '../modules/umccr/nextflow_modules/amber/main'
include { COBALT           } from '../modules/umccr/nextflow_modules/cobalt/main'
include { CUPPA_CLASSIFIER } from '../modules/umccr/nextflow_modules/cuppa/classifier/main'
include { CUPPA_VISUALISER } from '../modules/umccr/nextflow_modules/cuppa/visualiser/main'
include { ISOFOX           } from '../modules/umccr/nextflow_modules/isofox/main'
include { LINX_REPORT      } from '../modules/umccr/nextflow_modules/gpgr/linx_report/main'
include { PURPLE           } from '../modules/umccr/nextflow_modules/purple/main'
include { TEAL             } from '../modules/umccr/nextflow_modules/teal/main'
include { VIRUSBREAKEND    } from '../modules/umccr/nextflow_modules/virusbreakend/main'
include { VIRUSINTERPRETER } from '../modules/umccr/nextflow_modules/virusinterpreter/main'

include { PICARD_COLLECTWGSMETRICS as COLLECTWGSMETRICS } from '../modules/nf-core/modules/picard/collectwgsmetrics/main'

//
// SUBWORKFLOWS
//
include { GRIDSS        } from '../subworkflows/local/gridss'
include { GRIDSS_SVPREP } from '../subworkflows/local/gridss_svprep'
include { GRIPSS        } from '../subworkflows/local/gripss'
include { LILAC         } from '../subworkflows/local/lilac'
include { LINX          } from '../subworkflows/local/linx'
include { PAVE          } from '../subworkflows/local/pave'
include { SAGE          } from '../subworkflows/local/sage'
include { TEAL          } from '../subworkflows/local/teal'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Get absolute file paths
def get_file_object(path) {
  return path ? file(path) : []
}

samplesheet                           = get_file_object(params.input)
// Reference genome
ref_data_genome_fa                    = get_file_object(params.ref_data_genome_fa)
ref_data_genome_version               = params.ref_data_genome_version
ref_data_genome_fai                   = get_file_object(params.ref_data_genome_fai)
ref_data_genome_dict                  = get_file_object(params.ref_data_genome_dict)
ref_data_genome_bwa_index             = get_file_object(params.ref_data_genome_bwa_index)
ref_data_genome_bwa_index_image       = get_file_object(params.ref_data_genome_bwa_index_image)
ref_data_genome_gridss_index          = get_file_object(params.ref_data_genome_gridss_index)
// AMBER and COBALT
ref_data_amber_loci                   = get_file_object(params.ref_data_amber_loci)
ref_data_cobalt_gc_profile            = get_file_object(params.ref_data_cobalt_gc_profile)
// CUPPA
ref_data_cuppa                        = get_file_object(params.ref_data_cuppa)
// SVPREP, GRIDSS, GRIPSS
gridss_config                         = get_file_object(params.gridss_config)
ref_data_sv_prep_blacklist            = get_file_object(params.ref_data_sv_prep_blacklist)
ref_data_gridss_blacklist             = get_file_object(params.ref_data_gridss_blacklist)
ref_data_gridss_breakend_pon          = get_file_object(params.ref_data_gridss_breakend_pon)
ref_data_gridss_breakpoint_pon        = get_file_object(params.ref_data_gridss_breakpoint_pon)
ref_data_repeat_masker_file           = get_file_object(params.ref_data_repeat_masker_file)
// Isofox
ref_data_rna_exp_counts               = get_file_object(params.ref_data_rna_exp_counts)
ref_data_rna_exp_gc_ratios            = get_file_object(params.ref_data_rna_exp_gc_ratios)
// LINX
ref_data_linx_fragile_sites           = get_file_object(params.ref_data_linx_fragile_sites)
ref_data_linx_lines                   = get_file_object(params.ref_data_linx_lines)
// SAGE, PAVE
ref_data_sage_blacklist_bed           = get_file_object(params.ref_data_sage_blacklist_bed)
ref_data_sage_blacklist_vcf           = get_file_object(params.ref_data_sage_blacklist_vcf)
ref_data_sage_coding_panel            = get_file_object(params.ref_data_sage_coding_panel)
ref_data_sage_high_confidence         = get_file_object(params.ref_data_sage_high_confidence)
ref_data_sage_known_hotspots_germline = get_file_object(params.ref_data_sage_known_hotspots_germline)
ref_data_sage_known_hotspots_somatic  = get_file_object(params.ref_data_sage_known_hotspots_somatic)
ref_data_sage_pon_file                = get_file_object(params.ref_data_sage_pon_file)
// LILAC
ref_data_lilac_resource_dir           = get_file_object(params.ref_data_lilac_resource_dir)
// VIRUSBreakend, Virus Interpreter
ref_data_virusbreakenddb              = get_file_object(params.ref_data_virusbreakenddb)
ref_data_virus_taxonomy               = get_file_object(params.ref_data_virus_taxonomy)
ref_data_virus_reporting              = get_file_object(params.ref_data_virus_reporting)
// Other
ref_data_purple_germline_del          = get_file_object(params.ref_data_purple_germline_del)
ref_data_clinvar_vcf                  = get_file_object(params.ref_data_clinvar_vcf)
ref_data_driver_gene_panel            = get_file_object(params.ref_data_driver_gene_panel)
ref_data_ensembl_data_dir             = get_file_object(params.ref_data_ensembl_data_dir)
ref_data_known_fusion_data            = get_file_object(params.ref_data_known_fusion_data)
ref_data_known_fusions                = get_file_object(params.ref_data_known_fusions)
ref_data_mappability_bed              = get_file_object(params.ref_data_mappability_bed)

// Set stages to run
stages = Stages.set_stages(params.mode, log)

workflow HMFTOOLS {
  // Create channel for versions
  // channel: [versions.yml]
  ch_versions = Channel.empty()

  // Check samplesheet and prepare input channel
  CHECK_SAMPLESHEET(
    samplesheet,
    params.mode,
  )
  // channel: [val(meta)]
  ch_inputs = WorkflowHmftools.prepare_inputs(CHECK_SAMPLESHEET.out, workflow.stubRun, log)

  // Set up channel with common inputs for several stages
  def run_amber = Constants.Stage.AMBER in stages
  def run_cobalt = Constants.Stage.COBALT in stages
  def run_lilac = Constants.Stage.LILAC in stages
  def run_pave = Constants.Stage.PAVE in stages
  def run_teal = Constants.Stage.TEAL in stages
  if (run_amber || run_cobalt || run_pave || run_lilac || run_teal) {
    // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
    ch_bams_and_indices = ch_inputs
      .map { meta ->
        def tumor_bam = meta.get(['bam_wgs', 'tumor'])
        def normal_bam = meta.get(['bam_wgs', 'normal'])
        [meta, tumor_bam, normal_bam, "${tumor_bam}.bai", "${normal_bam}.bai"]
      }
  }

  // channel (present): [val(meta)]
  // channel (absent): [val(meta)]
  ch_inputs_wts = ch_inputs
    .branch { meta ->
      def key = ['bam_wts', 'tumor']
      present: meta.containsKey(key)
        return meta
      absent: ! meta.containsKey(key)
        return meta
    }


  //
  // MODULE: Run Isofox to analyse WTS data
  //
  ch_isofox_out = Channel.empty()
  if (Constants.Stage.ISOFOX in stages) {
    // channel: [meta, tumor_bam_wts]
    ch_isofox_inputs = ch_inputs_wts.present
      .map { meta ->
        return [meta, meta.get(['bam_wts', 'tumor'])]
      }
    ISOFOX(
      ch_isofox_inputs,
      ref_data_rna_exp_counts,
      ref_data_rna_exp_gc_ratios,
    )
    ch_versions = ch_versions.mix(ISOFOX.out.versions)
    ch_isofox_out = ch_isofox_out.mix(ISOFOX.out.isofox_dir)
  }

  //
  // MODULE: Run COLLECTWGSMETRICS to generate stats required for downstream processes
  //
  def run_cuppa = Constants.Stage.TEAL in stages
  if (run_cuppa || run_teal) {

    // NOTE(SW): CUPPA only requires collectwgsmetrics for the tumor sample in
    // the upstream process Virus Interpreter but TEAL currently requires
    // collectwgsmetrics for both tumor and normal sample
    // channel: [val(meta_cwm), bam]
    ch_cwm_inputs_all = ch_inputs
      .flatMap { meta ->
        def sample_types = run_teal ? ['tumor', 'normal'] : ['tumor']
        return sample_types
          .collect { sample_type ->
            def bam = meta.get(['bam_wgs', sample_type])
            def sample_name = meta.get(['sample_name', sample_type])
            def meta_cwm = [
              id: sample_name,
              sample_type: sample_type,
              meta_full: meta,
            ]
            return [meta_cwm, bam]
          }
      }

    // Gather duplicate files e.g. repeated normal BAMs for multiple tumor samples
    // NOTE(SW): no effective blocking by .groupTuple() as we're not dependent
    // on any process
    // channel: [val(meta_cwm), bam]
    ch_cwm_inputs = ch_cwm_inputs_all
      .map { [it[1..-1], it[0]] }
      .groupTuple()
      .map { filepaths, meta_cwm ->
        def (meta_fulls, sample_types) = meta_cwm
          .collect {
            [it.meta_full, it.sample_type]
          }
          .transpose()

        def sample_type = sample_types.unique(false)
        assert sample_type.size() == 1

        def id = meta_fulls.collect { it.id }.join('__')
        def meta_cwm_new = [
          id: "${id}_${sample_type[0]}",
          id_simple: id,
          metas_full: meta_fulls,
          sample_type: sample_type[0],
        ]
        return [meta_cwm_new, *filepaths]
      }

    COLLECTWGSMETRICS(
      ch_cwm_inputs,
      ref_data_genome_fa,
    )
    ch_versions = ch_versions.mix(COLLECTWGSMETRICS.out.versions)

    // Replicate outputs to undo unique operation
    // channel (tumor): [val(meta), metrics]
    // channel (normal): [val(meta), metrics]
    ch_cwm_output = COLLECTWGSMETRICS.out.metrics
      .flatMap { meta_cwm, metrics ->
        meta_cwm.metas_full.collect { meta -> [meta, meta_cwm.sample_type, metrics] }
      }
      .branch { meta, sample_type, metrics ->
        tumor: sample_type == 'tumor'
          return [meta, metrics]
        normal: sample_type == 'normal'
          return [meta, metrics]
      }
  }

  //
  // MODULE: Run AMBER to obtain b-allele frequencies
  //
  // channel: [val(meta), amber_dir]
  ch_amber_out = Channel.empty()
  if (run_amber) {
    AMBER(
      ch_bams_and_indices,
      ref_data_amber_loci,
    )
    ch_versions = ch_versions.mix(AMBER.out.versions)
    ch_amber_out = ch_amber_out.mix(AMBER.out.amber_dir)
  }

  //
  // MODULE: Run COBALT to obtain read ratios
  //
  // channel: [val(meta), cobalt_dir]
  ch_cobalt_out = Channel.empty()
  if (run_cobalt) {
    COBALT(
      ch_bams_and_indices,
      ref_data_cobalt_gc_profile,
    )
    ch_versions = ch_versions.mix(COBALT.out.versions)
    ch_cobalt_out = ch_cobalt_out.mix(COBALT.out.cobalt_dir)
  }

  //
  // SUBWORKFLOW: Call structural variants with GRIDSS
  //
  // channel: [val(meta), gridss_vcf]
  ch_gridss_out = Channel.empty()
  if (Constants.Stage.GRIDSS in stages) {
    if (params.no_svprep) {
      gridss_inputs = ch_inputs
        .map { meta ->
          def tumor_bam = meta.get(['bam', 'tumor'])
          def normal_bam = meta.get(['bam', 'normal'])
          [meta, tumor_bam, normal_bam]
        }
      GRIDSS(
        gridss_inputs,
        gridss_config,
        ref_data_genome_fa,
        ref_data_genome_fai,
        ref_data_genome_dict,
        ref_data_genome_bwa_index,
        ref_data_genome_bwa_index_image,
        ref_data_genome_gridss_index,
        ref_data_gridss_blacklist,
      )
      ch_versions = ch_versions.mix(GRIDSS.out.versions)
      ch_gridss_out = ch_gridss_out.mix(GRIDSS.out.results)
    } else {
      GRIDSS_SVPREP(
        ch_inputs,
        gridss_config,
        ref_data_genome_fa,
        ref_data_genome_version,
        ref_data_genome_fai,
        ref_data_genome_dict,
        ref_data_genome_bwa_index,
        ref_data_genome_bwa_index_image,
        ref_data_genome_gridss_index,
        ref_data_gridss_blacklist,
        ref_data_sv_prep_blacklist,
        ref_data_known_fusions,
      )
      ch_versions = ch_versions.mix(GRIDSS_SVPREP.out.versions)
      ch_gridss_out = ch_gridss_out.mix(GRIDSS_SVPREP.out.results)
    }
  }

  //
  // MODULE: Run GRIPSS to filter GRIDSS SV calls
  //
  // channel: [val(meta), hard_vcf, hard_tbi, soft_vcf, soft_tbi]
  ch_gripss_germline_out = Channel.empty()
  // channel: [val(meta), hard_vcf, hard_tbi, soft_vcf, soft_tbi]
  ch_gripss_somatic_out = Channel.empty()
  if (Constants.Stage.GRIPSS in stages) {
    GRIPSS(
      ch_gridss_out,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_version,
      ref_data_gridss_breakend_pon,
      ref_data_gridss_breakpoint_pon,
      ref_data_known_fusions,
      ref_data_repeat_masker_file,
    )
    ch_versions = ch_versions.mix(GRIPSS.out.versions)
    ch_gripss_germline_out = ch_gripss_germline_out.mix(GRIPSS.out.germline)
    ch_gripss_somatic_out = ch_gripss_somatic_out.mix(GRIPSS.out.somatic)
  }

  //
  // SUBWORKFLOW: call SNV, MNV, and small INDELS with SAGE
  //
  // channel: [val(meta), sage_vcf]
  ch_sage_germline_out = Channel.empty()
  // channel: [val(meta), sage_vcf]
  ch_sage_somatic_out = Channel.empty()
  if (Constants.Stage.SAGE in stages) {
    SAGE(
      ch_bams_and_indices,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_dict,
      ref_data_genome_version,
      ref_data_sage_known_hotspots_germline,
      ref_data_sage_known_hotspots_somatic,
      ref_data_sage_coding_panel,
      ref_data_sage_high_confidence,
      ref_data_sage_pon_file,
      ref_data_mappability_bed,
      ref_data_driver_gene_panel,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(SAGE.out.versions)
    ch_sage_germline_out = ch_sage_germline_out.mix(SAGE.out.germline)
    ch_sage_somatic_out = ch_sage_somatic_out.mix(SAGE.out.somatic)
  }

  //
  // SUBWORKFLOW: Annotate variants with PAVE
  //
  // channel: [val(meta), pave_vcf]
  ch_pave_germline_out = Channel.empty()
  // channel: [val(meta), pave_vcf]
  ch_pave_somatic_out = Channel.empty()
  if (run_pave) {
    PAVE(
      ch_sage_germline_out,
      ch_sage_somatic_out,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_version,
      ref_data_sage_pon_file,
      ref_data_sage_blacklist_bed,
      ref_data_sage_blacklist_vcf,
      ref_data_clinvar_vcf,
      ref_data_mappability_bed,
      ref_data_driver_gene_panel,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(PAVE.out.versions)
    ch_pave_germline_out = ch_pave_germline_out.mix(PAVE.out.germline)
    ch_pave_somatic_out = ch_pave_somatic_out.mix(PAVE.out.somatic)
  }

  //
  // MODULE: Run PURPLE for CNV calling, purity and ploidy inference, SV recovery
  //
  // channel: [val(meta), purple_dir]
  ch_purple_out = Channel.empty()
  if (Constants.Stage.PURPLE in stages) {

    // Mode: full
    if (Constants.Stage.PAVE in stages) {
      ch_purple_inputs = WorkflowHmftools.group_by_meta(
        ch_amber_out,
        ch_cobalt_out,
        ch_gripss_somatic_out,
        ch_pave_somatic_out,
        ch_pave_germline_out,
      )
    // Mode: gridss-purple-linx
    } else if (Constants.Stage.AMBER in stages && Constants.Stage.COBALT in stages) {
      ch_purple_inputs = WorkflowHmftools.group_by_meta(
        ch_amber_out,
        ch_cobalt_out,
        ch_gripss_somatic_out,
      )
        .map { data ->
          def meta = data[0]
          def fps = data[1..-1]
          return [
            meta,
            *fps,
            meta.get(['vcf_smlv', 'tumor'], []),
            meta.get(['vcf_smlv', 'normal'], []),
          ]
        }
    // Mode: purple
    } else {
      ch_purple_inputs = ch_inputs
        .map { meta ->
          def sv_hard_vcf = meta[['vcf_sv_gripss_hard', 'tumor']]
          def sv_soft_vcf = meta[['vcf_sv_gripss_soft', 'tumor']]
          return [
            meta,
            meta[['amber_dir', 'tumor']],
            meta[['cobalt_dir', 'tumor']],
            sv_hard_vcf,
            "${sv_hard_vcf}.tbi",
            sv_soft_vcf,
            "${sv_soft_vcf}.tbi",
            meta.get(['vcf_smlv', 'tumor'], []),
            meta.get(['vcf_smlv', 'normal'], []),
          ]
        }
    }

    PURPLE(
      ch_purple_inputs,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_dict,
      ref_data_genome_version,
      ref_data_cobalt_gc_profile,
      ref_data_sage_known_hotspots_somatic,
      ref_data_sage_known_hotspots_germline,
      ref_data_driver_gene_panel,
      ref_data_ensembl_data_dir,
      ref_data_purple_germline_del,
    )
    ch_versions = ch_versions.mix(PURPLE.out.versions)
    ch_purple_out = ch_purple_out.mix(PURPLE.out.purple_dir)
  }

  //
  // MODULE: Run TEAL to characterise teleomeres
  //
  if (run_teal) {

    // Mode: full
    // channel: [val(meta), cobalt_dir, purple_dir]
    if (run_cobalt && Constants.Stage.PURPLE in stages) {
      ch_teal_inputs_other = WorkflowHmftools.group_by_meta(
        ch_cobalt_out,
        ch_purple_out,
      )
    // Mode: teal
    } else {
      ch_teal_inputs_other = ch_inputs
        .map { meta ->
          return [
            meta,
            meta.get(['cobalt_dir', 'tumor']),
            meta.get(['purple_dir', 'tumor']),
          ]
        }
    }

    // channel: [val(meta), tumor_wgs_metrics, normal_wgs_metrics]
    // NOTE(SW): assuming here that TEAL is being run in tumor/normal mode and
    // so we expect a tumor metrics file and normal metrics file
    ch_teal_inputs_metrics = WorkflowHmftools.group_by_meta(
      ch_cwm_output.tumor,
      ch_cwm_output.normal,
    )

    // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai, tumor_wgs_metrics, normal_wgs_metrics, cobalt_dir, purple_dir]
    ch_teal_inputs = WorkflowHmftools.group_by_meta(
      ch_bams_and_indices,
      ch_teal_inputs_metrics,
      ch_teal_inputs_other,
    )

    TEAL(
      ch_teal_inputs,
    )
    ch_versions = ch_versions.mix(TEAL.out.versions)
  }

  //
  // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
  //
  if (run_lilac) {

    // Mode: full
    if (Constants.Stage.PURPLE in stages) {
      ch_lilac_inputs_purple_dir = ch_purple_out
    // Mode: lilac
    } else {
      ch_lilac_inputs_purple_dir = ch_inputs
        .map { meta ->
          return [meta, meta.get(['purple_dir', 'tumor'])]
        }
    }

    LILAC(
      ch_bams_and_indices,
      ch_lilac_inputs_purple_dir,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_version,
      ref_data_lilac_resource_dir,
    )
    ch_versions = ch_versions.mix(LILAC.out.versions)
  }

  //
  // MODULE: Run VIRUSBreakend and Virus Interpreter to quantify viral content
  //
  // NOTE(SW): kept separate from CUPPA conditional block since we'll allow users to run this independently
  ch_virusinterpreter_out = Channel.empty()
  if (Constants.Stage.CUPPA in stages) {

    // TODO(SW): PURPLE from difference sources

    // channel: [val(meta), tumor_bam]
    ch_virusbreakend_inputs = ch_inputs.map { meta -> [meta, meta.get(['bam_wgs', 'tumor'])] }

    VIRUSBREAKEND(
      ch_virusbreakend_inputs,
      ref_data_virusbreakenddb,
      gridss_config,
      ref_data_genome_fa,
      ref_data_genome_fai,
      ref_data_genome_dict,
      ref_data_genome_bwa_index,
      ref_data_genome_bwa_index_image,
      ref_data_genome_gridss_index,
    )
    ch_versions = ch_versions.mix(VIRUSBREAKEND.out.versions)

    // channel: [val(meta), purple_purity, purple_qc]
    ch_virusinterpreter_inputs_purple = PURPLE.out.purple_dir
      .map { meta, purple_dir ->
        def purple_purity = file(purple_dir).resolve("${meta.get(['sample_name', 'tumor'])}.purple.purity.tsv")
        def purple_qc = file(purple_dir).resolve("${meta.get(['sample_name', 'tumor'])}.purple.qc")
        return [meta, purple_purity, purple_qc]
      }

    // channel: [val(meta), virus_tsv, purple_purity, purple_qc, wgs_metrics]
    ch_virusinterpreter_inputs = WorkflowHmftools.group_by_meta(
      VIRUSBREAKEND.out.tsv,
      ch_virusinterpreter_inputs_purple,
      ch_cwm_output.tumor,
    )

    VIRUSINTERPRETER(
      ch_virusinterpreter_inputs,
      ref_data_virus_taxonomy,
      ref_data_virus_reporting,
    )
    ch_versions = ch_versions.mix(VIRUSINTERPRETER.out.versions)
    ch_virusinterpreter_out = ch_virusinterpreter_out.mix(VIRUSINTERPRETER.out.virusinterpreter_dir)
  }

  //
  // SUBWORKFLOW: Group structural variants into higher order events with LINX
  //
  // channel: [val(meta), linx_annotation_dir, linx_visuliaser_dir]
  ch_linx_somatic_out = Channel.empty()
  if (Constants.Stage.LINX in stages) {

    // Mode: full
    if (Constants.Stage.PURPLE in stages && Constants.Stage.GRIPSS in stages) {
      ch_linx_germline_inputs = ch_gripss_germline_out.map { meta, h, hi, s, si -> [meta, h] }
      ch_linx_somatic_inputs = ch_purple_out
    // Mode: linx
    } else {
      ch_linx_germline_inputs = ch_inputs
        .map { meta ->
          def sv = meta.get(['vcf_sv_gripss_hard', 'normal'], false)
          return sv ? [meta, sv] : false
        }
        .filter { it }
      ch_linx_somatic_inputs = ch_inputs.map { meta -> [meta, meta.get(['purple_dir', 'tumor'], []) ] }
    }

    LINX(
      ch_linx_germline_inputs,
      ch_linx_somatic_inputs,
      ref_data_genome_version,
      ref_data_linx_fragile_sites,
      ref_data_linx_lines,
      ref_data_ensembl_data_dir,
      ref_data_known_fusion_data,
      ref_data_driver_gene_panel,
    )
    ch_versions = ch_versions.mix(LINX.out.versions)
    ch_linx_somatic_out = ch_linx_somatic_out.mix(LINX.out.somatic)

    //
    // MODULE: Run gpgr to generate a LINX report
    //
    LINX_REPORT(
      ch_linx_somatic_out,
    )
    ch_versions = ch_versions.mix(LINX_REPORT.out.versions)
  }

  //
  // MODULE: Run CUPPA predict tissue of origin
  //
  if (Constants.Stage.CUPPA in stages) {

    // channel: [val(meta), isofox_dir]
    ch_cuppa_inputs_wts = Channel.empty()
      .mix(
        ch_inputs_wts.absent.map { meta -> [meta, []] },
        ch_isofox_out,
      )

    // channel: [val(meta), isofox_dir, purple_dir, linx_dir, virusinterpreter_dir]
    // NOTE(SW): the Groovy Collection.flatten method used in
    // WorkflowHmftools.group_by_meta removes optional Isofox input; flattening
    // done manually below to preserve
    ch_cuppa_inputs = WorkflowHmftools.group_by_meta(
      ch_cuppa_inputs_wts,
      ch_purple_out,
      ch_linx_somatic_out.map { meta, anno_dir, vis_dir -> [meta, anno_dir]},
      ch_virusinterpreter_out,
      flatten: false,
    )
      .map { data ->
        def meta = data[0]
        def inputs = data[1..-1].collect { it[0] }
        return [meta, *inputs]
      }

    CUPPA_CLASSIFIER(
      ch_cuppa_inputs,
      ref_data_cuppa,
    )
    ch_versions = ch_versions.mix(CUPPA_CLASSIFIER.out.versions)

    CUPPA_VISUALISER(
      CUPPA_CLASSIFIER.out.csv,
    )
    ch_versions = ch_versions.mix(CUPPA_VISUALISER.out.versions)
  }

  //
  // MODULE: Pipeline reporting
  //
  CUSTOM_DUMPSOFTWAREVERSIONS(
    ch_versions.unique().collectFile(name: 'collated_versions.yml')
  )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
  if (params.email || params.email_on_fail) {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
  }
  NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
