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
  // SVPREP, GRIDSS, GRIPSS
  params.ref_data_sv_prep_blacklist,
  params.gridss_config,
  params.ref_data_gridss_blacklist,
  params.ref_data_gridss_breakend_pon,
  params.ref_data_gridss_breakpoint_pon,
  params.ref_data_repeat_masker_file,
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

include { AMBER       } from '../modules/umccr/nextflow_modules/amber/main'
include { COBALT      } from '../modules/umccr/nextflow_modules/cobalt/main'
include { LINX_REPORT } from '../modules/umccr/nextflow_modules/gpgr/linx_report/main'
include { PURPLE      } from '../modules/umccr/nextflow_modules/purple/main'

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
// SVPREP, GRIDSS, GRIPSS
ref_data_sv_prep_blacklist            = get_file_object(params.ref_data_sv_prep_blacklist)
gridss_config                         = get_file_object(params.gridss_config)
ref_data_gridss_blacklist             = get_file_object(params.ref_data_gridss_blacklist)
ref_data_gridss_breakend_pon          = get_file_object(params.ref_data_gridss_breakend_pon)
ref_data_gridss_breakpoint_pon        = get_file_object(params.ref_data_gridss_breakpoint_pon)
ref_data_repeat_masker_file           = get_file_object(params.ref_data_repeat_masker_file)
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
// Other
ref_data_purple_germline_del          = get_file_object(params.ref_data_purple_germline_del)
ref_data_clinvar_vcf                  = get_file_object(params.ref_data_clinvar_vcf)
ref_data_driver_gene_panel            = get_file_object(params.ref_data_driver_gene_panel)
ref_data_ensembl_data_dir             = get_file_object(params.ref_data_ensembl_data_dir)
ref_data_known_fusion_data            = get_file_object(params.ref_data_known_fusion_data)
ref_data_known_fusions                = get_file_object(params.ref_data_known_fusions)
ref_data_mappability_bed              = get_file_object(params.ref_data_mappability_bed)

// Set stages to run
stages = WorkflowHmftools.set_stages(params.mode, log)

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
  def run_amber = WorkflowHmftools.Stage.AMBER in stages
  def run_cobalt = WorkflowHmftools.Stage.COBALT in stages
  def run_pave = WorkflowHmftools.Stage.PAVE in stages
  def run_lilac = WorkflowHmftools.Stage.LILAC in stages
  if (run_amber || run_cobalt || run_pave || run_lilac) {
    // channel: [val(meta), tumor_bam, normal_bam, tumor_bai, normal_bai]
    ch_bams_and_indices = ch_inputs
      .map { meta ->
        def tumor_bam = meta.get(['bam', 'tumor'])
        def normal_bam = meta.get(['bam', 'normal'])
        [meta, tumor_bam, normal_bam, "${tumor_bam}.bai", "${normal_bam}.bai"]
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
  if (WorkflowHmftools.Stage.GRIDSS in stages) {
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
  if (WorkflowHmftools.Stage.GRIPSS in stages) {
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
  if (WorkflowHmftools.Stage.SAGE in stages) {
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
  if (WorkflowHmftools.Stage.PURPLE in stages) {

    // Mode: full
    if (WorkflowHmftools.Stage.PAVE in stages) {
      ch_purple_inputs = WorkflowHmftools.group_by_meta(
        ch_amber_out,
        ch_cobalt_out,
        ch_gripss_somatic_out,
        ch_pave_somatic_out,
        ch_pave_germline_out,
      )
    // Mode: gridss-purple-linx
    } else if (WorkflowHmftools.Stage.AMBER in stages && WorkflowHmftools.Stage.COBALT in stages) {
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
            meta.get(['smlv', 'tumor'], []),
            meta.get(['smlv', 'normal'], []),
          ]
        }
    // Mode: purple
    } else {
      ch_purple_inputs = ch_inputs
        .map { meta ->
          def sv_hard_vcf = meta[['gripss_hard_sv', 'tumor']]
          def sv_soft_vcf = meta[['gripss_soft_sv', 'tumor']]
          return [
            meta,
            meta[['amber_dir', 'tumor']],
            meta[['cobalt_dir', 'tumor']],
            sv_hard_vcf,
            "${sv_hard_vcf}.tbi",
            sv_soft_vcf,
            "${sv_soft_vcf}.tbi",
            meta.get(['smlv', 'tumor'], []),
            meta.get(['smlv', 'normal'], []),
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
  // SUBWORKFLOW: Characterise teleomeres with TEAL
  //
  if (WorkflowHmftools.Stage.TEAL in stages) {

    ch_teal_inputs_bams = ch_inputs
      .map { meta ->
        def tumor_bam = meta.get(['bam', 'tumor'])
        def normal_bam = meta.get(['bam', 'normal'])
        [meta, tumor_bam, normal_bam, "${tumor_bam}.bai", "${normal_bam}.bai"]
      }

    // Mode: full
    if (run_cobalt && WorkflowHmftools.Stage.PURPLE in stages) {
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

    TEAL(
      ch_teal_inputs_bams,
      ch_teal_inputs_other,
      ref_data_genome_fa,
    )
    ch_versions = ch_versions.mix(TEAL.out.versions)
  }

  //
  // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
  //
  if (run_lilac) {

    // Mode: full
    if (WorkflowHmftools.Stage.PURPLE in stages) {
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
  // SUBWORKFLOW: Group structural variants into higher order events with LINX
  //
  // channel: [val(meta), linx_annotation_dir, linx_visuliaser_dir]
  ch_linx_somatic_out = Channel.empty()
  if (WorkflowHmftools.Stage.LINX in stages) {

    // Mode: full
    if (WorkflowHmftools.Stage.PURPLE in stages && WorkflowHmftools.Stage.GRIPSS in stages) {
      ch_linx_germline_inputs = ch_gripss_germline_out.map { meta, h, hi, s, si -> [meta, h] }
      ch_linx_somatic_inputs = ch_purple_out
    // Mode: linx
    } else {
      ch_linx_germline_inputs = ch_inputs
        .map { meta ->
          def sv = meta.get(['gripss_hard_sv', 'normal'], false)
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
