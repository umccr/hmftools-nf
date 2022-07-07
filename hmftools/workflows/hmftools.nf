/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowHmftools.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
  params.input,
  params.gridss_config,

  params.ref_data_genome,
  params.ref_data_amber_loci,
  params.ref_data_cobalt_gc_profile,
  params.ref_data_gridss_blacklist,
  params.ref_data_gridss_breakend_pon,
  params.ref_data_gridss_breakpoint_pon,
  params.ref_data_linx_fragile_sites,
  params.ref_data_linx_line_elements,
  params.ref_data_ensembl_data_dir,

  params.ref_data_sage_known_hotspots_somatic,
  params.ref_data_sage_known_hotspots_germline,

  params.ref_data_sage_coding_panel_germline,
  params.ref_data_sage_coding_panel_somatic,
  params.ref_data_sage_high_confidence,
  params.ref_data_sage_pon_file,
  params.ref_data_mappability_bed,

  params.ref_data_sage_blacklist_bed,
  params.ref_data_sage_blacklist_vcf,
  params.ref_data_clinvar_vcf,

  params.ref_data_lilac_resource_dir,

  params.ref_data_known_fusions,
  params.ref_data_known_fusion_data,
  params.ref_data_driver_gene_panel,
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
// MODULE: Loaded from modules/local/
//
include { CHECK_SAMPLESHEET } from '../modules/local/check_samplesheet'

//
// MODULE: Loaded from modules/scwatts/
//
include { AMBER       } from '../modules/scwatts/nextflow_modules/amber/main'
include { COBALT      } from '../modules/scwatts/nextflow_modules/cobalt/main'
include { LINX_REPORT } from '../modules/scwatts/nextflow_modules/gpgr/linx_report/main'
include { PURPLE      } from '../modules/scwatts/nextflow_modules/purple/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { GRIDSS } from '../subworkflows/local/gridss'
include { GRIPSS } from '../subworkflows/local/gripss'
include { LILAC  } from '../subworkflows/local/lilac'
include { LINX   } from '../subworkflows/local/linx'
include { PAVE   } from '../subworkflows/local/pave'
include { SAGE   } from '../subworkflows/local/sage'
include { TEAL   } from '../subworkflows/local/teal'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ref_data_genome = new File(params.ref_data_genome).absolutePath
ref_data_genome_dir = file(ref_data_genome).parent
ref_data_genome_fn = file(ref_data_genome).name

ref_data_amber_loci            = new File(params.ref_data_amber_loci).absolutePath
ref_data_cobalt_gc_profile     = new File(params.ref_data_cobalt_gc_profile).absolutePath
ref_data_gridss_blacklist      = new File(params.ref_data_gridss_blacklist).absolutePath
ref_data_gridss_breakend_pon   = new File(params.ref_data_gridss_breakend_pon).absolutePath
ref_data_gridss_breakpoint_pon = new File(params.ref_data_gridss_breakpoint_pon).absolutePath
ref_data_linx_fragile_sites    = new File(params.ref_data_linx_fragile_sites).absolutePath
ref_data_linx_line_elements    = new File(params.ref_data_linx_line_elements).absolutePath
ref_data_ensembl_data_dir      = new File(params.ref_data_ensembl_data_dir).absolutePath

ref_data_sage_known_hotspots_somatic = new File(params.ref_data_sage_known_hotspots_somatic).absolutePath
ref_data_sage_known_hotspots_germline = new File(params.ref_data_sage_known_hotspots_germline).absolutePath

ref_data_sage_coding_panel_germline = new File(params.ref_data_sage_coding_panel_germline).absolutePath
ref_data_sage_coding_panel_somatic = new File(params.ref_data_sage_coding_panel_somatic).absolutePath
ref_data_sage_high_confidence = new File(params.ref_data_sage_high_confidence).absolutePath
ref_data_sage_pon_file = new File(params.ref_data_sage_pon_file).absolutePath
ref_data_mappability_bed = new File(params.ref_data_mappability_bed).absolutePath

ref_data_sage_blacklist_bed = new File(params.ref_data_sage_blacklist_bed).absolutePath
ref_data_sage_blacklist_vcf = new File(params.ref_data_sage_blacklist_vcf).absolutePath
ref_data_clinvar_vcf = new File(params.ref_data_clinvar_vcf).absolutePath

ref_data_lilac_resource_dir = new File(params.ref_data_lilac_resource_dir).absolutePath

ref_data_known_fusions         = new File(params.ref_data_known_fusions).absolutePath
ref_data_known_fusion_data     = new File(params.ref_data_known_fusion_data).absolutePath
ref_data_driver_gene_panel     = new File(params.ref_data_driver_gene_panel).absolutePath

gridss_config         = params.gridss_config ? new File(params.gridss_config).absolutePath : []


// NOTE(SW): move to initialise
stages =  WorkflowHmftools.set_stages(params.mode, log)


def check_labels(d, log) {
  def labels_seen = []
  def labels_duplicated = []
  d.collect { label, fp ->
    labels_seen.contains(label) ? labels_duplicated << label : labels_seen << label
  }
  if (labels_duplicated) {
      log.error "\nERROR: got duplicate label(s) ${labels_duplicated.join(', ')}"
      System.exit(1)
  }
}


workflow HMFTOOLS {
  // Create channel for versions
  ch_versions = Channel.empty()

  // Check samplesheet and prepare input channel
  CHECK_SAMPLESHEET(
    params.input,
    params.mode,
  )
  ch_inputs = WorkflowHmftools.prepare_inputs(CHECK_SAMPLESHEET.out, workflow.stubRun, log)

  // Set up channel for AMBER and COBALT
  def run_amber = WorkflowHmftools.Stage.AMBER in stages
  def run_cobalt = WorkflowHmftools.Stage.COBALT in stages
  def run_pave = WorkflowHmftools.Stage.PAVE in stages
  def run_lilac = WorkflowHmftools.Stage.LILAC in stages
  if (run_amber || run_cobalt || run_pave || run_lilac) {
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
  ch_gridss_out = Channel.empty()
  if (WorkflowHmftools.Stage.GRIDSS in stages) {
    GRIDSS(
      ch_inputs,
      gridss_config,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_gridss_blacklist,
    )
    ch_versions = ch_versions.mix(GRIDSS.out.versions)
    ch_gridss_out = ch_gridss_out.mix(GRIDSS.out.results)
  }

  //
  // MODULE: Run GRIPSS to filter GRIDSS SV calls
  //
  ch_gripss_germline_out = Channel.empty()
  ch_gripss_somatic_out = Channel.empty()
  if (WorkflowHmftools.Stage.GRIPSS in stages) {
    GRIPSS(
      ch_gridss_out,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_gridss_breakend_pon,
      ref_data_gridss_breakpoint_pon,
      ref_data_known_fusions,
    )
    ch_versions = ch_versions.mix(GRIPSS.out.versions)
    ch_gripss_germline_out = ch_gripss_germline_out.mix(GRIPSS.out.germline)
    ch_gripss_somatic_out = ch_gripss_somatic_out.mix(GRIPSS.out.somatic)
  }

  //
  // SUBWORKFLOW: call SNV, MNV, and small INDELS with SAGE
  //
  ch_sage_germline_out = Channel.empty()
  ch_sage_somatic_out = Channel.empty()
  if (WorkflowHmftools.Stage.SAGE in stages) {
    SAGE(
      ch_bams_and_indices,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_sage_known_hotspots_germline,
      ref_data_sage_known_hotspots_somatic,
      ref_data_sage_coding_panel_germline,
      ref_data_sage_coding_panel_somatic,
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
  ch_pave_germline_out = Channel.empty()
  ch_pave_somatic_out = Channel.empty()
  if (run_pave) {
    PAVE(
      ch_sage_germline_out,
      ch_sage_somatic_out,
      ref_data_genome_dir,
      ref_data_genome_fn,
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
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_cobalt_gc_profile,
      ref_data_sage_known_hotspots_somatic,
      ref_data_sage_known_hotspots_germline,
      ref_data_driver_gene_panel,
      ref_data_ensembl_data_dir,
    )
    ch_versions = ch_versions.mix(PURPLE.out.versions)
    ch_purple_out = ch_purple_out.mix(PURPLE.out.purple_dir)
  }

  //
  // SUBWORKFLOW: Characterise teleomeres with TEAL
  //
  if (WorkflowHmftools.Stage.TEAL in stages) {

    ch_teal_inputs_base = ch_inputs
      .map { meta ->
        [meta, meta.get(['bam', 'tumor']), meta.get(['bam', 'normal'])]
      }

    // Mode: full
    if (run_cobalt && WorkflowHmftools.Stage.PURPLE in stages) {
      ch_teal_inputs = WorkflowHmftools.group_by_meta(
        ch_teal_inputs_base,
        ch_cobalt_out,
        ch_purple_out,
      )
    // Mode: teal
    } else {
      ch_teal_inputs = ch_teal_inputs_base
        .map { data ->
          def meta = data[0]
          def fps = data[1..-1]
          return [
            meta,
            *fps,
            meta.get(['cobalt_dir', 'tumor']),
            meta.get(['purple_dir', 'tumor']),
          ]
        }
    }

    TEAL(
      ch_teal_inputs,
      ref_data_genome_dir,
      ref_data_genome_fn,
    )
  }

  //
  // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
  //
  if (run_lilac) {

    // Mode: full
    if (WorkflowHmftools.Stage.PURPLE in stages) {
      ch_lilac_inputs = WorkflowHmftools.group_by_meta(
        ch_bams_and_indices,
        ch_purple_out,
      )
    // Mode: lilac
    } else {
      ch_lilac_inputs = ch_bams_and_indices
        .map { data ->
          def meta = data[0]
          def fps = data[1..-1]
          return [
            meta,
            *fps,
            meta.get(['purple_dir', 'tumor']),
          ]
        }
    }

    LILAC(
      ch_lilac_inputs,
      ref_data_genome_dir,
      ref_data_genome_fn,
      ref_data_lilac_resource_dir,
    )
    ch_versions = ch_versions.mix(LILAC.out.versions)
  }

  //
  // SUBWORKFLOW: Group structural variants into higher order events with LINX
  //
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
      ref_data_linx_fragile_sites,
      ref_data_linx_line_elements,
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
