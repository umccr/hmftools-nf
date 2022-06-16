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
  params.ref_data_genome,
  params.ref_data_amber_loci,
  params.ref_data_cobalt_gc_profile,
  params.ref_data_gridss_blacklist,
  params.ref_data_gridss_breakend_pon,
  params.ref_data_gridss_breakpoint_pon,
  params.ref_data_linx_fragile_sites,
  params.ref_data_linx_line_elements,
  params.ref_data_ensembl_data_dir,
  params.ref_data_known_hotspots,
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
include { AMBER       } from '../modules/scwatts/nextflow_modules/amber/main'
include { COBALT      } from '../modules/scwatts/nextflow_modules/cobalt/main'
include { GRIPSS      } from '../modules/scwatts/nextflow_modules/gripss/main'
include { LINX_REPORT } from '../modules/scwatts/nextflow_modules/gpgr/linx_report/main'
include { PURPLE      } from '../modules/scwatts/nextflow_modules/purple/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { GRIDSS      } from '../subworkflows/local/gridss'
include { LINX        } from '../subworkflows/local/linx'

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
ref_data_known_hotspots        = new File(params.ref_data_known_hotspots).absolutePath
ref_data_known_fusions         = new File(params.ref_data_known_fusions).absolutePath
ref_data_known_fusion_data     = new File(params.ref_data_known_fusion_data).absolutePath
ref_data_driver_gene_panel     = new File(params.ref_data_driver_gene_panel).absolutePath

workflow HMFTOOLS {
  // Set input objects
  ch_inputs = WorkflowHmftools.prepare_inputs(params.input, workflow.stubRun, log)

  // Create channel for modules and subworkflows
  ch_versions = Channel.empty()

  // Set up channel for AMBER and COBALT
  ch_bams_and_indices = ch_inputs
    .map { meta ->
      def tumor_bam = meta.get(['tumor', 'bam'])
      def normal_bam = meta.get(['normal', 'bam'])
      [meta, tumor_bam, normal_bam, "${tumor_bam}.bai", "${normal_bam}.bai"]
    }

  //
  // MODULE: Run AMBER to obtain b-allele frequencies
  //
  AMBER(
    ch_bams_and_indices,
    ref_data_amber_loci,
  )
  ch_versions = ch_versions.mix(AMBER.out.versions)

  //
  // MODULE: Run COBALT to obtain read ratios
  //
  COBALT(
    ch_bams_and_indices,
    ref_data_cobalt_gc_profile,
  )
  ch_versions = ch_versions.mix(COBALT.out.versions)

  //
  // SUBWORKFLOW: Call structural variants with GRIDSS
  //
  ch_gridss_reads_input = WorkflowHmftools.get_gridss_reads_input(ch_inputs)
  ch_gridss_svs_input = WorkflowHmftools.get_gridss_svs_input(ch_inputs)
  GRIDSS(
    ch_gridss_reads_input,
    ch_gridss_svs_input,
    ref_data_genome_dir,
    ref_data_genome_fn,
    ref_data_gridss_blacklist,
  )
  ch_versions = ch_versions.mix(GRIDSS.out.versions)

  //
  // MODULE: Run GRIPSS to filter GRIDSS SV calls
  //
  GRIPSS(
    GRIDSS.out.results,
    ref_data_genome_dir,
    ref_data_genome_fn,
    ref_data_gridss_breakend_pon,
    ref_data_gridss_breakpoint_pon,
    ref_data_known_fusions,
  )
  ch_versions = ch_versions.mix(GRIPSS.out.versions)

  //
  // SUBWORKFLOW: Run PURPLE for CNV calling, purity and ploidy inference, SV recovery
  //
  ch_smlv_vcfs = ch_inputs
    .map { meta ->
      [meta, meta.get(['tumor', 'smlv_vcf']), meta.get(['normal', 'smlv_vcf'])]
    }
  ch_purple_inputs = WorkflowHmftools.group_by_meta(
    AMBER.out.amber_dir,
    COBALT.out.cobalt_dir,
    GRIPSS.out.vcf_soft,
    GRIPSS.out.vcf_hard,
    ch_smlv_vcfs,
  )
  PURPLE(
    ch_purple_inputs,
    ref_data_genome_dir,
    ref_data_genome_fn,
    ref_data_cobalt_gc_profile,
    ref_data_known_hotspots,
    ref_data_driver_gene_panel,
    ref_data_ensembl_data_dir,
  )
  ch_versions = ch_versions.mix(PURPLE.out.versions)

  //
  // SUBWORKFLOW: Group structural variants into higher order events with LINX
  //
  LINX(
    PURPLE.out.purple_dir,
    ref_data_linx_fragile_sites,
    ref_data_linx_line_elements,
    ref_data_ensembl_data_dir,
    ref_data_known_fusion_data,
    ref_data_driver_gene_panel,
  )
  ch_versions = ch_versions.mix(LINX.out.versions)

  //
  // MODULE: Run gpgr to generate a LINX report
  //
  LINX_REPORT(
    LINX.out.results,
  )
  ch_versions = ch_versions.mix(LINX_REPORT.out.versions)

  //
  // MODULE: Pipeline reporting
  //
  CUSTOM_DUMPSOFTWAREVERSIONS (
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
