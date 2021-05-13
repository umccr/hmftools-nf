package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.patientreporter.algo.AnalysedReportData;
import com.hartwig.hmftools.patientreporter.algo.ImmutableAnalysedReportData;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReportData;
import com.hartwig.hmftools.patientreporter.summary.SummaryFile;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlackListModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlacklistFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryModel;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PatientReporterTestFactory {

    private static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PIPELINE_VERSION_FILE = RUN_DIRECTORY + "/pipeline.version";
    private static final String PURPLE_PURITY_TSV = RUN_DIRECTORY + "/purple/sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = RUN_DIRECTORY + "/purple/sample.purple.qc";
    private static final String PURPLE_SOMATIC_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/purple/sample.driver.catalog.somatic.tsv";
    private static final String PURPLE_GERMLINE_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/purple/sample.driver.catalog.germline.tsv";
    private static final String PURPLE_SOMATIC_VARIANT_VCF = RUN_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String PURPLE_GERMLINE_VARIANT_VCF = RUN_DIRECTORY + "/purple/sample.purple.germline.vcf";
    private static final String PURPLE_SOMATIC_COPYNUMBER_TSV = RUN_DIRECTORY + "/purple/sample.purple.cnv.somatic.tsv";
    private static final String PURPLE_CIRCOS_FILE = RUN_DIRECTORY + "/purple/plot/sample.circos.png";
    private static final String LINX_FUSIONS_TSV = RUN_DIRECTORY + "/linx/sample.linx.fusion.tsv";
    private static final String LINX_BREAKEND_TSV = RUN_DIRECTORY + "/linx/sample.linx.breakend.tsv";
    private static final String LINX_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/linx/sample.linx.driver.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = RUN_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = RUN_DIRECTORY + "/cuppa/sample.cuppa.conclusion.txt";
    private static final String MOLECULAR_TISSUE_ORIGIN_PLOT = RUN_DIRECTORY + "/cuppa/sample.cuppa.chart.png";
    private static final String VIRUS_BREAKEND_TSV = RUN_DIRECTORY + "/virusbreakend/sample.virusbreakend.vcf.summary.tsv";
    private static final String PEACH_GENOTYPE_TSV = RUN_DIRECTORY + "/peach/sample.peach.genotype.tsv";
    private static final String PROTECT_EVIDENCE_TSV = RUN_DIRECTORY + "/protect/sample.protect.tsv";

    private static final String SIGNATURE_PATH = Resources.getResource("signature/signature_test.png").getPath();
    private static final String RVA_LOGO_PATH = Resources.getResource("rva_logo/rva_logo_test.jpg").getPath();
    private static final String COMPANY_LOGO_PATH = Resources.getResource("company_logo/hartwig_logo_test.jpg").getPath();

    private static final String SAMPLE_SUMMARY_TSV = Resources.getResource("sample_summary/sample_summary.tsv").getPath();
    private static final String VIRUS_DB_TSV = Resources.getResource("virusbreakend/virusdb.tsv").getPath();
    private static final String VIRUS_SUMMARY_TSV = Resources.getResource("virusbreakend/virus_summary.tsv").getPath();
    private static final String VIRUS_BLACKLIST_TSV = Resources.getResource("virusbreakend/virus_blacklist.tsv").getPath();

    private PatientReporterTestFactory() {
    }

    @NotNull
    public static LimsCohortConfig createCPCTCohortConfig() {
        return createCohortConfig("CPCT", true, false, false, false, false, true, false, false, true, false, false);
    }

    @NotNull
    public static LimsCohortConfig createCOLOCohortConfig() {
        return createCohortConfig("COLO", true, false, false, false, false, true, false, false, true, false, false);
    }

    @NotNull
    public static LimsCohortConfig createCORECohortConfig() {
        return createCohortConfig("CORE", true, true, false, true, true, true, true, true, false, true, true);
    }

    @NotNull
    public static LimsCohortConfig createWIDECohortConfig() {
        return createCohortConfig("WIDE", true, true, true, true, true, true, false, true, true, false, true);
    }

    @NotNull
    public static LimsCohortConfig createCOREDBCohortConfig() {
        return createCohortConfig("COREDB", true, true, true, false, true, true, true, true, true, false, true);
    }

    @NotNull
    private static LimsCohortConfig createCohortConfig(@NotNull String cohortId, boolean sampleContainsHospitalCenterId,
            boolean reportGermline, boolean reportGermlineFlag, boolean reportConclusion, boolean reportViral, boolean reportPeach,
            boolean requireHospitalId, boolean requireHospitalPAId, boolean requireHospitalPersonsStudy,
            boolean requireHospitalPersonsRequester, boolean requireAdditionalInformationForSidePanel) {
        return ImmutableLimsCohortConfig.builder()
                .cohortId(cohortId)
                .sampleContainsHospitalCenterId(sampleContainsHospitalCenterId)
                .reportGermline(reportGermline)
                .reportGermlineFlag(reportGermlineFlag)
                .reportConclusion(reportConclusion)
                .reportViral(reportViral)
                .reportPeach(reportPeach)
                .requireHospitalId(requireHospitalId)
                .requireHospitalPAId(requireHospitalPAId)
                .requireHospitalPersonsStudy(requireHospitalPersonsStudy)
                .requireHospitalPersonsRequester(requireHospitalPersonsRequester)
                .requireAdditionalInformationForSidePanel(requireAdditionalInformationForSidePanel)
                .build();
    }

    @NotNull
    public static PatientReporterConfig createTestReporterConfig() {
        return ImmutablePatientReporterConfig.builder()
                .tumorSampleId(Strings.EMPTY)
                .tumorSampleBarcode(Strings.EMPTY)
                .outputDirReport(Strings.EMPTY)
                .outputDirData(Strings.EMPTY)
                .reportingDbTsv(Strings.EMPTY)
                .primaryTumorTsv(Strings.EMPTY)
                .limsDir(Strings.EMPTY)
                .rvaLogo(RVA_LOGO_PATH)
                .companyLogo(COMPANY_LOGO_PATH)
                .signature(SIGNATURE_PATH)
                .qcFail(false)
                .pipelineVersionFile(PIPELINE_VERSION_FILE)
                .purplePurityTsv(PURPLE_PURITY_TSV)
                .purpleQcFile(PURPLE_QC_FILE)
                .purpleSomaticDriverCatalogTsv(PURPLE_SOMATIC_DRIVER_CATALOG_TSV)
                .purpleGermlineDriverCatalogTsv(PURPLE_GERMLINE_DRIVER_CATALOG_TSV)
                .purpleSomaticVariantVcf(PURPLE_SOMATIC_VARIANT_VCF)
                .purpleGermlineVariantVcf(PURPLE_GERMLINE_VARIANT_VCF)
                .purpleSomaticCopyNumberTsv(PURPLE_SOMATIC_COPYNUMBER_TSV)
                .purpleCircosPlot(PURPLE_CIRCOS_FILE)
                .linxFusionTsv(LINX_FUSIONS_TSV)
                .linxBreakendTsv(LINX_BREAKEND_TSV)
                .linxDriverCatalogTsv(LINX_DRIVER_CATALOG_TSV)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .molecularTissueOriginTxt(MOLECULAR_TISSUE_ORIGIN_TXT)
                .molecularTissueOriginPlot(MOLECULAR_TISSUE_ORIGIN_PLOT)
                .virusBreakendTsv(VIRUS_BREAKEND_TSV)
                .peachGenotypeTsv(PEACH_GENOTYPE_TSV)
                .protectEvidenceTsv(PROTECT_EVIDENCE_TSV)
                .germlineReportingTsv(Strings.EMPTY)
                .sampleSummaryTsv(SAMPLE_SUMMARY_TSV)
                .virusDbTsv(VIRUS_DB_TSV)
                .virusSummaryTsv(VIRUS_SUMMARY_TSV)
                .virusBlacklistTsv(VIRUS_BLACKLIST_TSV)
                .isCorrectedReport(false)
                .onlyCreatePDF(false)
                .build();
    }

    @NotNull
    public static ReportData loadTestReportData() {
        List<PatientPrimaryTumor> patientPrimaryTumors = Lists.newArrayList();
        Lims lims = LimsFactory.empty();

        return ImmutableQCFailReportData.builder()
                .patientPrimaryTumors(patientPrimaryTumors)
                .limsModel(lims)
                .signaturePath(SIGNATURE_PATH)
                .logoRVAPath(RVA_LOGO_PATH)
                .logoCompanyPath(COMPANY_LOGO_PATH)
                .build();
    }

    @NotNull
    public static AnalysedReportData loadTestAnalysedReportData() {
        try {
            SummaryModel summaryModel = SummaryFile.buildFromTsv(SAMPLE_SUMMARY_TSV);
            VirusDbModel virusDbModel = VirusDbFile.buildFromTsv(VIRUS_DB_TSV);
            VirusSummaryModel virusSummaryModel = VirusSummaryFile.buildFromTsv(VIRUS_SUMMARY_TSV);
            VirusBlackListModel virusBlackListModel = VirusBlacklistFile.buildFromTsv(VIRUS_BLACKLIST_TSV);

            return ImmutableAnalysedReportData.builder()
                    .from(loadTestReportData())
                    .germlineReportingModel(new GermlineReportingModel(Lists.newArrayList()))
                    .summaryModel(summaryModel)
                    .virusDbModel(virusDbModel)
                    .virusSummaryModel(virusSummaryModel)
                    .virusBlackListModel(virusBlackListModel)
                    .build();
        } catch (IOException exception) {
            throw new IllegalStateException("Could not load test analysed report data: " + exception.getMessage());
        }
    }
}
