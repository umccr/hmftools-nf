package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.ClonalityCutoffKernel;
import com.hartwig.hmftools.common.variant.ClonalityFactory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.BaseReporterData;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.ImmutableAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.copynumber.ImmutablePurpleAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.SomaticVariantAnalyzer;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientReporter {
    private static final Logger LOGGER = LogManager.getLogger(PatientReporter.class);

    @NotNull
    public abstract BaseReporterData baseReporterData();

    @NotNull
    public abstract HmfReporterData reporterData();

    @NotNull
    public abstract StructuralVariantAnalyzer structuralVariantAnalyzer();

    @NotNull
    public AnalysedPatientReport run(@NotNull final String runDirectory, @Nullable final String comments) throws IOException {
        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDirectory);
        assert run.isSomaticRun();

        final PurpleAnalysis purpleAnalysis = analyzePurpleCopyNumbers(run, reporterData().panelGeneModel().panel());

        final SomaticVariantAnalysis somaticVariantAnalysis = analyzeSomaticVariants(run,
                purpleAnalysis,
                reporterData().panelGeneModel().panel(),
                reporterData().highConfidenceRegions(),
                reporterData().refGenomeFastaFile());

        final StructuralVariantAnalysis structuralVariantAnalysis =
                analyzeStructuralVariants(run, purpleAnalysis, structuralVariantAnalyzer());
        final List<GeneFusion> reportableFusions = structuralVariantAnalysis.reportableFusions();
        final List<GeneDisruption> reportableDisruptions = structuralVariantAnalysis.reportableDisruptions();

        LOGGER.info("Printing analysis results:");
        LOGGER.info(" Number of somaticVariants to report : " + Integer.toString(somaticVariantAnalysis.variantsToReport().size()));
        LOGGER.info(" Microsatellite analysis results: " + Double.toString(somaticVariantAnalysis.indelsPerMb()) + " indels per MB");
        LOGGER.info(" Mutational load results: " + Integer.toString(somaticVariantAnalysis.mutationalLoad()));
        LOGGER.info(" Number of copy number events to report: " + Integer.toString(purpleAnalysis.reportableGeneCopyNumbers().size()));
        LOGGER.info(" Number of gene fusions to report : " + Integer.toString(reportableFusions.size()));
        LOGGER.info(" Number of gene disruptions to report : " + Integer.toString(reportableDisruptions.size()));

        final String tumorSample = run.tumorSample();
        final SampleReport sampleReport = ImmutableSampleReport.of(tumorSample,
                PatientReporterHelper.extractPatientTumorLocation(baseReporterData().patientTumorLocations(), tumorSample),
                baseReporterData().limsModel().tumorPercentageForSample(tumorSample),
                baseReporterData().limsModel().arrivalDateForSample(tumorSample),
                baseReporterData().limsModel().arrivalDateForSample(run.refSample()),
                baseReporterData().limsModel().labProceduresForSample(tumorSample),
                baseReporterData().centerModel().getAddresseeStringForSample(tumorSample));

        return ImmutableAnalysedPatientReport.of(sampleReport,
                purpleAnalysis.status(),
                purpleAnalysis.fittedPurity().purity(),
                somaticVariantAnalysis.variantsToReport(),
                somaticVariantAnalysis.mutationalLoad(),
                somaticVariantAnalysis.indelsPerMb(),
                purpleAnalysis.reportableGeneCopyNumbers(),
                reportableFusions,
                reportableDisruptions,
                PatientReporterHelper.findCircosPlotPath(runDirectory, tumorSample),
                Optional.ofNullable(comments),
                baseReporterData().signaturePath());
    }

    @NotNull
    private static PurpleAnalysis analyzePurpleCopyNumbers(@NotNull RunContext run, @NotNull Set<String> genePanel) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        final PurityContext purityContext = PatientReporterHelper.loadPurity(runDirectory, sample);

        final List<PurpleCopyNumber> purpleCopyNumbers = PatientReporterHelper.loadPurpleCopyNumbers(runDirectory, sample);
        LOGGER.info(" " + purpleCopyNumbers.size() + " purple copy number regions loaded for sample " + sample);

        final List<GeneCopyNumber> panelGeneCopyNumbers = PatientReporterHelper.loadPurpleGeneCopyNumbers(runDirectory, sample)
                .stream()
                .filter(geneCopyNumber -> genePanel.contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        return ImmutablePurpleAnalysis.builder()
                .gender(purityContext.gender())
                .status(purityContext.status())
                .fittedPurity(purityContext.bestFit())
                .fittedScorePurity(purityContext.score())
                .copyNumbers(purpleCopyNumbers)
                .panelGeneCopyNumbers(panelGeneCopyNumbers)
                .build();
    }

    @NotNull
    private static SomaticVariantAnalysis analyzeSomaticVariants(@NotNull RunContext run, @NotNull PurpleAnalysis purpleAnalysis,
            @NotNull Set<String> genePanel, @NotNull Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull IndexedFastaSequenceFile refGenomeFastaFile) throws IOException {
        final String runDirectory = run.runDirectory();
        final String sample = run.tumorSample();

        LOGGER.info("Loading somatic somaticVariants...");
        final List<SomaticVariant> variants = PatientReporterHelper.loadPassedSomaticVariants(sample, runDirectory);
        LOGGER.info(" " + variants.size() + " PASS somatic somaticVariants loaded for sample " + sample);

        LOGGER.info("Enriching somatic somaticVariants");
        final List<EnrichedSomaticVariant> enrichedSomaticVariants =
                enrich(variants, purpleAnalysis, highConfidenceRegions, refGenomeFastaFile);

        LOGGER.info("Analyzing somatic somaticVariants....");
        return SomaticVariantAnalyzer.run(enrichedSomaticVariants, genePanel);
    }

    @NotNull
    private static List<EnrichedSomaticVariant> enrich(@NotNull List<SomaticVariant> variants, @NotNull PurpleAnalysis purpleAnalysis,
            @NotNull Multimap<String, GenomeRegion> highConfidenceRegions, @NotNull IndexedFastaSequenceFile refGenomeFastaFile) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(purpleAnalysis.gender(), purpleAnalysis.fittedPurity());
        final PurityAdjustedSomaticVariantFactory purityAdjustedFactory =
                new PurityAdjustedSomaticVariantFactory(purityAdjuster, purpleAnalysis.copyNumbers(), Collections.emptyList());
        final List<PurityAdjustedSomaticVariant> purityAdjustedSomaticVariants = purityAdjustedFactory.create(variants);

        final double clonalPloidy = ClonalityCutoffKernel.clonalCutoff(purityAdjustedSomaticVariants);
        final ClonalityFactory clonalityFactory = new ClonalityFactory(purityAdjuster, clonalPloidy);

        final EnrichedSomaticVariantFactory enrichedSomaticFactory =
                new EnrichedSomaticVariantFactory(highConfidenceRegions, refGenomeFastaFile, clonalityFactory);

        return enrichedSomaticFactory.enrich(purityAdjustedSomaticVariants);
    }

    @NotNull
    private static StructuralVariantAnalysis analyzeStructuralVariants(@NotNull RunContext run, @NotNull PurpleAnalysis purpleAnalysis,
            @NotNull StructuralVariantAnalyzer structuralVariantAnalyzer) throws IOException {
        final Path structuralVariantVCF = PatientReporterHelper.findStructuralVariantVCF(run.runDirectory());
        LOGGER.info("Loading structural somaticVariants...");
        final List<StructuralVariant> structuralVariants = StructuralVariantFileLoader.fromFile(structuralVariantVCF.toString(), true);

        LOGGER.info("Enriching structural somaticVariants with purple data.");
        final PurityAdjuster purityAdjuster = new PurityAdjuster(purpleAnalysis.gender(), purpleAnalysis.fittedPurity());
        final Multimap<String, PurpleCopyNumber> copyNumberMap =
                Multimaps.index(purpleAnalysis.copyNumbers(), PurpleCopyNumber::chromosome);

        final List<EnrichedStructuralVariant> enrichedStructuralVariants =
                EnrichedStructuralVariantFactory.enrich(structuralVariants, purityAdjuster, copyNumberMap);

        LOGGER.info("Analysing structural somaticVariants...");
        return structuralVariantAnalyzer.run(enrichedStructuralVariants);
    }
}
