package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaData;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaInterpretation;
import com.hartwig.hmftools.orange.algo.cuppa.CuppaPrediction;
import com.hartwig.hmftools.orange.algo.purple.DriverInterpretation;
import com.hartwig.hmftools.orange.algo.purple.PurpleCharacteristics;
import com.hartwig.hmftools.orange.algo.purple.PurpleGainLoss;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;
import com.hartwig.hmftools.orange.cohort.datamodel.Evaluation;
import com.hartwig.hmftools.orange.cohort.mapping.CohortConstants;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Drivers;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Images;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FrontPageChapter implements ReportChapter {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat TWO_DIGITS = ReportResources.decimalFormat("#.##");
    private static final DecimalFormat PERCENTAGE = ReportResources.decimalFormat("#'%'");

    private static final String NONE = "None";

    @NotNull
    private final OrangeReport report;
    @NotNull
    private final PlotPathResolver plotPathResolver;

    public FrontPageChapter(@NotNull final OrangeReport report, @NotNull final PlotPathResolver plotPathResolver) {
        this.report = report;
        this.plotPathResolver = plotPathResolver;
    }

    @NotNull
    @Override
    public String name() {
        return "Front Page";
    }

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull Document document) {
        addSummaryTable(document);
        addDetailsAndPlots(document);
    }

    private void addSummaryTable(@NotNull Document document) {
        Table table = Tables.createContent(contentWidth(),
                new float[] { 3, 2, 1 },
                new Cell[] { Cells.createHeader("Configured Primary Tumor"), Cells.createHeader("Cuppa Cancer Type"),
                        Cells.createHeader("QC") });

        table.addCell(Cells.createContent(configuredPrimaryTumor(report.configuredPrimaryTumor())));
        table.addCell(Cells.createContent(cuppaCancerType(report.cuppa())));
        table.addCell(Cells.createContent(purpleQCString()));
        document.add(Tables.createWrapping(table));
    }

    @NotNull
    private static String cuppaCancerType(@Nullable CuppaData cuppa) {
        if (cuppa == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        CuppaPrediction best = CuppaInterpretation.best(cuppa);
        return best.cancerType() + " (" + PERCENTAGE.format(best.likelihood() * 100) + ")";
    }

    @NotNull
    private static String configuredPrimaryTumor(@NotNull Set<DoidNode> nodes) {
        Set<String> configured = Sets.newHashSet();
        for (DoidNode node : nodes) {
            configured.add(node.doidTerm() + " (DOID " + node.doid() + ")");
        }

        return concat(configured);
    }

    @NotNull
    private String purpleQCString() {
        Set<String> purpleStatuses = Sets.newHashSet();
        for (PurpleQCStatus status : report.purple().fit().qc().status()) {
            purpleStatuses.add(status.toString());
        }
        return concat(purpleStatuses);
    }

    private void addDetailsAndPlots(@NotNull Document document) {
        Table topTable = new Table(UnitValue.createPercentArray(new float[] { 1, 1 })).setWidth(contentWidth() - 5);

        Table summary = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        summary.addCell(Cells.createKey("Purity:"));
        summary.addCell(Cells.createValue(purityString()));
        summary.addCell(Cells.createKey("Ploidy:"));
        summary.addCell(Cells.createValue(ploidyString()));
        summary.addCell(Cells.createKey("Somatic variant drivers:"));
        summary.addCell(Cells.createValue(somaticVariantDriverString()));
        summary.addCell(Cells.createKey("Germline variant drivers:"));
        summary.addCell(Cells.createValue(germlineVariantDriverString()));
        summary.addCell(Cells.createKey("Somatic copy number drivers:"));
        summary.addCell(Cells.createValue(somaticCopyNumberDriverString()));
        summary.addCell(Cells.createKey("Germline copy number drivers:"));
        summary.addCell(Cells.createValue(germlineCopyNumberDriverString()));
        summary.addCell(Cells.createKey("Somatic disruption drivers:"));
        summary.addCell(Cells.createValue(somaticDisruptionDriverString()));
        summary.addCell(Cells.createKey("Germline disruption drivers:"));
        summary.addCell(Cells.createValue(germlineDisruptionDriverString()));
        summary.addCell(Cells.createKey("Fusion drivers:"));
        summary.addCell(Cells.createValue(fusionDriverString()));
        summary.addCell(Cells.createKey("Viral presence:"));
        summary.addCell(Cells.createValue(virusString()));
        summary.addCell(Cells.createKey("Whole genome duplicated:"));
        summary.addCell(Cells.createValue(report.purple().characteristics().wholeGenomeDuplication() ? "Yes" : "No"));
        summary.addCell(Cells.createKey("Microsatellite indels per Mb:"));
        summary.addCell(Cells.createValue(msiString()));
        summary.addCell(Cells.createKey("Tumor mutations per Mb:"));
        summary.addCell(Cells.createValue(SINGLE_DIGIT.format(report.purple().characteristics().tumorMutationalBurdenPerMb())));
        summary.addCell(Cells.createKey("Tumor mutational load:"));
        summary.addCell(Cells.createValue(tmlString()));
        summary.addCell(Cells.createKey("HR deficiency score:"));
        summary.addCell(Cells.createValue(hrDeficiencyString()));
        summary.addCell(Cells.createKey("DPYD status:"));
        summary.addCell(Cells.createValue(dpydStatus()));
        summary.addCell(Cells.createKey("Number of SVs:"));
        summary.addCell(Cells.createValue(svTmbString()));
        summary.addCell(Cells.createKey("Max complex cluster size:"));
        summary.addCell(Cells.createValue(maxComplexSizeString()));
        summary.addCell(Cells.createKey("Telomeric SGLs:"));
        summary.addCell(Cells.createValue(telomericSGLString()));
        summary.addCell(Cells.createKey("Number of LINE insertions:"));
        summary.addCell(Cells.createValue(lineCountString()));

        Image circosImage = Images.build(plotPathResolver.resolve(report.plots().purpleFinalCircosPlot()));
        circosImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        circosImage.setMaxHeight(290);

        topTable.addCell(summary);
        topTable.addCell(circosImage);

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(contentWidth()).setPadding(0);
        table.addCell(topTable);

        Image clonalityImage = Images.build(plotPathResolver.resolve(report.plots().purpleClonalityPlot()));
        clonalityImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        clonalityImage.setMaxHeight(270);

        table.addCell(clonalityImage);
        document.add(table);
    }

    @NotNull
    private String purityString() {
        return String.format("%s (%s-%s)",
                PERCENTAGE.format(report.purple().fit().purity() * 100),
                PERCENTAGE.format(report.purple().fit().minPurity() * 100),
                PERCENTAGE.format(report.purple().fit().maxPurity() * 100));
    }

    @NotNull
    private String ploidyString() {
        return String.format("%s (%s-%s)",
                TWO_DIGITS.format(report.purple().fit().ploidy()),
                TWO_DIGITS.format(report.purple().fit().minPloidy()),
                TWO_DIGITS.format(report.purple().fit().maxPloidy()));
    }

    @NotNull
    private String somaticVariantDriverString() {
        return variantDriverString(report.purple().reportableSomaticVariants(), report.purple().somaticDrivers());
    }

    @NotNull
    private String germlineVariantDriverString() {
        List<PurpleVariant> reportableGermlineVariants = report.purple().reportableGermlineVariants();
        List<DriverCatalog> germlineDrivers = report.purple().germlineDrivers();
        if (reportableGermlineVariants != null && germlineDrivers != null) {
            return variantDriverString(reportableGermlineVariants, germlineDrivers);
        } else {
            return ReportResources.NOT_AVAILABLE;
        }
    }

    @NotNull
    private static String variantDriverString(@NotNull List<PurpleVariant> variants, @NotNull List<DriverCatalog> drivers) {
        if (variants.isEmpty()) {
            return NONE;
        }

        Set<String> highDriverGenes = Sets.newTreeSet(Comparator.naturalOrder());
        for (PurpleVariant variant : variants) {
            DriverCatalog driver = Drivers.canonicalMutationEntryForGene(drivers, variant.gene());
            if (driver != null && DriverInterpretation.interpret(driver.driverLikelihood()) == DriverInterpretation.HIGH) {
                highDriverGenes.add(variant.gene());
            }
        }

        return !highDriverGenes.isEmpty() ? variants.size() + " (" + concat(highDriverGenes) + ")" : String.valueOf(variants.size());
    }

    @NotNull
    private String somaticCopyNumberDriverString() {
        return copyNumberDriverString(report.purple().reportableSomaticGainsLosses());
    }

    @NotNull
    private String germlineCopyNumberDriverString() {
        List<PurpleGainLoss> germlineGainsLosses = report.purple().reportableGermlineFullLosses();
        if (germlineGainsLosses == null) {
            return ReportResources.NOT_AVAILABLE;
        }
        return copyNumberDriverString(germlineGainsLosses);
    }

    @NotNull
    private static String copyNumberDriverString(@NotNull List<PurpleGainLoss> gainsLosses) {
        if (gainsLosses.isEmpty()) {
            return NONE;
        }

        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        for (PurpleGainLoss gainLoss : gainsLosses) {
            genes.add(gainLoss.gene());
        }
        return gainsLosses.size() + " (" + concat(genes) + ")";
    }

    @NotNull
    private String somaticDisruptionDriverString() {
        return disruptionDriverString(report.linx().somaticHomozygousDisruptions());
    }

    @NotNull
    private String germlineDisruptionDriverString() {
        List<HomozygousDisruption> germlineHomozygousDisruptions = report.linx().germlineHomozygousDisruptions();
        if (germlineHomozygousDisruptions == null) {
            return ReportResources.NOT_AVAILABLE;
        }
        return disruptionDriverString(germlineHomozygousDisruptions);
    }

    @NotNull
    private static String disruptionDriverString(@NotNull List<HomozygousDisruption> homozygousDisruptions) {
        if (homozygousDisruptions.isEmpty()) {
            return NONE;
        }

        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        for (HomozygousDisruption homozygousDisruption : homozygousDisruptions) {
            genes.add(homozygousDisruption.gene());
        }
        return homozygousDisruptions.size() + " (" + concat(genes) + ")";
    }

    @NotNull
    private String fusionDriverString() {
        if (report.linx().reportableSomaticFusions().isEmpty()) {
            return NONE;
        }

        Set<String> fusions = Sets.newTreeSet(Comparator.naturalOrder());
        for (LinxFusion fusion : report.linx().reportableSomaticFusions()) {
            fusions.add(fusion.name());
        }
        return report.linx().reportableSomaticFusions().size() + " (" + concat(fusions) + ")";
    }

    @NotNull
    private String virusString() {
        VirusInterpreterData virusInterpreter = report.virusInterpreter();
        if (virusInterpreter == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        if (virusInterpreter.reportableViruses().isEmpty()) {
            return NONE;
        }

        Set<String> viruses = Sets.newTreeSet(Comparator.naturalOrder());
        for (AnnotatedVirus virus : virusInterpreter.reportableViruses()) {
            if (virus.interpretation() != null) {
                viruses.add(virus.interpretation());
            } else {
                viruses.add(virus.name());
            }
        }

        return virusInterpreter.reportableViruses().size() + " (" + concat(viruses) + ")";
    }

    @NotNull
    private String msiString() {
        PurpleCharacteristics characteristics = report.purple().characteristics();
        return SINGLE_DIGIT.format(characteristics.microsatelliteIndelsPerMb()) + " (" + characteristics.microsatelliteStatus().display()
                + ")";
    }

    @NotNull
    private String tmlString() {
        PurpleCharacteristics characteristics = report.purple().characteristics();
        return characteristics.tumorMutationalLoad() + " (" + characteristics.tumorMutationalLoadStatus().display() + ")";
    }

    @NotNull
    private String hrDeficiencyString() {
        ChordData chord = report.chord();
        if (chord == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        if (chord.hrStatus() == ChordStatus.CANNOT_BE_DETERMINED) {
            return ChordStatus.CANNOT_BE_DETERMINED.display();
        }

        String addon = Strings.EMPTY;
        if (chord.hrStatus() == ChordStatus.HR_DEFICIENT) {
            if (chord.hrdType().contains("BRCA1")) {
                addon = " - BRCA1 (" + TWO_DIGITS.format(chord.BRCA1Value()) + ")";
            } else if (chord.hrdType().contains("BRCA2")) {
                addon = " - BRCA2 (" + TWO_DIGITS.format(chord.BRCA2Value()) + ")";
            } else {
                addon = chord.hrdType();
            }
        }
        return SINGLE_DIGIT.format(chord.hrdValue()) + " (" + chord.hrStatus().display() + addon + ")";
    }

    @NotNull
    private String dpydStatus() {
        List<PeachGenotype> genotypes = report.peach();
        if (genotypes == null ) {
            return ReportResources.NOT_AVAILABLE;
        }

        Set<String> haplotypes = Sets.newHashSet();
        for (PeachGenotype genotype : genotypes) {
            if (genotype.gene().equals("DPYD")) {
                haplotypes.add(genotype.haplotype() + " (" + genotype.function() + ")");
            }
        }
        return !haplotypes.isEmpty() ? concat(haplotypes) : NONE;
    }

    @NotNull
    private String svTmbString() {
        String svTmb = String.valueOf(report.purple().characteristics().svTumorMutationalBurden());

        Evaluation evaluation = report.cohortEvaluations().get(PercentileType.SV_TMB);
        String addon = Strings.EMPTY;
        if (evaluation != null) {
            String panCancerPercentile = PERCENTAGE.format(evaluation.panCancerPercentile() * 100);
            addon = " (Pan " + panCancerPercentile;
            String cancerType = evaluation.cancerType();
            if (cancerType != null && !cancerType.equals(CohortConstants.COHORT_OTHER)
                    && !cancerType.equals(CohortConstants.COHORT_UNKNOWN)) {
                Double percentile = evaluation.cancerTypePercentile();
                String cancerTypePercentile = percentile != null ? PERCENTAGE.format(evaluation.cancerTypePercentile() * 100) : "NA";
                addon = addon + " | " + evaluation.cancerType() + " " + cancerTypePercentile;
            }
            addon = addon + ")";
        }

        return svTmb + addon;
    }

    @NotNull
    private String maxComplexSizeString() {
        CuppaData cuppa = report.cuppa();
        return cuppa != null ? Integer.toString(cuppa.maxComplexSize()) : ReportResources.NOT_AVAILABLE;
    }

    @NotNull
    private String telomericSGLString() {
        CuppaData cuppa = report.cuppa();
        return cuppa != null ? Integer.toString(cuppa.telomericSGLs()) : ReportResources.NOT_AVAILABLE;
    }

    @NotNull
    private String lineCountString() {
        CuppaData cuppa = report.cuppa();
        return cuppa != null ? Integer.toString(cuppa.LINECount()) : ReportResources.NOT_AVAILABLE;
    }

    @NotNull
    private static String concat(@NotNull Iterable<String> strings) {
        StringJoiner joiner = new StringJoiner(", ");
        for (String string : strings) {
            joiner.add(string);
        }

        return joiner.toString();
    }
}
