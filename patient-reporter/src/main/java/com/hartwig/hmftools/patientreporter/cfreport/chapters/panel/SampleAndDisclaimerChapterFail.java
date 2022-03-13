package com.hartwig.hmftools.patientreporter.cfreport.chapters.panel;

import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.ReportChapter;
import com.hartwig.hmftools.patientreporter.cfreport.components.ReportSignature;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.panel.PanelFailReport;
import com.itextpdf.io.IOException;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.UnitValue;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SampleAndDisclaimerChapterFail implements ReportChapter {

    @NotNull
    private final PanelFailReport report;

    public SampleAndDisclaimerChapterFail(@NotNull PanelFailReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String pdfTitle() {
        return Strings.EMPTY;
    }

    @Override
    @NotNull
    public String name() {
        return "Sample details & disclaimers";
    }

    @Override
    public boolean isFullWidth() {
        return false;
    }

    @Override
    public void render(@NotNull Document reportDocument) throws IOException {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 0.1f, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(createSampleDetailsColumn()));
        table.addCell(TableUtil.createLayoutCell());
        table.addCell(TableUtil.createLayoutCell().add(createDisclaimerColumn()));
        reportDocument.add(table);

        reportDocument.add(ReportSignature.createSignatureDiv(report.logoRVAPath(), report.signaturePath(), true));
        reportDocument.add(ReportSignature.createEndOfReportIndication());
    }

    @NotNull
    private Div createSampleDetailsColumn() {
        SampleReport sampleReport = report.sampleReport();

        Div div = createSampleDetailsDiv();
        div.add(createContentParagraph("The samples have been sequenced at ", ReportResources.HARTWIG_ADDRESS));
        div.add(createContentParagraph("The sample(s) have been analyzed by Next Generation Sequencing using targeted enrichment."));
        div.add(generateHMFSampleIDParagraph(report.sampleReport()));

        String earliestArrivalDate = sampleReport.earliestArrivalDate();
        div.add(createContentParagraphTwice("The results in this report have been obtained between ",
                DataUtil.formatNullableString(earliestArrivalDate),
                " and ",
                report.reportDate()));

        div.add(createContentParagraphTwice("This experiment is performed on the tumor sample which arrived on ",
                DataUtil.formatDate(sampleReport.tumorArrivalDate()),
                " with internal tumor barcode ",
                sampleReport.tumorSampleBarcode()));
        div.add(createContentParagraph("This experiment is performed according to lab procedures: ", sampleReport.labProcedures()));
        String whoVerified = "This report was generated " + report.user();

        div.add(createContentParagraph(whoVerified));
        div.add(createContentParagraph("This report is addressed to: ", sampleReport.addressee()));
        report.comments().ifPresent(comments -> div.add(createContentParagraphRed("Comments: " + comments)));

        return div;
    }

    @NotNull
    private static Paragraph createContentParagraphRed(@NotNull String text) {
        return new Paragraph(text).addStyle(ReportResources.smallBodyTextStyleRed()).setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private static Paragraph generateHMFSampleIDParagraph(@NotNull SampleReport sampleReport) {
        return createContentParagraph("The HMF sample ID is: ", sampleReport.tumorSampleId());
    }

    @NotNull
    private Div createDisclaimerColumn() {
        Div div = createDisclaimerDiv();

        String pipelineVersion = "" == null ? "No pipeline version is known" : "";
        div.add(createContentParagraph("This report is based on pipeline version ", pipelineVersion + "."));
        div.add(createContentParagraph("No check is performed to verify the ‘primary tumor location’ and ‘primary tumor type’ information."));
        div.add(createContentParagraph("The results of this report and result file is based solely on the results of the DNA sequencing of "
                + "the received tumor material"));
        div.add(createContentParagraph("Any clinical Interpretation of the result file is the responsibility of the hospital."));
        div.add(createContentParagraph("Based on a tumor purity of at least 5%, the test has a sensitivity of >95% for detection of "
                + "somatic variants and >95% for detec=on of translocations and gene copy number changes."));
        div.add(createContentParagraph("For feedback or complaints please contact ", ReportResources.CONTACT_EMAIL_QA + "."));
        div.add(createContentParagraph("For questions about the contents of this report, please contact ",
                ReportResources.CONTACT_EMAIL_GENERAL + "."));

        return div;
    }

    @NotNull
    private static Div createDisclaimerDiv() {
        Div div = new Div();
        div.add(new Paragraph("Disclaimer").addStyle(ReportResources.smallBodyHeadingStyle()));
        return div;
    }

    @NotNull
    private static Div createSampleDetailsDiv() {
        Div div = new Div();
        div.add(new Paragraph("Sample details").addStyle(ReportResources.smallBodyHeadingStyle()));
        return div;
    }

    @NotNull
    private static Paragraph createContentParagraph(@NotNull String regularPart, @NotNull String boldPart) {
        return createContentParagraph(regularPart).add(new Text(boldPart).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private static Paragraph createContentParagraphTwice(@NotNull String regularPart, @NotNull String boldPart,
            @NotNull String regularPart2, @NotNull String boldPart2) {
        return createContentParagraph(regularPart).add(new Text(boldPart).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .add(regularPart2)
                .add(new Text(boldPart2).addStyle(ReportResources.smallBodyBoldTextStyle()))
                .setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }

    @NotNull
    private static Paragraph createContentParagraph(@NotNull String text) {
        return new Paragraph(text).addStyle(ReportResources.smallBodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING);
    }
}

