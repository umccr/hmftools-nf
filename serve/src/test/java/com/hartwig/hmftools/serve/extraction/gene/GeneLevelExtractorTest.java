package com.hartwig.hmftools.serve.extraction.gene;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.DriverGeneTestFactory;
import com.hartwig.hmftools.serve.extraction.catalog.DealWithDriverInconsistentModeAnnotation;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResourceTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneLevelExtractorTest {

    private static final GeneChecker GENE_CHECKER =
            new GeneChecker(Sets.newHashSet("KIT", "NTRK3", "STK11", "MET", "TP53", "KRAS", "NOTCH1", "BRCA1", "NTRK3", "BRAF"));

    @Test
    public void canFilterInCatalogFusion() {
        GeneLevelExtractor geneLevelExtractorIgnore = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertNotNull(geneLevelExtractorIgnore.extract("NTRK3", EventType.PROMISCUOUS_FUSION, "NTRK3  fusion"));

        GeneLevelExtractor geneLevelExtractorWarn = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.WARN_ONLY);
        assertNotNull(geneLevelExtractorWarn.extract("NTRK3", EventType.PROMISCUOUS_FUSION, "NTRK3  fusion"));

        GeneLevelExtractor geneLevelExtractorFilter = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.FILTER);
        assertNotNull(geneLevelExtractorFilter.extract("NTRK3", EventType.PROMISCUOUS_FUSION, "NTRK3  fusion"));
    }

    @Test
    public void canFilterNotInCatalogFusion() {
        GeneLevelExtractor geneLevelExtractorIgnore = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertNotNull(geneLevelExtractorIgnore.extract("KIT", EventType.PROMISCUOUS_FUSION, "KIT  fusion"));

        GeneLevelExtractor geneLevelExtractorWarn = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.WARN_ONLY);
        assertNotNull(geneLevelExtractorWarn.extract("KIT", EventType.PROMISCUOUS_FUSION, "KIT  fusion"));

        GeneLevelExtractor geneLevelExtractorFilter = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.FILTER);
        assertNull(geneLevelExtractorFilter.extract("KIT", EventType.PROMISCUOUS_FUSION, "KIT  fusion"));
    }

    @Test
    public void canFilterInCatalogWildType() {
        GeneLevelExtractor geneLevelExtractorIgnore = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertNotNull(geneLevelExtractorIgnore.extract("KIT", EventType.WILD_TYPE, "KIT  wild type"));

        GeneLevelExtractor geneLevelExtractorWarn = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.WARN_ONLY);
        assertNotNull(geneLevelExtractorWarn.extract("KIT", EventType.WILD_TYPE, "KIT  wild type"));

        GeneLevelExtractor geneLevelExtractorFilter = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.FILTER);
        assertNotNull(geneLevelExtractorFilter.extract("KIT", EventType.WILD_TYPE, "KIT  wild type"));
    }

    @Test
    public void canFilterNotInCatalogWildType() {
        GeneLevelExtractor geneLevelExtractorIgnore = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KRAS"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertNotNull(geneLevelExtractorIgnore.extract("BRAF", EventType.WILD_TYPE, "BRAF  wild type"));

        GeneLevelExtractor geneLevelExtractorWarn = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KRAS"),
                DealWithDriverInconsistentModeAnnotation.WARN_ONLY);
        assertNotNull(geneLevelExtractorWarn.extract("BRAF", EventType.WILD_TYPE, "BRAF  wild type"));

        GeneLevelExtractor geneLevelExtractorFilter = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KRAS"),
                DealWithDriverInconsistentModeAnnotation.FILTER);
        assertNull(geneLevelExtractorFilter.extract("BRAF", EventType.WILD_TYPE, "BRAF  wild type"));
    }

    @Test
    public void canFilterInCatalogGeneLevel() {
        GeneLevelExtractor geneLevelExtractorIgnore = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertNotNull(geneLevelExtractorIgnore.extract("KIT", EventType.GENE_LEVEL, "KIT  mutant"));

        GeneLevelExtractor geneLevelExtractorWarnMis = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.WARN_ONLY);
        assertNotNull(geneLevelExtractorWarnMis.extract("STK11", EventType.GENE_LEVEL, "STK11  act mut"));

        GeneLevelExtractor geneLevelExtractorFilter = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.FILTER);
        assertNotNull(geneLevelExtractorFilter.extract("KIT", EventType.GENE_LEVEL, "KIT  mutant"));

        GeneLevelExtractor geneLevelExtractorFilterMismatch = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.FILTER);
        assertNull(geneLevelExtractorFilterMismatch.extract("STK11", EventType.GENE_LEVEL, "STK11  act mut"));
    }

    @Test
    public void canFilterNotInCatalogGeneLevel() {
        GeneLevelExtractor geneLevelExtractorIgnore = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertNotNull(geneLevelExtractorIgnore.extract("MET", EventType.GENE_LEVEL, "MET  activation"));

        GeneLevelExtractor geneLevelExtractorWarn = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.WARN_ONLY);
        assertNotNull(geneLevelExtractorWarn.extract("MET", EventType.GENE_LEVEL, "MET  activation"));

        GeneLevelExtractor geneLevelExtractorFilter = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.FILTER);
        assertNull(geneLevelExtractorFilter.extract("MET", EventType.GENE_LEVEL, "MET  activation"));
    }

    @Test
    public void canExtractGeneLevelEventWiltType() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        GeneLevelAnnotation geneLevelEvent = geneLevelExtractor.extract("KIT", EventType.WILD_TYPE, "KIT  wild type");

        assertNotNull(geneLevelEvent);
        assertEquals("KIT", geneLevelEvent.gene());
        assertEquals(GeneLevelEvent.WILD_TYPE, geneLevelEvent.event());
    }

    @Test
    public void canExtractGeneLevelEventOnco() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        GeneLevelAnnotation geneLevelEvent = geneLevelExtractor.extract("KIT", EventType.GENE_LEVEL, "KIT  positive");

        assertNotNull(geneLevelEvent);
        assertEquals("KIT", geneLevelEvent.gene());
        assertEquals(GeneLevelEvent.ACTIVATION, geneLevelEvent.event());
    }

    @Test
    public void canExtractGeneLevelEventTsg() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("TP53", "KIT"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        GeneLevelAnnotation geneLevelEvent = geneLevelExtractor.extract("TP53", EventType.GENE_LEVEL, "TP53  negative");

        assertNotNull(geneLevelEvent);
        assertEquals("TP53", geneLevelEvent.gene());
        assertEquals(GeneLevelEvent.INACTIVATION, geneLevelEvent.event());
    }

    @Test
    public void pickEventClassificationOnConflict() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "KIT"),
                DealWithDriverInconsistentModeAnnotation.FILTER);

        GeneLevelAnnotation conflictingGeneLevelEvent = geneLevelExtractor.extract("STK11", EventType.GENE_LEVEL, "STK11 positive");
        assertNull(conflictingGeneLevelEvent);
    }

    @Test
    public void canExtractGeneLevelEventGeneral() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "MET"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        GeneLevelAnnotation geneLevelEvent = geneLevelExtractor.extract("STK11", EventType.GENE_LEVEL, "Truncating Mutations");

        assertNotNull(geneLevelEvent);
        assertEquals("STK11", geneLevelEvent.gene());
        assertEquals(GeneLevelEvent.INACTIVATION, geneLevelEvent.event());
    }

    @Test
    public void canExtractGeneLevelEventFusion() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "MET"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        GeneLevelAnnotation geneLevelEvent = geneLevelExtractor.extract("NTRK3", EventType.PROMISCUOUS_FUSION, "NTRK3 fusion");

        assertNotNull(geneLevelEvent);
        assertEquals("NTRK3", geneLevelEvent.gene());
        assertEquals(GeneLevelEvent.FUSION, geneLevelEvent.event());
    }

    @Test
    public void filtersNonExistingGenes() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("STK11", "MET"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertNull(geneLevelExtractor.extract("NOT-A-GENE", EventType.PROMISCUOUS_FUSION, "NTRK3 fusion"));
    }

    @Test
    public void canExtractGeneLevelEvent() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("NOTCH1", "MET"),
                DealWithDriverInconsistentModeAnnotation.FILTER);

        assertEquals(ImmutableGeneLevelAnnotation.builder().gene("MET").event(GeneLevelEvent.ACTIVATION).build(),
                geneLevelExtractor.extractGeneLevelEvent("MET", "MET activating mutation"));
        assertEquals(ImmutableGeneLevelAnnotation.builder().gene("MET").event(GeneLevelEvent.ACTIVATION).build(),
                geneLevelExtractor.extractGeneLevelEvent("MET", "MET act mut"));
        assertEquals(ImmutableGeneLevelAnnotation.builder().gene("NOTCH1").event(GeneLevelEvent.INACTIVATION).build(),
                geneLevelExtractor.extractGeneLevelEvent("NOTCH1", "LOSS-OF-FUNCTION"));
        assertEquals(ImmutableGeneLevelAnnotation.builder().gene("NOTCH1").event(GeneLevelEvent.INACTIVATION).build(),
                geneLevelExtractor.extractGeneLevelEvent("NOTCH1", "inact mut"));

        assertEquals(ImmutableGeneLevelAnnotation.builder().gene("MET").event(GeneLevelEvent.ACTIVATION).build(),
                geneLevelExtractor.extractGeneLevelEvent("MET", "MUTATION"));
        assertEquals(ImmutableGeneLevelAnnotation.builder().gene("NOTCH1").event(GeneLevelEvent.INACTIVATION).build(),
                geneLevelExtractor.extractGeneLevelEvent("NOTCH1", "MUTATION"));
        assertEquals(ImmutableGeneLevelAnnotation.builder().gene("NOTCH1").event(GeneLevelEvent.INACTIVATION).build(),
                geneLevelExtractor.extractGeneLevelEvent("NOTCH1", "NOTCH1 "));
        assertNull(geneLevelExtractor.extractGeneLevelEvent("BRCA1", "BRCA1"));
        assertNull(geneLevelExtractor.extractGeneLevelEvent("KRAS", "not a gene level event"));
    }

    @Test
    public void canExtractWildTypeEvents() {
        GeneLevelExtractor geneLevelExtractor = createWithDriverGenes(DriverGeneTestFactory.createDriverGenes("NOTCH1", "MET"),
                DealWithDriverInconsistentModeAnnotation.IGNORE);
        assertEquals(ImmutableGeneLevelAnnotation.builder().gene("MET").event(GeneLevelEvent.WILD_TYPE).build(),
                geneLevelExtractor.extractWildTypeEvents("MET", EventType.WILD_TYPE));
        assertEquals(ImmutableGeneLevelAnnotation.builder().gene("AB").event(GeneLevelEvent.WILD_TYPE).build(),
                geneLevelExtractor.extractWildTypeEvents("AB", EventType.WILD_TYPE));
        assertNotEquals(ImmutableGeneLevelAnnotation.builder().gene("AB").event(GeneLevelEvent.WILD_TYPE).build(),
                geneLevelExtractor.extractWildTypeEvents("TP53", EventType.WILD_TYPE));
    }

    @Test
    public void canDetermineGeneLevelFromDriverGenes() {
        List<DriverGene> driverGenes = DriverGeneTestFactory.createDriverGenes("STK11", "MET");

        assertEquals(GeneLevelEvent.ACTIVATION, GeneLevelExtractor.determineGeneLevelEventFromDriverGenes(driverGenes, "MET"));
        assertEquals(GeneLevelEvent.INACTIVATION, GeneLevelExtractor.determineGeneLevelEventFromDriverGenes(driverGenes, "STK11"));
        assertEquals(GeneLevelEvent.ANY_MUTATION, GeneLevelExtractor.determineGeneLevelEventFromDriverGenes(driverGenes, "MAP1K1"));
    }

    @NotNull
    private static GeneLevelExtractor createWithDriverGenes(@NotNull List<DriverGene> driverGenes,
            @NotNull DealWithDriverInconsistentModeAnnotation annotation) {
        return new GeneLevelExtractor(GENE_CHECKER,
                GENE_CHECKER,
                driverGenes,
                RefGenomeResourceTestFactory.buildTestResource37().knownFusionCache(),
                Sets.newHashSet("positive", "activating mutation", "act mut"),
                Sets.newHashSet("negative", "LOSS-OF-FUNCTION", "inact mut"),
                annotation);
    }
}