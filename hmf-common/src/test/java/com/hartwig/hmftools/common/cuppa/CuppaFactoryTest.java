package com.hartwig.hmftools.common.cuppa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CuppaFactoryTest {

    private static final String CUPPA_DATA_CSV = Resources.getResource("cuppa/sample.cup.data.csv").getPath();

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canInterpretFromTestFile() throws IOException {
        List<CuppaEntry> entries = CuppaDataFile.read(CUPPA_DATA_CSV);

        CuppaData cuppa = CuppaFactory.create(entries);

        List<CuppaPrediction> predictions = cuppa.predictions();
        assertEquals(1, predictions.size());
        assertEquals("Esophagus/Stomach", predictions.get(0).cancerType());
        assertEquals(0.000132, predictions.get(0).likelihood(), EPSILON);

        assertEquals(10, cuppa.simpleDups32To200B());
        assertEquals(8, cuppa.maxComplexSize());
        assertEquals(5, cuppa.LINECount());
        assertEquals(0, cuppa.telomericSGLs());
    }

    @Test
    public void respectOrderingOfCombinedDataTypes() {
        CuppaEntry rna = combined().dataType(CuppaFactory.RNA_COMBINED_DATATYPE).refCancerType("rna").build();
        CuppaEntry dna = combined().dataType(CuppaFactory.DNA_COMBINED_DATATYPE).refCancerType("dna").build();
        CuppaEntry overall = combined().dataType(CuppaFactory.OVERALL_COMBINED_DATATYPE).refCancerType("overall").build();

        CuppaData fromRna = CuppaFactory.create(Lists.newArrayList(rna));
        assertEquals("rna", fromRna.predictions().get(0).cancerType());

        CuppaData fromRnaDna = CuppaFactory.create(Lists.newArrayList(rna, dna));
        assertEquals("dna", fromRnaDna.predictions().get(0).cancerType());

        CuppaData fromAll = CuppaFactory.create(Lists.newArrayList(rna, dna, overall));
        assertEquals("overall", fromAll.predictions().get(0).cancerType());
    }

    @Test
    public void doNotCrashOnMissingEntries() {
        assertNotNull(CuppaFactory.create(Lists.newArrayList()));
    }

    @NotNull
    private static ImmutableCuppaEntry.Builder combined() {
        return CuppaTestFactory.builder().category(CuppaFactory.COMBINED_CATEGORY);
    }
}