package com.hartwig.hmftools.serve.extraction.characteristic;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.immuno.ImmunoHLA;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TumorCharacteristicExtractorTest {

    private static final double EPSILON = 1e-10;

    private static final String MSI = "msi";
    private static final String MSS = "mss";
    private static final String HIGH_TML = "high_tml";
    private static final String LOW_TML = "low_tml";
    private static final String HIGH_TMB = "high_tmb";
    private static final String LOW_TMB = "low_tmb";
    private static final String HRD = "hrd";
    private static final String HPV = "hpv";
    private static final String EBV = "ebv";
    private static final String IMMUNO_HLA = "hla";

    @Test
    public void canDetermineCutOffMSI() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSI, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, characteristic.tumorCharacteristicAnnotation());

        assertEquals(TumorCharacteristicsAtLeast.EQUALS_GREATHER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, "MSI >= 4"));
        assertEquals(4,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, "MSI >= 4"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffMSS() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSS, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, characteristic.tumorCharacteristicAnnotation());

        assertEquals(TumorCharacteristicsAtLeast.LESSER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, "MSS < 4"));
        assertEquals(4,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, "MSS < 4"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffTMLLow() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, LOW_TML, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());

        assertEquals(TumorCharacteristicsAtLeast.LESSER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, "TML < 140"));
        assertEquals(140,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, "TML < 140"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffTMLHigh() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TML, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());

        assertEquals(TumorCharacteristicsAtLeast.EQUALS_GREATHER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, "TML >= 140"));
        assertEquals(140,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, "TML >= 140"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffTMBLow() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, LOW_TMB, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN, characteristic.tumorCharacteristicAnnotation());

        assertEquals(TumorCharacteristicsAtLeast.LESSER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN, "TMB < 3"));
        assertEquals(3,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_BURDEN, "TMB < 3"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffTMBHigh() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TMB, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN, characteristic.tumorCharacteristicAnnotation());

        assertEquals(TumorCharacteristicsAtLeast.EQUALS_GREATHER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN, "TMB >= 14"));
        assertEquals(14,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_BURDEN, "TMB >= 14"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffHRD() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HRD, Strings.EMPTY);
        assertEquals(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT, characteristic.tumorCharacteristicAnnotation());

        assertEquals(TumorCharacteristicsAtLeast.EQUALS_GREATHER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT,
                        "HRD >= 0.5"));
        assertEquals(0.5,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT,
                        "HRD >= 0.5"),
                EPSILON);
    }

    @Test
    public void canDetermineCutOffSource() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TML, "TML >= 200");
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());

        assertEquals(TumorCharacteristicsAtLeast.EQUALS_GREATHER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, "TML >= 200"));
        assertEquals(200.0,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, "TML >= 200"),
                EPSILON);

    }

    @Test
    public void canExtractMicrosatelliteUnstableCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSI, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, characteristic.tumorCharacteristicAnnotation());
        assertEquals(TumorCharacteristicsAtLeast.EQUALS_GREATHER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, Strings.EMPTY));
        assertEquals(4.0,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE, Strings.EMPTY),
                EPSILON);
    }

    @Test
    public void canExtractMicrosatelliteStableCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, MSS, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, characteristic.tumorCharacteristicAnnotation());
        assertEquals(TumorCharacteristicsAtLeast.LESSER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, Strings.EMPTY));
        assertEquals(4.0,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.MICROSATELLITE_STABLE, Strings.EMPTY),
                EPSILON);
    }

    @Test
    public void canExtractHighTumorMutationalLoadCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HIGH_TML, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());
        assertEquals(TumorCharacteristicsAtLeast.EQUALS_GREATHER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, Strings.EMPTY));
        assertEquals(140,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HIGH_TUMOR_MUTATIONAL_LOAD, Strings.EMPTY),
                EPSILON);
    }

    @Test
    public void canExtractLowTumorMutationalLoadCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, LOW_TML, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, characteristic.tumorCharacteristicAnnotation());
        assertEquals(TumorCharacteristicsAtLeast.LESSER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, Strings.EMPTY));
        assertEquals(140,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.LOW_TUMOR_MUTATIONAL_LOAD, Strings.EMPTY),
                EPSILON);
    }

    @Test
    public void canExtractHrDeficientCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HRD, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT, characteristic.tumorCharacteristicAnnotation());
        assertEquals(TumorCharacteristicsAtLeast.EQUALS_GREATHER,
                tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT,
                        Strings.EMPTY));
        assertEquals(0.5,
                tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HOMOLOGOUS_RECOMBINATION_DEFICIENT,
                        Strings.EMPTY),
                EPSILON);
    }

    @Test
    public void canExtractHPVPositiveCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, HPV, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.HPV_POSITIVE, characteristic.tumorCharacteristicAnnotation());
        assertNull(tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.HPV_POSITIVE, Strings.EMPTY));
        assertNull(tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HPV_POSITIVE, Strings.EMPTY));
    }

    @Test
    public void canExtractEBVPositiveCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();
        TumorCharacteristic characteristic = tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, EBV, Strings.EMPTY);

        assertNotNull(characteristic);
        assertEquals(TumorCharacteristicAnnotation.EBV_POSITIVE, characteristic.tumorCharacteristicAnnotation());
        assertNull(tumorCharacteristicExtractor.determineAtLeast(TumorCharacteristicAnnotation.HPV_POSITIVE, Strings.EMPTY));
        assertNull(tumorCharacteristicExtractor.determineCutoff(TumorCharacteristicAnnotation.HPV_POSITIVE, Strings.EMPTY));
    }

    @Test
    public void canFilterUnknownCharacteristic() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();

        assertNull(tumorCharacteristicExtractor.extract(EventType.CHARACTERISTIC, "Not a tumor characteristic", Strings.EMPTY));
    }

    @Test
    public void canFilterWrongTypes() {
        TumorCharacteristicExtractor tumorCharacteristicExtractor = buildTestExtractor();

        assertNull(tumorCharacteristicExtractor.extract(EventType.COMPLEX, MSI, Strings.EMPTY));
    }

    @NotNull
    private static TumorCharacteristicExtractor buildTestExtractor() {
        return new TumorCharacteristicExtractor(Sets.newHashSet(MSI),
                Sets.newHashSet(MSS),
                Sets.newHashSet(HIGH_TML),
                Sets.newHashSet(LOW_TML),
                Sets.newHashSet(HIGH_TMB),
                Sets.newHashSet(LOW_TMB),
                Sets.newHashSet(HRD),
                Sets.newHashSet(HPV),
                Sets.newHashSet(EBV));
    }
}