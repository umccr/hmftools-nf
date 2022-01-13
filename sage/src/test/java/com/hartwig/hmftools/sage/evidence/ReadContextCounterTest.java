package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.read.ReadContextTest.makeDefaultBaseQualitities;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.variant.VariantTier;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextCounterTest
{
    private static final String SAMPLE = "sample";
    private static final int MAX_COVERAGE = 1000;
    private static final VariantTier TIER = VariantTier.PANEL;
    private static final SageConfig CONFIG = new SageConfig();
    private static final QualityRecalibrationMap RECALIBRATION = new QualityRecalibrationMap(Collections.emptyList());

    private static final Map<String,QualityRecalibrationMap> RECALIBRATION_MAP = Maps.newHashMap();

    // convert to using MockRefGenome
    private static final IndexedBases REF_BASES = new IndexedBases(550, 0, "TGTTTCTGTTTC".getBytes());

    static
    {
        RECALIBRATION_MAP.put(SAMPLE, RECALIBRATION);
    }

    private static final QualityCalculator QUALITY_CALCULATOR = new QualityCalculator(CONFIG.Quality, RECALIBRATION_MAP, REF_BASES);

    @Test
    public void testInsertInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("GT").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 5, "TGTTTC", Strings.EMPTY);

        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(555, "3S3M", "TGTTTC", "######");
        victim.accept(record, CONFIG, QUALITY_CALCULATOR,1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testDeleteInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("GT").alt("G").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 4, "TGTTC", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(556, "2S3M", "TGTTC", "#####");
        victim.accept(record, CONFIG, QUALITY_CALCULATOR,1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testSnvInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("A").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 2, "CAT", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(555, "2S1M", "CAT", "#####");
        victim.accept(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testRefInLeftSoftClipDoesNotContributeToDepth()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("A").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 2,"CAT", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(555, "2S1M", "CGT", "#####");
        victim.accept(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(0, victim.depth());
        assertEquals(0, victim.altSupport());
    }

    @Test
    public void testMnvInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("TCG").alt("ATC").position(552).build();
        final ReadContext readContext = createReadContext(552, 2, 0, 6, "GAAAAAT", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(555, "5S3M", "GAAAAATC", "########");
        victim.accept(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testInsertInRightSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("GT").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 5, "TGTTTC", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(553, "2M4S", "TGTTTC", "######");
        victim.accept(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testDeleteInRightSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("GT").alt("G").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 4, "TGTTC", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(553, "2M3S", "TGTTC", "#####");
        victim.accept(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testMnvInRightSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("TCG").alt("ATC").position(552).build();

        final ReadContext readContext = createReadContext(552, 2, 0, 6, "GAAAAAT", Strings.EMPTY);

        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(550, "2M6S", "GAAAAATC", "########");
        victim.accept(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    public static ReadContext createReadContext(
            int refPosition, int readIndex, int leftCentreIndex, int rightCentreIndex, String readBases, String microhomology)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readBases.length() - 1);
        boolean incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        IndexedBases readBasesIndexed = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, 0, readBases.getBytes());
        int[] baseQualities = makeDefaultBaseQualitities(readBases.length());

        return new ReadContext(refPosition, "", 0, microhomology, readBasesIndexed, baseQualities, incompleteCore);
    }

    public static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString,
            @NotNull final String qualities)
    {
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);
        return record;
    }
}