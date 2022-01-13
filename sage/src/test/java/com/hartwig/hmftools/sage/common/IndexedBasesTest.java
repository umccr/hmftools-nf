package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.sage.read.ReadContextTest.makeDefaultBaseQualitities;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextMatch;

import org.apache.logging.log4j.util.Strings;
import org.junit.Assert;
import org.junit.Test;

public class IndexedBasesTest
{

    private final IndexedBases victim = new IndexedBases(1000, 5, 4, 6, 3, "GATCTCCTCA".getBytes());

    @Test
    public void testMaxFlankLength()
    {
        IndexedBases constrainedOnLeft = new IndexedBases(1000, 5, 2, 6, 3, "GATCTCCTCA".getBytes());
        IndexedBases constrainedOnRight = new IndexedBases(1000, 5, 4, 6, 3, "GATCTCCTCA".getBytes());

        assertEquals(2, constrainedOnLeft.maxFlankLength());
        assertEquals(3, constrainedOnRight.maxFlankLength());
    }

    @Test
    public void testRightFlankMatchingBases()
    {
        assertEquals(-1, victim.rightFlankMatchingBases(3, "TCTCCTCG".getBytes()));

        assertEquals(3, victim.rightFlankMatchingBases(3, "TCTCCTCAG".getBytes()));
        assertEquals(3, victim.rightFlankMatchingBases(3, "TCTCCTCA".getBytes()));
        assertEquals(2, victim.rightFlankMatchingBases(3, "TCTCCTC".getBytes()));
        assertEquals(1, victim.rightFlankMatchingBases(3, "TCTCCT".getBytes()));
        assertEquals(0, victim.rightFlankMatchingBases(3, "TCTCC".getBytes()));
    }

    @Test
    public void testLeftFlankMatchingBases()
    {
        assertEquals(-1, victim.leftFlankMatchingBases(5, "TTCTCCTCA".getBytes()));

        assertEquals(3, victim.leftFlankMatchingBases(5, "GATCTCCTCA".getBytes()));
        assertEquals(3, victim.leftFlankMatchingBases(4, "ATCTCCTCA".getBytes()));
        assertEquals(2, victim.leftFlankMatchingBases(3, "TCTCCTCA".getBytes()));
        assertEquals(1, victim.leftFlankMatchingBases(2, "CTCCTCA".getBytes()));
        assertEquals(0, victim.leftFlankMatchingBases(1, "TCCTCA".getBytes()));
    }

    @Test
    public void testCoreMatch()
    {
        assertTrue(victim.coreMatch(true, 5, "GATCT.CTCA".getBytes()));
        assertFalse(victim.coreMatch(false, 5, "GATCT.CTCA".getBytes()));

        assertTrue(victim.coreMatch(false, 5, "GATCTCCTCA".getBytes()));
        assertTrue(victim.coreMatch(false, 1, "TCC".getBytes()));

        assertFalse(victim.coreMatch(false, 1, "CCC".getBytes()));
        assertFalse(victim.coreMatch(false, 1, "TTC".getBytes()));
        assertFalse(victim.coreMatch(false, 1, "TCT".getBytes()));
        assertFalse(victim.coreMatch(false, 1, "TC".getBytes()));
        assertFalse(victim.coreMatch(false, 0, "CC".getBytes()));
    }

    @Test
    public void testPartialMatchMustHaveAtLeastOneFullSide()
    {
        ReadContext victim = createReadContext(1000, 2, 2, 2, 2, "GGTAA", Strings.EMPTY);
        Assert.assertEquals(ReadContextMatch.FULL, victim.matchAtPosition(false, 2, "GGTAA".getBytes()));

        assertEquals(ReadContextMatch.PARTIAL, victim.matchAtPosition(false, 2, "GGTA".getBytes()));
        assertEquals(ReadContextMatch.PARTIAL, victim.matchAtPosition(false, 2, "GGT".getBytes()));
        assertEquals(ReadContextMatch.CORE, victim.matchAtPosition(false, 1, "GT".getBytes()));

        assertEquals(ReadContextMatch.PARTIAL, victim.matchAtPosition(false, 1, "GTAA".getBytes()));
        assertEquals(ReadContextMatch.PARTIAL, victim.matchAtPosition(false, 0, "TAA".getBytes()));
        assertEquals(ReadContextMatch.CORE, victim.matchAtPosition(false, 0, "TA".getBytes()));
        assertEquals(ReadContextMatch.CORE, victim.matchAtPosition(false, 0, "T".getBytes()));
    }

    @Test
    public void testNegativeReadIndex()
    {
        ReadContext victim = createReadContext(1000, 2, 2, 2, 2, "GGTAA", Strings.EMPTY);
        assertEquals(ReadContextMatch.FULL, victim.matchAtPosition(false, 2, "GGTAA".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.matchAtPosition(false, -1, "GGTAA".getBytes()));
    }

    @Test
    public void testPhasedMNV()
    {
        ReadContext victim1 = createReadContext(1000, 4, 4, 4, 4, "GATCTTGAT", Strings.EMPTY);
        ReadContext victim2 = createReadContext(1001, 5, 5, 5, 4, "GATCTTGATC", Strings.EMPTY);

        assertTrue(victim1.phased(-1, victim2));
        assertTrue(victim2.phased(1, victim1));
    }

    @Test
    public void testPhasedReadLongEnoughOnAtLeastOneSide()
    {
        ReadContext victim1 = createReadContext(1000, 4, 4, 4, 4, "GATCTTGA", Strings.EMPTY);
        ReadContext victim2 = createReadContext(1001, 5, 5, 5, 4, "GATCTTGATCT", Strings.EMPTY);

        assertTrue(victim1.phased(-1, victim2));
        assertTrue(victim2.phased(1, victim1));
    }

    @Test
    public void testLongReadShortFlanks()
    {
        final String read = "TACCACAAATACATATACGTGTATCTGTCTGTGTGTTATGAACTTATATAAACCATCAC";
        ReadContext victim1 = createReadContext(1010, 10, 9, 11, 3, read, Strings.EMPTY);
        ReadContext victim2 = createReadContext(1030, 40, 39, 41, 3, read, Strings.EMPTY);

        assertTrue(victim1.phased(-30, victim2));
        assertTrue(victim2.phased(30, victim1));
    }

    @Test
    public void testBothCentreMatches()
    {
        ReadContext victim1 = createReadContext(1000, 4, 4, 4, 4, "AAAATGGGG", Strings.EMPTY);
        ReadContext victim2 = createReadContext(1005, 5, 5, 5, 4, "TGGGGACCCC", Strings.EMPTY);
        assertFalse(victim1.phased(-5, victim2));
        assertFalse(victim2.phased(5, victim1));
    }

    @Test
    public void testStrings()
    {
        ReadContext victim = createReadContext(1000, 4, 3, 5, 2, "AACATGAGG", Strings.EMPTY);
        assertEquals("ATG", victim.centerBases());
        assertEquals("AC", victim.leftFlankString());
        assertEquals("AG", victim.rightFlankString());
    }

    public static ReadContext createReadContext(
            int refPosition, int readIndex, int leftCentreIndex, int rightCentreIndex, int flankSize, String readBases, String microhomology)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readBases.length() - 1);
        boolean incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        IndexedBases readBasesIndexed = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, flankSize, readBases.getBytes());
        int[] baseQualities = makeDefaultBaseQualitities(readBases.length());

        return new ReadContext(refPosition, "", 0, microhomology, readBasesIndexed, baseQualities, incompleteCore);
    }


    @Test
    public void testCreate()
    {
        IndexedBases victimWithExtra = createIndexedBases(1000, 1, "AA", "TA", "ATG", "CG", "TT");
        assertEquals(1, victimWithExtra.indexInCore());
        assertEquals(5, victimWithExtra.Index);
        assertEquals("TA", victimWithExtra.leftFlankString());
        assertEquals("ATG", victimWithExtra.centerString());
        assertEquals("CG", victimWithExtra.rightFlankString());

        IndexedBases victimWithoutExtra = createIndexedBases(1000, 1, Strings.EMPTY, "TA", "ATG", "CG", Strings.EMPTY);
        assertEquals(1, victimWithoutExtra.indexInCore());
        assertEquals("TA", victimWithoutExtra.leftFlankString());
        assertEquals("ATG", victimWithoutExtra.centerString());
        assertEquals("CG", victimWithoutExtra.rightFlankString());
    }

    public static IndexedBases createIndexedBases(
            int position, int indexInCore, String leftExtra, String leftFlank, String core, String rightFlank, String rightExtra)
    {
        assertTrue(indexInCore <= core.length());

        int flankSize = Math.max(leftFlank.length(), rightFlank.length());
        int totalLength = leftExtra.length() + leftFlank.length() + core.length() + rightFlank.length() + rightExtra.length();
        byte[] bases = new byte[totalLength];

        int destPos = 0;
        System.arraycopy(leftExtra.getBytes(), 0, bases, destPos, leftExtra.length());
        System.arraycopy(leftFlank.getBytes(), 0, bases, destPos += leftExtra.length(), leftFlank.length());
        System.arraycopy(core.getBytes(), 0, bases, destPos += leftFlank.length(), core.length());
        System.arraycopy(rightFlank.getBytes(), 0, bases, destPos += core.length(), rightFlank.length());
        System.arraycopy(rightExtra.getBytes(), 0, bases, destPos += rightFlank.length(), rightExtra.length());

        int leftCoreIndex = leftExtra.length() + leftFlank.length();
        int rightCoreIndex = leftCoreIndex + core.length() - 1;
        return new IndexedBases(position,
                leftExtra.length() + leftFlank.length() + indexInCore,
                leftCoreIndex,
                rightCoreIndex,
                flankSize,
                bases);
    }

}