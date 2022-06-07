package com.hartwig.hmftools.sage.common;

import static java.util.Arrays.fill;

import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class IndexedBasesTest
{
    private final IndexedBases mIndexedBases = new IndexedBases(
            1000, 5, 4, 6, 3, "GATCTCCTCA".getBytes());

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
        assertEquals(-1, mIndexedBases.rightFlankMatchingBases(3, "TCTCCTCG".getBytes()));

        assertEquals(3, mIndexedBases.rightFlankMatchingBases(3, "TCTCCTCAG".getBytes()));
        assertEquals(3, mIndexedBases.rightFlankMatchingBases(3, "TCTCCTCA".getBytes()));
        assertEquals(2, mIndexedBases.rightFlankMatchingBases(3, "TCTCCTC".getBytes()));
        assertEquals(1, mIndexedBases.rightFlankMatchingBases(3, "TCTCCT".getBytes()));
        assertEquals(0, mIndexedBases.rightFlankMatchingBases(3, "TCTCC".getBytes()));
    }

    @Test
    public void testLeftFlankMatchingBases()
    {
        assertEquals(-1, mIndexedBases.leftFlankMatchingBases(5, "TTCTCCTCA".getBytes()));

        assertEquals(3, mIndexedBases.leftFlankMatchingBases(5, "GATCTCCTCA".getBytes()));
        assertEquals(3, mIndexedBases.leftFlankMatchingBases(4, "ATCTCCTCA".getBytes()));
        assertEquals(2, mIndexedBases.leftFlankMatchingBases(3, "TCTCCTCA".getBytes()));
        assertEquals(1, mIndexedBases.leftFlankMatchingBases(2, "CTCCTCA".getBytes()));
        assertEquals(0, mIndexedBases.leftFlankMatchingBases(1, "TCCTCA".getBytes()));
    }

    @Test
    public void testCoreMatch()
    {
        assertTrue(mIndexedBases.coreMatch(true, 5, "GATCT.CTCA".getBytes()));
        assertFalse(mIndexedBases.coreMatch(false, 5, "GATCT.CTCA".getBytes()));

        assertTrue(mIndexedBases.coreMatch(false, 5, "GATCTCCTCA".getBytes()));
        assertTrue(mIndexedBases.coreMatch(false, 1, "TCC".getBytes()));

        assertFalse(mIndexedBases.coreMatch(false, 1, "CCC".getBytes()));
        assertFalse(mIndexedBases.coreMatch(false, 1, "TTC".getBytes()));
        assertFalse(mIndexedBases.coreMatch(false, 1, "TCT".getBytes()));
        assertFalse(mIndexedBases.coreMatch(false, 1, "TC".getBytes()));
        assertFalse(mIndexedBases.coreMatch(false, 0, "CC".getBytes()));
    }

    @Test
    public void testPartialMatchMustHaveAtLeastOneFullSide()
    {
        IndexedBases indexedBases = new IndexedBases(1000, 2, 2, 2, 2, "GGTAA".getBytes());

        IndexedBases testBases = new IndexedBases(1000, 2, 2, 2, 2, "GGTAA".getBytes());
        assertEquals(ReadContextMatch.FULL, indexedBases.matchAtPosition(testBases));

        testBases = new IndexedBases(1000, 2, 2, 2, 2, "GGTA".getBytes());
        assertEquals(ReadContextMatch.PARTIAL, indexedBases.matchAtPosition(testBases));

        testBases = new IndexedBases(1000, 2, 2, 2, 2, "GGT".getBytes());
        assertEquals(ReadContextMatch.PARTIAL, indexedBases.matchAtPosition(testBases));

        testBases = new IndexedBases(1000, 1, 2, 2, 2, "GT".getBytes());
        assertEquals(ReadContextMatch.CORE, indexedBases.matchAtPosition(testBases));

        testBases = new IndexedBases(1000, 1, 2, 2, 2, "GTAA".getBytes());
        assertEquals(ReadContextMatch.PARTIAL, indexedBases.matchAtPosition(testBases));

        testBases = new IndexedBases(1000, 0, 2, 2, 2, "TAA".getBytes());
        assertEquals(ReadContextMatch.PARTIAL, indexedBases.matchAtPosition(testBases));

        testBases = new IndexedBases(1000, 0, 2, 2, 2, "TA".getBytes());
        assertEquals(ReadContextMatch.CORE, indexedBases.matchAtPosition(testBases));

        testBases = new IndexedBases(1000, 0, 2, 2, 2, "T".getBytes());
        assertEquals(ReadContextMatch.CORE, indexedBases.matchAtPosition(testBases));
    }

    @Test
    public void testNegativeReadIndex()
    {
        IndexedBases indexedBases = new IndexedBases(1000, 2, 2, 2, 2, "GGTAA".getBytes());

        IndexedBases testBases = new IndexedBases(1000, 2, 2, 2, 2, "GGTAA".getBytes());
        assertEquals(ReadContextMatch.FULL, indexedBases.matchAtPosition(testBases));

        testBases = new IndexedBases(1000, -1, 2, 2, 2, "GGTAA".getBytes());
        assertEquals(ReadContextMatch.NONE, indexedBases.matchAtPosition(testBases));
    }

    @Test
    public void testBaseQualityMatch()
    {
        String preLeftFlank = "TTTTTTTTTT";
        String leftFlank = "GGGGG";
        String rightFlank = "AAAAA";
        String core = "AATCC";
        String refBases = preLeftFlank + leftFlank + core + rightFlank;

        int position = 1000;

        int flankSize = leftFlank.length();
        int leftFlankIndex = 10;
        int leftCoreIndex = leftFlankIndex + flankSize;
        int readIndex = leftCoreIndex + 2;
        int rightCoreIndex = readIndex + 2;

        IndexedBases indexedBases = new IndexedBases(position, readIndex, leftCoreIndex, rightCoreIndex, flankSize, refBases.getBytes());

        byte[] readQuals = new byte[refBases.length()];
        byte highQual = MATCHING_BASE_QUALITY + 1;
        fill(readQuals, highQual);

        IndexedBases testBases = new IndexedBases(position, readIndex, leftCoreIndex, rightCoreIndex, flankSize, refBases.getBytes());
        assertEquals(ReadContextMatch.FULL, indexedBases.matchAtPosition(testBases, false, readQuals));

        // base diff in core - permitted if not in an actual SNV or MNV
        String readBases = preLeftFlank + leftFlank + "ATTCC" + rightFlank;
        testBases = new IndexedBases(position, readIndex, leftCoreIndex, rightCoreIndex, flankSize, readBases.getBytes());
        assertEquals(ReadContextMatch.NONE, indexedBases.matchAtPosition(testBases, false, readQuals));

        // with low qual at that base
        readQuals[leftCoreIndex + 1] = MATCHING_BASE_QUALITY - 1;
        //assertEquals(ReadContextMatch.FULL, indexedBases.matchAtPosition(testBases, false, readQuals));

        // not permitted in the alt itself
        readBases = preLeftFlank + leftFlank + "AAGCC" + rightFlank;
        testBases = new IndexedBases(position, readIndex, leftCoreIndex, rightCoreIndex, flankSize, readBases.getBytes());
        assertEquals(ReadContextMatch.NONE, indexedBases.matchAtPosition(testBases, false, readQuals));

        // with low qual at that base
        readQuals[leftCoreIndex + 2] = MATCHING_BASE_QUALITY - 1;
        assertEquals(ReadContextMatch.NONE, indexedBases.matchAtPosition(testBases, false, readQuals));

        // low qual mismatch in left flank
        readBases = preLeftFlank + "GGGAG" + core + rightFlank;
        testBases = new IndexedBases(position, readIndex, leftCoreIndex, rightCoreIndex, flankSize, readBases.getBytes());
        assertEquals(ReadContextMatch.CORE, indexedBases.matchAtPosition(testBases, false, readQuals));

        readQuals[leftFlankIndex + 3] = MATCHING_BASE_QUALITY - 1;
        assertEquals(ReadContextMatch.FULL, indexedBases.matchAtPosition(testBases, false, readQuals));

        // low qual mismstch in right flank
        readBases = preLeftFlank + leftFlank + core + "GAAAA";
        testBases = new IndexedBases(position, readIndex, leftCoreIndex, rightCoreIndex, flankSize, readBases.getBytes());
        assertEquals(ReadContextMatch.CORE, indexedBases.matchAtPosition(testBases, false, readQuals));

        readQuals[rightCoreIndex + 1] = MATCHING_BASE_QUALITY - 1;
        assertEquals(ReadContextMatch.FULL, indexedBases.matchAtPosition(testBases, false, readQuals));

        readBases = preLeftFlank + leftFlank + core + "AAAAG";
        testBases = new IndexedBases(position, readIndex, leftCoreIndex, rightCoreIndex, flankSize, readBases.getBytes());
        assertEquals(ReadContextMatch.CORE, indexedBases.matchAtPosition(testBases, false, readQuals));

        readQuals[rightCoreIndex + 5] = MATCHING_BASE_QUALITY - 1;
        assertEquals(ReadContextMatch.FULL, indexedBases.matchAtPosition(testBases, false, readQuals));
    }

    @Test
    public void testStrings()
    {
        ReadContext victim = createReadContext(1000, 4, 3, 5, 2, "AACATGAGG", Strings.EMPTY);
        assertEquals("ATG", victim.coreString());
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
        // int[] baseQualities = makeDefaultBaseQualitities(readBases.length());

        return new ReadContext(refPosition, "", 0, microhomology, readBasesIndexed, incompleteCore);
    }


    @Test
    public void testCreate()
    {
        IndexedBases victimWithExtra = createIndexedBases(1000, 1, "AA", "TA", "ATG", "CG", "TT");
        assertEquals(1, victimWithExtra.indexInCore());
        assertEquals(5, victimWithExtra.Index);
        assertEquals("TA", victimWithExtra.leftFlankString());
        assertEquals("ATG", victimWithExtra.coreString());
        assertEquals("CG", victimWithExtra.rightFlankString());

        IndexedBases victimWithoutExtra = createIndexedBases(1000, 1, Strings.EMPTY, "TA", "ATG", "CG", Strings.EMPTY);
        assertEquals(1, victimWithoutExtra.indexInCore());
        assertEquals("TA", victimWithoutExtra.leftFlankString());
        assertEquals("ATG", victimWithoutExtra.coreString());
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
