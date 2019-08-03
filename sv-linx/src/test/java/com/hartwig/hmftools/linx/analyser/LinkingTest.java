package com.hartwig.hmftools.linx.analyser;

import static Utils.SvTestRoutines.createDel;
import static Utils.SvTestRoutines.createDup;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_DB;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

public class LinkingTest
{
    @Test
    public void testLinkedPairs()
    {
        // test linked pair switching
        final SvVarData var1 = createDup(1, "1", 100, 200);
        final SvVarData var2 = createDup(2, "1", 300, 400);
        SvLinkedPair lp1 = new SvLinkedPair(var1, var2, LINK_TYPE_TI, false, true);

        assertEquals(lp1.first(), var1);
        assertEquals(lp1.second(), var2);
        assertEquals(lp1.firstLinkOnStart(), false);
        assertEquals(lp1.secondLinkOnStart(), true);
        assertEquals(lp1.firstUnlinkedOnStart(), true);
        assertEquals(lp1.secondUnlinkedOnStart(), false);

        lp1.switchSVs();
        assertEquals(lp1.first(), var2);
        assertEquals(lp1.second(), var1);
        assertEquals(lp1.firstLinkOnStart(), true);
        assertEquals(lp1.secondLinkOnStart(), false);

        // test short TIs converted to DBs
        final SvVarData var3 = createDel(3, "1", 100, 200);
        final SvVarData var4 = createDel(4, "1", 210, 400);
        SvLinkedPair lp2 = new SvLinkedPair(var3, var4, LINK_TYPE_TI, false, true);
        assertEquals(lp2.linkType(), LINK_TYPE_DB);
        assertEquals(lp2.length(), -11);

        final SvVarData var5 = createDel(4, "1", 250, 400);

        SvLinkedPair lp3 = new SvLinkedPair(var3, var5, LINK_TYPE_TI, false, true);

        lp3.sameVariants(lp2);
        lp3.hasLinkClash(lp2);
    }

}
