package com.hartwig.hmftools.serve.sources.ckb;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ActionableEntryFactoryTest {

    @Test
    public void testToActionableEntry() {
        CkbEntry entryDeletion =
                CkbTestFactory.createEntry("KRAS", "deletion", "KRAS deletion", "sensitive", "Emerging", "AB", "AB", "A", "DOID:162");
        Set<ActionableEntry> entryDeletionSet = ActionableEntryFactory.toActionableEntries(entryDeletion, "KRAS deletion");
        assertEquals(0, entryDeletionSet.size());

        CkbEntry entryCharacteristics =
                CkbTestFactory.createEntry("-", "MSI neg", "MSI neg", "sensitive", "Actionable", "AB", "AB", "A", "DOID:162");
        Set<ActionableEntry> entryCharacteristicsSet = ActionableEntryFactory.toActionableEntries(entryCharacteristics, "MSI neg");
        assertEquals(1, entryCharacteristicsSet.size());
        ActionableEntry characteristics = entryCharacteristicsSet.iterator().next();
        assertEquals("MSI neg", characteristics.source().sourceEvent());
        assertEquals(Knowledgebase.CKB, characteristics.source().source());
        assertEquals("AB", characteristics.treatment());
        assertEquals("AB", characteristics.cancerType());
        assertEquals("162", characteristics.doid());
        assertEquals("Hematologic cancer,2531", characteristics.tumorLocationBlacklisting());
        assertEquals(EvidenceLevel.A, characteristics.level());
        assertEquals(EvidenceDirection.RESPONSIVE, characteristics.direction());

        CkbEntry entryAmplification = CkbTestFactory.createEntry("KRAS",
                "amplification",
                "KRAS amplification",
                "sensitive",
                "Actionable",
                "AB",
                "AB",
                "A",
                "DOID:163");
        Set<ActionableEntry> entryAmplificationSet = ActionableEntryFactory.toActionableEntries(entryAmplification, "KRAS amplification");
        assertEquals(1, entryAmplificationSet.size());
        ActionableEntry amplification = entryAmplificationSet.iterator().next();
        assertEquals("KRAS amplification", amplification.source().sourceEvent());
        assertEquals(Knowledgebase.CKB, amplification.source().source());
        assertEquals("AB", amplification.treatment());
        assertEquals("AB", amplification.cancerType());
        assertEquals("163", amplification.doid());
        assertEquals(Strings.EMPTY, amplification.tumorLocationBlacklisting());
        assertEquals(EvidenceLevel.A, amplification.level());
        assertEquals(EvidenceDirection.RESPONSIVE, amplification.direction());

        CkbEntry entryHotspot =
                CkbTestFactory.createEntry("BRAF", "V600", "BRAF V600E", "sensitive", "Actionable", "AB", "AB", "A", "DOID:162");
        Set<ActionableEntry> entryHotspotSet = ActionableEntryFactory.toActionableEntries(entryHotspot, "BRAF V600E");
        assertEquals(1, entryHotspotSet.size());
        ActionableEntry hotspot = entryHotspotSet.iterator().next();
        assertEquals("BRAF V600E", hotspot.source().sourceEvent());
        assertEquals(Knowledgebase.CKB, hotspot.source().source());
        assertEquals("AB", hotspot.treatment());
        assertEquals("AB", hotspot.cancerType());
        assertEquals("162", hotspot.doid());
        assertEquals("Hematologic cancer,2531", hotspot.tumorLocationBlacklisting());
        assertEquals(EvidenceLevel.A, characteristics.level());
        assertEquals(EvidenceDirection.RESPONSIVE, characteristics.direction());
    }

    @Test
    public void canExtractAndMapDoid() {
        assertNull(ActionableEntryFactory.extractDoid(null));
        assertNull(ActionableEntryFactory.extractDoid("not a doid"));

        assertEquals("0060463", ActionableEntryFactory.extractDoid("DOID:0060463"));
        assertEquals("162", ActionableEntryFactory.extractDoid("JAX:10000003"));
        assertEquals("1749", ActionableEntryFactory.extractDoid("JAX:10000009"));
        assertEquals("299", ActionableEntryFactory.extractDoid("JAX:10000008"));
        assertNull(ActionableEntryFactory.extractDoid("JAX:10000004"));
    }

    @Test
    public void canExtractAndMapDoidKB() {
        assertNull(ActionableEntryFactory.extractDoidKB(null));
        assertNull(ActionableEntryFactory.extractDoidKB("not a doid"));

        assertEquals("0060463", ActionableEntryFactory.extractDoidKB("DOID:0060463"));
        assertEquals("10000003", ActionableEntryFactory.extractDoidKB("JAX:10000003"));
    }

    @Test
    public void canExtractResponseType() {
        assertEquals("predicted+-+sensitive", ActionableEntryFactory.extractResponseType("predicted - sensitive"));
        assertEquals("predicted+-+resistant", ActionableEntryFactory.extractResponseType("predicted - resistant"));
        assertEquals("resistant", ActionableEntryFactory.extractResponseType("resistant"));
        assertEquals("sensitive", ActionableEntryFactory.extractResponseType("sensitive"));
    }

    @Test
    public void canTestHasUsableEvidenceType() {
        assertTrue(ActionableEntryFactory.hasUsableEvidenceType("Actionable"));
        assertFalse(ActionableEntryFactory.hasUsableEvidenceType("Prognostic"));
        assertFalse(ActionableEntryFactory.hasUsableEvidenceType("Emerging"));
        assertFalse(ActionableEntryFactory.hasUsableEvidenceType("Risk Factor"));
        assertFalse(ActionableEntryFactory.hasUsableEvidenceType("Diagnostic"));
    }

    @Test
    public void canTestResolveLevel() {
        assertNull(ActionableEntryFactory.resolveLevel("NA"));
        assertEquals(EvidenceLevel.A, ActionableEntryFactory.resolveLevel("A"));
        assertEquals(EvidenceLevel.B, ActionableEntryFactory.resolveLevel("B"));
        assertEquals(EvidenceLevel.C, ActionableEntryFactory.resolveLevel("C"));
        assertEquals(EvidenceLevel.D, ActionableEntryFactory.resolveLevel("D"));
    }

    @Test
    public void canTestResolveDirection() {
        assertNull(ActionableEntryFactory.resolveDirection(null));
        assertNull(ActionableEntryFactory.resolveDirection("unknown"));
        assertNull(ActionableEntryFactory.resolveDirection("not applicable"));
        assertNull(ActionableEntryFactory.resolveDirection("conflicting"));
        assertNull(ActionableEntryFactory.resolveDirection("no benefit"));
        assertNull(ActionableEntryFactory.resolveDirection("not predictive"));
        assertNull(ActionableEntryFactory.resolveDirection("decreased response"));
        assertEquals(EvidenceDirection.RESPONSIVE, ActionableEntryFactory.resolveDirection("sensitive"));
        assertEquals(EvidenceDirection.PREDICTED_RESPONSIVE, ActionableEntryFactory.resolveDirection("predicted - sensitive"));
        assertEquals(EvidenceDirection.RESISTANT, ActionableEntryFactory.resolveDirection("resistant"));
        assertEquals(EvidenceDirection.PREDICTED_RESISTANT, ActionableEntryFactory.resolveDirection("predicted - resistant"));
    }
}