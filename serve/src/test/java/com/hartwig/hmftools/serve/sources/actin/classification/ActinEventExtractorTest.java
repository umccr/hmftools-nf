package com.hartwig.hmftools.serve.sources.actin.classification;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ActinEventExtractorTest {

    @Test
    public void canGenerateEventsForAllRules() {
        ImmutableActinEntry.Builder builder =
                ImmutableActinEntry.builder().trial(Strings.EMPTY).gene(Strings.EMPTY).mutation(Strings.EMPTY);

        for (ActinRule rule : ActinRule.values()) {
            assertNotNull(ActinEventExtractor.extractEvents(builder.rule(rule).build()));
        }
    }

    @Test
    public void canExtractEvent() {
        ActinEntry trial1 = ImmutableActinEntry.builder()
                .trial("trial1")
                .rule(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y)
                .gene("A")
                .mutation("mut")
                .build();

        assertEquals(Sets.newHashSet("mut"), ActinEventExtractor.extractEvents(trial1));

        ActinEntry trial2 = ImmutableActinEntry.builder()
                .trial("trial1")
                .rule(ActinRule.WILDTYPE_OF_GENE_X)
                .gene("A")
                .mutation("wildtype")
                .build();

        assertEquals(Sets.newHashSet("wildtype"), ActinEventExtractor.extractEvents(trial2));
    }
}