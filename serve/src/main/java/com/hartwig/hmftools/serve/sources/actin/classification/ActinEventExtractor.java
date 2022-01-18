package com.hartwig.hmftools.serve.sources.actin.classification;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.jetbrains.annotations.NotNull;

public final class ActinEventExtractor {

    private ActinEventExtractor() {
    }

    @NotNull
    public static Set<String> extractEvents(@NotNull ActinEntry entry) {
        switch (entry.rule()) {
            case ACTIVATION_OR_AMPLIFICATION_OF_GENE_X:
                return Sets.newHashSet(ActinKeywords.ACTIVATION, ActinKeywords.AMPLIFICATION);
            case ACTIVATING_MUTATION_IN_GENE_X:
                return Sets.newHashSet(ActinKeywords.ACTIVATION);
            case INACTIVATION_OF_GENE_X:
                return Sets.newHashSet(ActinKeywords.INACTIVATION);
            case MUTATION_IN_GENE_X_OF_TYPE_Y: {
                String mutation = entry.mutation();
                if (mutation == null) {
                    throw new IllegalStateException("No mutation provided in ACTIN entry: " + entry);
                }
                return Sets.newHashSet(mutation);
            }
            case AMPLIFICATION_OF_GENE_X:
                return Sets.newHashSet(ActinKeywords.AMPLIFICATION);
            case DELETION_OF_GENE_X:
                return Sets.newHashSet(ActinKeywords.DELETION);
            case ACTIVATING_FUSION_IN_GENE_X:
                return Sets.newHashSet(ActinKeywords.PROMISCUOUS_FUSION);
            case SPECIFIC_FUSION_X:
                return Sets.newHashSet(entry.gene() + " fusion");
            default: {
                throw new IllegalStateException("Unrecognized event: " + entry.rule());
            }
        }
    }
}