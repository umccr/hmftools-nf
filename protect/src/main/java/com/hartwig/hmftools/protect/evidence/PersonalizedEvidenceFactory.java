package com.hartwig.hmftools.protect.evidence;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ImmutableProtectSource;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.cancertype.CancerType;
import com.hartwig.hmftools.serve.cancertype.CancerTypeFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PersonalizedEvidenceFactory {

    private static final Logger LOGGER = LogManager.getLogger(PersonalizedEvidenceFactory.class);

    @NotNull
    private final Set<String> patientTumorDoids;

    public PersonalizedEvidenceFactory(@NotNull final Set<String> patientTumorDoids) {
        this.patientTumorDoids = patientTumorDoids;
    }

    @NotNull
    public ImmutableProtectEvidence.Builder somaticEvidence(@NotNull ActionableEvent event) {
        return evidenceBuilder(event).reported(false).germline(false);
    }

    @NotNull
    public ImmutableProtectEvidence.Builder somaticReportableEvidence(@NotNull ActionableEvent event) {
        return evidenceBuilder(event).reported(true).germline(false);
    }

    @NotNull
    public ImmutableProtectEvidence.Builder evidenceBuilder(@NotNull ActionableEvent actionable) {
        StringJoiner sourceUrlJoiner = new StringJoiner(",");
        for (String url : actionable.sourceUrls()) {
            sourceUrlJoiner.add(url);
        }

        return ImmutableProtectEvidence.builder()
                .treatment(actionable.treatment())
                .onLabel(determineOnLabel(actionable.applicableCancerType(), actionable.blacklistCancerTypes(), actionable.treatment()))
                .level(actionable.level())
                .direction(actionable.direction())
                .sources(determineProtectSources(actionable));
    }

    @NotNull
    public Set<ProtectSource> determineProtectSources(@NotNull ActionableEvent actionable) {
        Set<ProtectSource> sources = Sets.newHashSet();
        String sourceEvent = actionable.sourceEvent();
        Set<String> sourceUrls = actionable.sourceUrls();
        Integer rank = determineRangeRank(actionable);
        ProtectEvidenceType evidenceType = determineEvidenceType(actionable);
        Set<String> evidenceUrls = actionable.evidenceUrls();

        ProtectSource source = ImmutableProtectSource.builder()
                .name(actionable.source())
                .sourceEvent(sourceEvent)
                .sourceUrls(sourceUrls)
                .evidenceType(evidenceType)
                .rangeRank(rank)
                .evidenceUrls(evidenceUrls)
                .build();
        sources.add(source);
        return sources;
    }

    public boolean determineOnLabel(@NotNull CancerType applicableCancerType, @NotNull Set<CancerType> blacklistCancerTypes,
            @NotNull String treatment) {
        //TODO filter for blacklisting in v2.2. Should be analyzed in more depth
        return patientTumorDoids.contains(applicableCancerType.doid());
    }

    public boolean determineBlacklistedEvidence(@NotNull Set<CancerType> blacklistCancerTypes, @NotNull String treatment) {
        Set<String> blacklistDoids = CancerTypeFactory.doidStrings(blacklistCancerTypes);

        if (!blacklistDoids.isEmpty()) {
            LOGGER.info(" Starting doid resolving for blacklisting evidence  '{}' for treatment '{}'", blacklistDoids, treatment);
        }

        for (String result : blacklistDoids) {
            for (String doidPatient : patientTumorDoids) {
                if (doidPatient.equals(result)) {
                    return true;
                }
            }
        }
        return false;
    }

    @VisibleForTesting
    @NotNull
    static ProtectEvidenceType determineEvidenceType(@NotNull ActionableEvent actionable) {
        if (actionable instanceof ActionableHotspot) {
            return ProtectEvidenceType.HOTSPOT_MUTATION;
        } else if (actionable instanceof ActionableRange) {
            return fromActionableRange((ActionableRange) actionable);
        } else if (actionable instanceof ActionableGene) {
            return fromActionableGene((ActionableGene) actionable);
        } else if (actionable instanceof ActionableFusion) {
            return ProtectEvidenceType.FUSION_PAIR;
        } else if (actionable instanceof ActionableCharacteristic) {
            return fromActionableCharacteristic((ActionableCharacteristic) actionable);
        } else {
            throw new IllegalStateException("Unexpected actionable event detected in variant evidence: " + actionable);
        }
    }

    @NotNull
    private static ProtectEvidenceType fromActionableRange(@NotNull ActionableRange range) {
        switch (range.rangeType()) {
            case EXON:
                return ProtectEvidenceType.EXON_MUTATION;
            case CODON:
                return ProtectEvidenceType.CODON_MUTATION;
            default: {
                throw new IllegalStateException("Unsupported range type: " + range.rangeType());
            }
        }
    }

    @NotNull
    private static ProtectEvidenceType fromActionableGene(@NotNull ActionableGene gene) {
        switch (gene.event()) {
            case AMPLIFICATION:
                return ProtectEvidenceType.AMPLIFICATION;
            case DELETION:
                return ProtectEvidenceType.DELETION;
            case ACTIVATION:
                return ProtectEvidenceType.ACTIVATION;
            case INACTIVATION:
                return ProtectEvidenceType.INACTIVATION;
            case ANY_MUTATION:
                return ProtectEvidenceType.ANY_MUTATION;
            case FUSION:
                return ProtectEvidenceType.PROMISCUOUS_FUSION;
            case WILD_TYPE:
                return ProtectEvidenceType.WILD_TYPE;
            default: {
                throw new IllegalStateException("Unsupported gene level event: " + gene.event());
            }
        }
    }

    @NotNull
    private static ProtectEvidenceType fromActionableCharacteristic(@NotNull ActionableCharacteristic characteristic) {
        switch (characteristic.name()) {
            case MICROSATELLITE_UNSTABLE:
            case MICROSATELLITE_STABLE:
            case HIGH_TUMOR_MUTATIONAL_LOAD:
            case LOW_TUMOR_MUTATIONAL_LOAD:
            case HIGH_TUMOR_MUTATIONAL_BURDEN:
            case LOW_TUMOR_MUTATIONAL_BURDEN:
            case HOMOLOGOUS_RECOMBINATION_DEFICIENT:
                return ProtectEvidenceType.SIGNATURE;
            case HPV_POSITIVE:
            case EBV_POSITIVE:
                return ProtectEvidenceType.VIRAL_PRESENCE;
            default: {
                throw new IllegalStateException("Unsupported tumor characteristic: " + characteristic.name());
            }
        }
    }

    @VisibleForTesting
    @Nullable
    static Integer determineRangeRank(@NotNull ActionableEvent actionable) {
        return actionable instanceof ActionableRange ? ((ActionableRange) actionable).rank() : null;
    }
}