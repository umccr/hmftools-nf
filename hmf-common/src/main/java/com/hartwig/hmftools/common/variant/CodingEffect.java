package com.hartwig.hmftools.common.variant;

import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public enum CodingEffect {
    NONSENSE_OR_FRAMESHIFT,
    SPLICE,
    MISSENSE,
    NONE;

    @NotNull
    public static CodingEffect effect(@NotNull final List<VariantConsequence> consequences) {
        final List<CodingEffect> simplifiedEffects =
                consequences.stream().map(CodingEffect::effect).collect(Collectors.toList());
        if (simplifiedEffects.stream().anyMatch(x -> x.equals(NONSENSE_OR_FRAMESHIFT))) {
            return NONSENSE_OR_FRAMESHIFT;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(SPLICE))) {
            return SPLICE;
        }

        if (simplifiedEffects.stream().anyMatch(x -> x.equals(MISSENSE))) {
            return MISSENSE;
        }

        return NONE;
    }

    @NotNull
    private static CodingEffect effect(@NotNull final VariantConsequence consequence) {
        switch (consequence) {
            case FRAMESHIFT_VARIANT:
            case STOP_GAINED:
                return NONSENSE_OR_FRAMESHIFT;
            case MISSENSE_VARIANT:
            case PROTEIN_PROTEIN_CONTACT:
            case STRUCTURAL_INTERACTION_VARIANT:
            case INFRAME_DELETION:
            case INFRAME_INSERTION:
                return MISSENSE;
            case SPLICE_ACCEPTOR_VARIANT:
            case SPLICE_DONOR_VARIANT:
            case SPLICE_REGION_VARIANT:
                return SPLICE;
        }

        return NONE;
    }
}
