package com.hartwig.hmftools.common.variant.impact;

import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffUtils.SNPEFF_CANONICAL;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffUtils.SNPEFF_WORST;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.logging.log4j.util.Strings;

import htsjdk.variant.variantcontext.VariantContext;

// methods for reading from and writing to VCFs
public final class VariantImpactSerialiser
{
    public static final String VAR_IMPACT_WORST = "VI_WORST";
    public static final String VAR_IMPACT_CANONICAL = "VI_CANON";

    public static VariantImpact fromVariantContext(final VariantContext context)
    {
        if(context.hasAttribute(VAR_IMPACT_WORST) && context.hasAttribute(VAR_IMPACT_CANONICAL))
        {
            final List<String> worst = context.getAttributeAsStringList(VAR_IMPACT_WORST, "");
            final List<String> canonical = context.getAttributeAsStringList(VAR_IMPACT_CANONICAL, "");
            return fromVcfAnnotation(worst, canonical);
        }

        // revert to SnpEff until migration is complete
        final List<String> worst = context.getAttributeAsStringList(SNPEFF_WORST, Strings.EMPTY);
        final List<String> canonical = context.getAttributeAsStringList(SNPEFF_CANONICAL, Strings.EMPTY);
        return fromVcfAnnotation(worst, canonical);
    }

    public static VariantImpact fromVcfAnnotation(final List<String> worst, final List<String> canonical)
    {
        int genesAffected = 0;
        String canonicalGene = "";
        String canonicalEffect = "";
        String canonicalTranscript = "";
        CodingEffect canonicalCodingEffect = UNDEFINED;
        String canonicalHgvsCodingImpact = "";
        String canonicalHgvsProteinImpact = "";
        String worstGene = "";
        String worstEffect = "";
        String worstTranscript = "";
        CodingEffect worstCodingEffect = UNDEFINED;

        if(worst.size() == 5)
        {
            worstGene = worst.get(0);
            worstTranscript = worst.get(1);
            worstEffect = readEffect(worst.get(2));
            worstCodingEffect = CodingEffect.valueOf(worst.get(3));
            genesAffected = Integer.parseInt(worst.get(4));
        }

        if(canonical.size() == 6)
        {
            canonicalGene = canonical.get(0);
            canonicalTranscript = canonical.get(1);
            canonicalEffect = canonical.get(2);
            canonicalCodingEffect = CodingEffect.valueOf(canonical.get(3));
            canonicalHgvsCodingImpact = canonical.get(4);
            canonicalHgvsProteinImpact = canonical.get(5);
        }

        return new VariantImpact(
                genesAffected, canonicalGene, canonicalEffect, canonicalTranscript, canonicalCodingEffect, canonicalHgvsCodingImpact,
                canonicalHgvsProteinImpact, worstGene, worstEffect, worstTranscript, worstCodingEffect);
    }

    public static List<String> worstDetails(final VariantImpact impact)
    {
        return Lists.newArrayList(
                impact.WorstGene,
                impact.WorstTranscript,
                writeEffect(impact.WorstEffect),
                impact.WorstCodingEffect.toString(),
                String.valueOf(impact.GenesAffected));
    }

    public static List<String> canonicalDetails(final VariantImpact impact)
    {
        return Lists.newArrayList(
                impact.CanonicalGene,
                impact.CanonicalTranscript,
                writeEffect(impact.CanonicalEffect),
                impact.CanonicalCodingEffect.toString(),
                impact.CanonicalHgvsCodingImpact,
                impact.CanonicalHgvsProteinImpact);
    }

    public static String writeEffect(final String effect)
    {
        return effect.replace("; ", "&").replace(" ", "_");
    }

    public static String readEffect(final String effect)
    {
        return effect.replace("&", "; ").replace("_", " ");
    }

}
