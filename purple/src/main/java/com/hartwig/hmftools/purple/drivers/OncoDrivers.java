package com.hartwig.hmftools.purple.drivers;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverImpact;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.ReportablePredicate;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;

public class OncoDrivers
{
    public static final int MAX_REPEAT_COUNT = 7;

    private final ReportablePredicate mReportablePredicate;
    private final Map<String, DndsDriverGeneLikelihood> mLikelihoodsByGene;

    public OncoDrivers(final DriverGenePanel genePanel)
    {
        mLikelihoodsByGene = genePanel.oncoLikelihood();
        mReportablePredicate = new ReportablePredicate(DriverCategory.ONCO, genePanel.driverGenes());
    }

    public List<DriverCatalog> drivers(
            final List<SomaticVariant> variants, final Map<String,List<GeneCopyNumber>> geneCopyNumberMap,
            final Map<VariantType,Integer> variantTypeCounts)
    {
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        int sampleSNVCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0);
        int sampleINDELCount = variantTypeCounts.getOrDefault(VariantType.INDEL, 0);

        final Map<String,List<SomaticVariant>> codingVariants = oncogenicVariantsByGene(variants);

        for(String gene : codingVariants.keySet())
        {
            final DndsDriverGeneLikelihood geneMissenseLikelihood = mLikelihoodsByGene.get(gene);

            final List<SomaticVariant> geneVariants = codingVariants.get(gene);

            List<GeneCopyNumber> geneCopyNumbers = geneCopyNumberMap.get(gene);

            if(geneCopyNumbers == null)
                continue;

            for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
            {
                driverCatalog.add(geneDriver(
                        sampleSNVCount, sampleINDELCount, geneMissenseLikelihood, geneVariants, geneCopyNumber));
            }

        }

        return driverCatalog;
    }

    private Map<String, List<SomaticVariant>> oncogenicVariantsByGene(final List<SomaticVariant> variants)
    {
        return variants.stream().filter(mReportablePredicate).collect(Collectors.groupingBy(SomaticVariant::gene));
    }

    public static DriverCatalog geneDriver(
            int sampleSNVCount, int sampleIndelCount, final DndsDriverGeneLikelihood geneLikelihood,
            final List<SomaticVariant> geneVariants, final GeneCopyNumber geneCopyNumber)
    {
        final Map<DriverImpact,Integer> variantCounts = DriverCatalogFactory.driverImpactCount(geneVariants);
        int missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0);
        int nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0);
        int spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0);
        int inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0);
        int frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber == null ? Strings.EMPTY : geneCopyNumber.chromosomeBand())
                .gene(geneLikelihood.gene())
                .transcript(geneCopyNumber != null ? geneCopyNumber.transName() : "")
                .isCanonical(geneCopyNumber != null ? geneCopyNumber.isCanonical() : true)
                .driver(DriverType.MUTATION)
                .category(DriverCategory.ONCO)
                .driverLikelihood(1)
                .missense(missenseVariants)
                .nonsense(nonsenseVariants)
                .splice(spliceVariants)
                .inframe(inframeVariants)
                .frameshift(frameshiftVariants)
                .biallelic(geneVariants.stream().anyMatch(SomaticVariant::biallelic))
                .minCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.maxCopyNumber())
                .likelihoodMethod(LikelihoodMethod.DNDS);

        if(geneVariants.stream().anyMatch(SomaticVariant::isHotspot))
        {
            return builder.likelihoodMethod(LikelihoodMethod.HOTSPOT).build();
        }

        if(geneVariants.stream().anyMatch(OncoDrivers::isKnownInframeIndel))
        {
            return builder.likelihoodMethod(LikelihoodMethod.INFRAME).build();
        }

        double driverLikelihood = 0;

        if(geneLikelihood != null)
        {
            for(SomaticVariant variant : geneVariants)
            {
                final DriverImpact impact = DriverImpact.select(variant);

                final DndsDriverImpactLikelihood likelihood = geneLikelihood.select(impact);

                final int sampleVariantCount =
                        impact == DriverImpact.FRAMESHIFT || impact == DriverImpact.INFRAME ? sampleIndelCount : sampleSNVCount;

                driverLikelihood = Math.max(driverLikelihood, DriverCatalogFactory.probabilityDriverVariant(sampleVariantCount, likelihood));
            }
        }

        return builder.driverLikelihood(driverLikelihood).build();
    }

    private static boolean isKnownInframeIndel(final SomaticVariant variant)
    {
        return variant.type() == VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.MISSENSE
                && variant.repeatCount() <= MAX_REPEAT_COUNT;
    }
}