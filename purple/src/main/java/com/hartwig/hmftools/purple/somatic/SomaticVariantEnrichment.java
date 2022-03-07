package com.hartwig.hmftools.purple.somatic;

import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.common.variant.enrich.SomaticPurityEnrichment;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichmentFactory;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantEnrichment implements VariantContextEnrichment
{
    private final VariantContextEnrichment mPurityEnrichment;
    private final VariantContextEnrichment mHotspotEnrichment;
    private final VariantContextEnrichment mKataegisEnrichment;
    private final VariantContextEnrichment mSomaticRefContextEnrichment;
    private final SubclonalLikelihoodEnrichment mSubclonalLikelihoodEnrichment;
    private final VariantContextEnrichment mSnpEffEnrichment;
    private final SomaticGenotypeEnrichment mGenotypeEnrichment;

    public SomaticVariantEnrichment(
            boolean hotspotEnabled, boolean snpEffEnrichmentEnabled, double clonalityBinWidth, final String purpleVersion,
            final String referenceId, final String tumorSample, final ReferenceData refData,
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, final List<FittedRegion> fittedRegions,
            final Multimap<Chromosome, VariantHotspot> hotspots, final List<PeakModel> peakModel, final Consumer<VariantContext> consumer)
    {
        mGenotypeEnrichment = new SomaticGenotypeEnrichment(referenceId, tumorSample, consumer);

        mSubclonalLikelihoodEnrichment = new SubclonalLikelihoodEnrichment(clonalityBinWidth, peakModel, mGenotypeEnrichment);

        mPurityEnrichment = new SomaticPurityEnrichment(
                purpleVersion, tumorSample, purityAdjuster, copyNumbers, fittedRegions, mSubclonalLikelihoodEnrichment);

        mKataegisEnrichment = new KataegisEnrichment(mPurityEnrichment);

        mSomaticRefContextEnrichment = new SomaticRefContextEnrichment(refData.RefGenome, mKataegisEnrichment);

        if(snpEffEnrichmentEnabled)
        {
            final Set<String> somaticGenes = refData.DriverGenes.driverGenes().stream()
                    .filter(DriverGene::reportSomatic).map(DriverGene::gene).collect(Collectors.toSet());

            mSnpEffEnrichment = new SnpEffEnrichment(
                    somaticGenes, refData.GeneTransCache, refData.OtherReportableTranscripts, mSomaticRefContextEnrichment);
        }
        else
        {
            mSnpEffEnrichment = null;
        }

        VariantContextEnrichment prevConsumer = mSnpEffEnrichment != null ? mSnpEffEnrichment : mSomaticRefContextEnrichment;

        if(hotspotEnabled)
        {
            mHotspotEnrichment = new VariantHotspotEnrichment(hotspots, prevConsumer);
        }
        else
        {
            mHotspotEnrichment = VariantContextEnrichmentFactory.noEnrichment().create(prevConsumer);
        }
    }

    @Override
    public void flush()
    {
        mHotspotEnrichment.flush();

        if(mSnpEffEnrichment != null)
            mSnpEffEnrichment.flush();

        mSomaticRefContextEnrichment.flush();
        mKataegisEnrichment.flush();
        mPurityEnrichment.flush();
        mSubclonalLikelihoodEnrichment.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template)
    {
        VCFHeader header = mSomaticRefContextEnrichment.enrichHeader(template);
        header = mKataegisEnrichment.enrichHeader(header);
        header = mSubclonalLikelihoodEnrichment.enrichHeader(header);
        header = mHotspotEnrichment.enrichHeader(header);

        if(mSnpEffEnrichment != null)
            header = mSnpEffEnrichment.enrichHeader(header);

        return mPurityEnrichment.enrichHeader(header);
    }

    @Override
    public void accept(@NotNull final VariantContext context)
    {
        mHotspotEnrichment.accept(context);
    }
}