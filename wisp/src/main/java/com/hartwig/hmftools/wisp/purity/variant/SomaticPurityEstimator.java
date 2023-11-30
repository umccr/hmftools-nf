package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.stats.PoissonCalcs.calcPoissonNoiseValue;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_COUNT_MODEL_MIN_AVG_DEPTH;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_COUNT_MODEL_MIN_FRAG_VARIANTS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.VAF_PEAK_MODEL_MIN_AVG_DEPTH;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.VAF_PEAK_MODEL_MIN_FRAG_VARIANTS;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityMethod.isRecomputed;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.LOW_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.expectedNoise;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.estimatedPurity;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityResult.INVALID_RESULT;

import java.util.List;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.wisp.common.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;

public class SomaticPurityEstimator
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;
    private final SampleData mSample;

    public SomaticPurityEstimator(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;
    }

    public SomaticPurityResult calculatePurity(
            final String sampleId, final PurityContext purityContext, final List<SomaticVariant> variants)
    {
        // TO-DO - pass in or set later
        int totalVariants = 0;

        FragmentTotals fragmentTotals = new FragmentTotals();

        UmiTypeCounts umiTypeCounts = new UmiTypeCounts();

        for(SomaticVariant variant : variants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            fragmentTotals.addVariantData(
                    variant.copyNumber(), tumorFragData.AlleleCount, sampleFragData.AlleleCount,
                    tumorFragData.Depth, sampleFragData.Depth, sampleFragData.QualTotal);

            umiTypeCounts.add(sampleFragData.UmiCounts);
        }

        int sampleDepthTotal = fragmentTotals.sampleAdTotal();
        if(sampleDepthTotal == 0)
            return INVALID_RESULT;

        double tumorPurity = purityContext.bestFit().purity();

        /*
        double tumorVaf;

        if(!mConfig.hasSyntheticTumor())
        {
            double tumorDepthTotal = tumorCounts.totalFragments();
            if(tumorDepthTotal == 0)
                return INVALID_RESULT;

            tumorVaf = tumorCounts.alleleFragments() / tumorDepthTotal;
        }
        else
        {
            tumorVaf = SYNTHETIC_TUMOR_VAF;
        }
        */

        double expectedSampleNoise = expectedNoise(fragmentTotals, mConfig.NoiseReadsPerMillion);

        PurityCalcData purityCalcData = new PurityCalcData();

        // firstly estimate raw purity without consideration of clonal peaks
        purityCalcData.RawPurityEstimate = estimatedPurity(
                tumorPurity, fragmentTotals.adjTumorVaf(), fragmentTotals.adjSampleVaf(), expectedSampleNoise);

        double weightedAvgDepth = fragmentTotals.weightedSampleDepth();

        ClonalityModel model = null;

        if(fragmentTotals.sampleTwoPlusCount() >= VAF_PEAK_MODEL_MIN_FRAG_VARIANTS
        && weightedAvgDepth > VAF_PEAK_MODEL_MIN_AVG_DEPTH)
        {
            model = new VafPeakModel(mConfig, mResultsWriter, mSample, variants);
        }
        else if(fragmentTotals.sampleOneFragmentCount() + fragmentTotals.sampleTwoPlusCount() >= LOW_COUNT_MODEL_MIN_FRAG_VARIANTS
        && weightedAvgDepth < LOW_COUNT_MODEL_MIN_AVG_DEPTH)
        {
            model = new LowCountModel(mConfig, mResultsWriter, mSample, variants);
        }

        if(model != null)
        {
            purityCalcData.Clonality = model.calculate(sampleId, fragmentTotals, purityCalcData.RawPurityEstimate);

            if(isRecomputed(purityCalcData.Clonality.Method))
            {
                purityCalcData.PurityEstimate = estimatedPurity(
                        tumorPurity, fragmentTotals.adjTumorVaf(), purityCalcData.Clonality.Vaf, expectedSampleNoise);

                purityCalcData.PurityRangeLow = estimatedPurity(
                        tumorPurity, fragmentTotals.adjTumorVaf(), purityCalcData.Clonality.VafLow, expectedSampleNoise);

                purityCalcData.PurityRangeHigh = estimatedPurity(
                        tumorPurity, fragmentTotals.adjTumorVaf(), purityCalcData.Clonality.VafHigh, expectedSampleNoise);
            }
        }

        // calculate a limit-of-detection (LOD), being the number of fragments that would return a 95% confidence of a tumor presence
        double lodFragments = calcPoissonNoiseValue((int)round(expectedSampleNoise), LOW_PROBABILITY);
        double lodSampleVaf = max(lodFragments - expectedSampleNoise, 0) / (double)fragmentTotals.sampleDepthTotal();

        purityCalcData.LodPurityEstimate = estimatedPurity(tumorPurity, fragmentTotals.adjTumorVaf(), lodSampleVaf, expectedSampleNoise);


        /*
        double dualFragsNoise = sampleCountsDual.totalFragments() / 1000000.0 * mConfig.NoiseReadsPerMillionDualStrand + lowQualNoiseFactor;

        FragmentCalcResult dualFragsResult = SomaticPurityCalcs.calc(
                tumorPloidy, adjustedTumorVaf, sampleCountsDual.totalFragments(), sampleCountsDual.alleleFragments(), dualFragsNoise);


        FragmentCalcResult lodFragsResult = SomaticPurityCalcs.calc(
                tumorPloidy, adjustedTumorVaf, sampleDepthTotal, (int)round(lodFragments), expectedSampleNoise);
        */

        // CT_LOGGER.info(format("patient(%s) sample(%s) sampleTotalFrags(%d) noise(%.1f) LOD(%.6f)",
        //        mSample.PatientId, sampleId, sampleDepthTotal, allFragsNoise, lodFragsResult.EstimatedPurity));

        return new SomaticPurityResult(true, totalVariants, fragmentTotals, umiTypeCounts, purityCalcData);
    }

    private void removeChipVariants(final FragmentTotals fragmentTotals)
    {

    }
}
