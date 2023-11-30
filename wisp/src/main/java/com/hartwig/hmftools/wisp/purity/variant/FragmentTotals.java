package com.hartwig.hmftools.wisp.purity.variant;

public class FragmentTotals
{
    private int mVariantCount;
    private int mTumorAdTotal;
    private int mSampleAdTotal;
    private double mSampleAllelelQualTotal;

    private int mTumorDepthTotal;
    private int mSampleDepthTotal;

    private int mSampleOneFragmentCount; // count of variants with 1 observed fragment
    private int mSampleTwoPlusCount;

    private double mSampleAdjustedAdTotal;
    private double mSampleAdjustedDepthTotal;

    private double mSampleWeightedDepthTotal;

    private double mTumorAdjustedAdTotal;
    private double mTumorAdjustedDepthTotal;

    public FragmentTotals()
    {
        mVariantCount = 0;
        mTumorAdTotal = 0;
        mSampleAdTotal = 0;
        mSampleAllelelQualTotal = 0;
        mTumorDepthTotal = 0;
        mSampleDepthTotal = 0;
        mSampleOneFragmentCount = 0;
        mSampleTwoPlusCount = 0;
        mSampleAdjustedAdTotal = 0;
        mSampleAdjustedDepthTotal = 0;
        mSampleWeightedDepthTotal = 0;
        mTumorAdjustedAdTotal = 0;
        mTumorAdjustedDepthTotal = 0;
    }

    public void addVariantData(
            double copyNumber, int tumorAlleleFrags, int sampleAlleleFrags, int tumorDepth, int sampleDepth, double sampleQualTotal)
    {
        ++mVariantCount;
        mTumorAdTotal += tumorAlleleFrags;
        mTumorDepthTotal += tumorDepth;

        mSampleAdTotal += sampleAlleleFrags;
        mSampleDepthTotal += sampleDepth;
        mSampleAllelelQualTotal += sampleQualTotal;

        if(sampleAlleleFrags >= 2)
            ++mSampleTwoPlusCount;
        else if(sampleAlleleFrags == 1)
            ++mSampleOneFragmentCount;

        // wVAF = Σ(i=1->n)[ADi /CNi] * Σ(i=1->n)[DPi /CNi]

        mTumorAdjustedAdTotal += tumorAlleleFrags / copyNumber;
        mSampleAdjustedAdTotal += sampleAlleleFrags / copyNumber;

        double tumorDpPerCn = tumorDepth / copyNumber;
        mTumorAdjustedDepthTotal += tumorDpPerCn;
        mSampleAdjustedDepthTotal += sampleDepth / copyNumber;

        mSampleWeightedDepthTotal += tumorDpPerCn * sampleDepth;
    }

    public int variantCount() { return mVariantCount; }
    public int tumorAdTotal() { return mTumorAdTotal; }
    public int sampleAdTotal() { return mSampleAdTotal; }
    public int tumorDepthTotal() { return mTumorDepthTotal; }
    public int sampleDepthTotal() { return mSampleDepthTotal; }
    public int sampleOneFragmentCount() { return mSampleOneFragmentCount; }
    public int sampleTwoPlusCount() { return mSampleTwoPlusCount; }

    public double rawTumorVaf() { return mTumorDepthTotal > 0 ? mTumorAdTotal / mTumorDepthTotal : 0; }

    public double adjTumorVaf() { return mTumorAdjustedDepthTotal > 0 ? mTumorAdjustedAdTotal / mTumorAdjustedDepthTotal : 0; }

    public double adjSampleVaf() { return mSampleAdjustedDepthTotal > 0 ? mSampleAdjustedAdTotal / mSampleAdjustedDepthTotal : 0; }

    public double weightedSampleDepth()
    {
        // wAD = Σ(i=1->n)[DPi_cfDNA * DPi_tissue/CNn] / Σ(i=1->n)[DPi_Tissue / CNi]
        return mTumorAdjustedDepthTotal > 0 ? mSampleWeightedDepthTotal / mTumorAdjustedDepthTotal : 0;
    }

    public double qualPerAllele()
    {
        return mSampleAdTotal > 0 ? mSampleAllelelQualTotal / mSampleAdTotal : 0;
    }
}
