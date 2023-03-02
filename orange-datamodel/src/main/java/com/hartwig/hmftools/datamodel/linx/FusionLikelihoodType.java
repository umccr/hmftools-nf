package com.hartwig.hmftools.datamodel.linx;

public enum FusionLikelihoodType
{
    HIGH,
    LOW,
    NA;

    public String displayStr()
    {
        switch(this)
        {
            case HIGH: return "High";
            case LOW: return "Low";
            case NA: return "NA";
            default: return "Invalid";
        }
    }
}
