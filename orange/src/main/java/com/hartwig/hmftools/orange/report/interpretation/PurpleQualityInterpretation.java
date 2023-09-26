package com.hartwig.hmftools.orange.report.interpretation;

import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;

import org.jetbrains.annotations.NotNull;

public final class PurpleQualityInterpretation
{
    public static boolean isQCFail(@NotNull PurpleQC purpleQC)
    {
        return purpleQC.status().contains(PurpleQCStatus.FAIL_NO_TUMOR) || purpleQC.status().contains(PurpleQCStatus.FAIL_CONTAMINATION);
    }
}
