package com.hartwig.hmftools.markdups.umi;

public enum ConsensusOutcome
{
    UNSET,
    ALIGNMENT_ONLY,
    INDEL_MATCH,
    INDEL_MISMATCH,
    INDEL_FAIL;

    public boolean valid() { return this != INDEL_FAIL; }
}