package com.hartwig.hmftools.neo.scorer;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.rna.RnaExpressionMatrix.EXPRESSION_SCOPE_TRANS;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.scorer.NeoScorerConfig.COHORT_SAMPLE_TPM_FILE;
import static com.hartwig.hmftools.neo.scorer.NeoScorerConfig.COHORT_TPM_MEDIANS_FILE;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.rna.RnaExpressionMatrix;
import com.hartwig.hmftools.neo.bind.BindScorer;
import com.hartwig.hmftools.neo.bind.ScoreConfig;

import org.apache.commons.cli.CommandLine;

public class ReferenceData
{
    public final BindScorer PeptideScorer;
    public final Map<String,TranscriptAminoAcids> TransAminoAcidMap;
    public final RnaExpressionMatrix TranscriptExpression;
    public final TpmMediansCache TpmMedians;

    public ReferenceData(final CommandLine cmd)
    {
        PeptideScorer = new BindScorer(new ScoreConfig(cmd));

        String cohortSampleTpmFile = cmd.getOptionValue(COHORT_SAMPLE_TPM_FILE);

        if(cohortSampleTpmFile != null)
        {
            NE_LOGGER.info("loading cohort transcript expression");
            TranscriptExpression = new RnaExpressionMatrix(cohortSampleTpmFile, EXPRESSION_SCOPE_TRANS);
        }
        else
        {
            TranscriptExpression = null;
        }

        NE_LOGGER.info("loading cohort transcript medians");

        String cohortTpmMediansFile = cmd.getOptionValue(COHORT_TPM_MEDIANS_FILE);
        TpmMedians = new TpmMediansCache(cohortTpmMediansFile);

        TransAminoAcidMap = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(
                cmd.getOptionValue(ENSEMBL_DATA_DIR), TransAminoAcidMap, Lists.newArrayList(), false);
    }

    public boolean peptideMatchesWildtype(final String peptide, final String transName)
    {
        TranscriptAminoAcids transAAs = TransAminoAcidMap.get(transName);

        if(transAAs == null)
            return false;

        return transAAs.AminoAcids.contains(peptide);
    }
}