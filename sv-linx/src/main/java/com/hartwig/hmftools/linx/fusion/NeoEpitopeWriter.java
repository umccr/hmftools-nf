package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.getExponent;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_NONSENSE_MED_DECAY;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.switchStream;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.DELIMITER;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_KNOWN;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.MAX_UPSTREAM_DISTANCE_OTHER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.fusion.TranscriptRegionType;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.linx.types.LinkedPair;

import org.apache.commons.compress.utils.Lists;

public class NeoEpitopeWriter
{
    private final String mOutputDir;

    private boolean mIsMultiSample;
    private BufferedWriter mFileWriter;
    private String mSampleId;

    private final EnsemblDataCache mGeneTransCache;
    private final KnownFusionCache mKnownFusionCache;
    private final List<NeoEpitopeFusion> mFusions;

    public NeoEpitopeWriter(
            final String outputDir, boolean isMultiSample, final EnsemblDataCache geneDataCache, final KnownFusionCache knownFusionCache)
    {
        mOutputDir = outputDir;
        mIsMultiSample = isMultiSample;
        mGeneTransCache = geneDataCache;
        mKnownFusionCache = knownFusionCache;

        mFileWriter = null;
        mSampleId = null;
        mFusions = Lists.newArrayList();
    }

    public void initialiseSample(final String sampleId)
    {
        mSampleId = sampleId;

        // clear any cache
        mFusions.clear();
    }

    public void processFusionCandidate(
            final List<BreakendGeneData> breakendGenes1, final List<BreakendGeneData> breakendGenes2,
            final List<LinkedPair> traversedPairs, final DisruptionFinder disruptionFinder,
            final LinkedPair lowerLink, final LinkedPair upperLink)
    {
        if(breakendGenes1.isEmpty() || breakendGenes2.isEmpty())
            return;

        int[] extensionLengths = new int[FS_PAIR];

        int chainLength = traversedPairs != null ? traversedPairs.stream().mapToInt(x -> x.length()).sum() : 0;

        boolean fusionAdded = false;

        for (final BreakendGeneData gene1 : breakendGenes1)
        {
            for (final BreakendGeneData gene2 : breakendGenes2)
            {
                if (gene1.isUpstream() == gene2.isUpstream())
                    continue;

                final BreakendGeneData upGene = gene1.isUpstream() ? gene1 : gene2;

                if(!isCandidateUpstreamGene(upGene))
                    continue;

                final BreakendGeneData downGene = upGene == gene1 ? gene2 : gene1;

                if(!hasValidTraversal(upGene, traversedPairs, gene1.isUpstream(), disruptionFinder))
                    continue;

                int geneStreamIndex = gene1.isUpstream() ? FS_UP : FS_DOWN;
                extensionLengths[geneStreamIndex] = lowerLink != null ? lowerLink.length() : -1;
                extensionLengths[switchStream(geneStreamIndex)] = upperLink != null ? upperLink.length() : -1;

                int maxUpstreamDistance = getMaxUpstreamDistance(upGene.GeneName, downGene.GeneName);

                final List<BreakendTransData> validUpTrans = findValidTranscripts(upGene, extensionLengths[FS_UP], 0);
                final List<BreakendTransData> validDownTrans = findValidTranscripts(downGene, extensionLengths[FS_DOWN], maxUpstreamDistance);

                if(validUpTrans.isEmpty() || validDownTrans.isEmpty())
                    continue;

                // avoid writing duplicates of the same junction, including in the reverse direction
                if(isDuplicate(breakendGenes1.get(0), breakendGenes2.get(0)) || isDuplicate(breakendGenes2.get(0), breakendGenes1.get(0)))
                    continue;

                double avgJcn = (upGene.jcn() + downGene.jcn()) * 0.5;

                final StringJoiner sjUp = new StringJoiner(ITEM_DELIM);
                validUpTrans.stream().map(x -> x.TransData.TransName).forEach(x -> sjUp.add(x));

                final StringJoiner sjDown = new StringJoiner(ITEM_DELIM);
                validDownTrans.stream().map(x -> x.TransData.TransName).forEach(x -> sjDown.add(x));

                NeoEpitopeFusion fusion = new NeoEpitopeFusion(
                        upGene.StableId, upGene.GeneName, upGene.chromosome(), upGene.position(), upGene.orientation(), upGene.id(),
                        downGene.StableId, downGene.GeneName, downGene.chromosome(), downGene.position(), downGene.orientation(), downGene.id(),
                        avgJcn, upGene.insertSequence(), chainLength, new String[] { sjUp.toString(), sjDown.toString()});

                writeData(fusion);

                if(!fusionAdded)
                {
                    mFusions.add(fusion);
                    fusionAdded = true;
                }
            }
        }
    }

    private final List<BreakendTransData> findValidTranscripts(final BreakendGeneData gene, int linkExtensionLength, int maxUpstreamDistance)
    {
        final List<BreakendTransData> validTrans = Lists.newArrayList();

        for(final BreakendTransData transcript : gene.transcripts())
        {
            if(gene.isUpstream())
            {
                if(!transcript.isDisruptive())
                    continue;

                if(transcript.TransData.CodingStart == null)
                    continue;
            }
            else
            {
                // must have a splice acceptor
                if(transcript.TransData.exons().size() <= 1)
                    continue;

                if(transcript.regionType() == EXONIC && transcript.ExonDownstream == transcript.exonCount())
                    continue;

                // cannot be NMD
                if(transcript.TransData.BioType.equals(BIOTYPE_NONSENSE_MED_DECAY))
                    continue;

                if(transcript.getDistanceUpstream() > maxUpstreamDistance)
                    continue;

                // check for a preceding splice acceptor from another transcript, whether same gene or not
                int preTransDistance = gene.Strand == POS_STRAND ?
                        transcript.transStart() - gene.position() : gene.position() - transcript.transEnd();

                if(preTransDistance > 0)
                {
                    int preTransSpliceAcceptorPos = mGeneTransCache.findPrecedingGeneSpliceAcceptorPosition(transcript.transId());

                    if(preTransSpliceAcceptorPos > 0)
                    {
                        int preTransSpliceAcceptorDistance = gene.Strand == POS_STRAND ?
                                transcript.transStart() - preTransSpliceAcceptorPos : preTransSpliceAcceptorPos - transcript.transEnd();

                        if(preTransSpliceAcceptorDistance < preTransDistance)
                            continue;
                    }
                }
            }

            // ensure the transcript is not interrupted by a chain link going elsewhere
            if(linkExtensionLength > 0)
            {
                int transLength = (gene.Strand == POS_STRAND) == gene.isUpstream() ?
                        abs(gene.position() - transcript.TransData.TransStart) : abs(gene.position() - transcript.TransData.TransEnd);

                if(linkExtensionLength < transLength)
                    continue;
            }

            validTrans.add(transcript);
        }

        return validTrans;
    }

    private int getMaxUpstreamDistance(final String upGene, final String downGene)
    {
        if(mKnownFusionCache.hasKnownFusion(upGene, downGene))
            return MAX_UPSTREAM_DISTANCE_KNOWN;

        if(mKnownFusionCache.hasPromiscuousFiveGene(upGene) && mKnownFusionCache.isHighImpactPromiscuous(PROMISCUOUS_5, upGene, downGene))
            return MAX_UPSTREAM_DISTANCE_KNOWN;

        if(mKnownFusionCache.hasPromiscuousThreeGene(downGene) && mKnownFusionCache.isHighImpactPromiscuous(PROMISCUOUS_3, upGene, downGene))
            return MAX_UPSTREAM_DISTANCE_KNOWN;

        return MAX_UPSTREAM_DISTANCE_OTHER;
    }

    private boolean hasValidTraversal(
            final BreakendGeneData upGene, final List<LinkedPair> traversedPairs,
            boolean fusionLowerToUpper, final DisruptionFinder disruptionFinder)
    {
        if(traversedPairs == null || traversedPairs.isEmpty())
            return true;

        int upGeneStrand = upGene.Strand;
        boolean isPrecodingUpstream = false;

        for(LinkedPair pair : traversedPairs)
        {
            // if going lower to upper, if the orientation of the first breakend in the pair is opposite to the strand of
            // the upstream gene, then the fusion direction for that pair is the same as a the upstream gene
            // otherwise it needs to be switched
            int fusionDirection = 0;

            if(fusionLowerToUpper)
            {
                fusionDirection = pair.firstBreakend().orientation() != upGeneStrand ? upGeneStrand : -upGeneStrand;
            }
            else
            {
                fusionDirection = pair.secondBreakend().orientation() != upGeneStrand ? upGeneStrand : -upGeneStrand;
            }

            // any invalid traversal causes this fusion to be entirely skipped from further analysis
            if(disruptionFinder.pairTraversesGene(pair, fusionDirection, isPrecodingUpstream))
                return false;
        }

        return true;
    }

    private boolean isCandidateUpstreamGene(final BreakendGeneData gene)
    {
        // must have at least 1 coding transcript
        return gene.transcripts().stream().anyMatch(x -> x.isCoding() && x.isDisruptive());
    }

    private boolean isDuplicate(final BreakendGeneData gene1, final BreakendGeneData gene2)
    {
        for(final NeoEpitopeFusion fusion : mFusions)
        {
            for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
            {
                if(fusion.SvIds[fs] == gene1.id() && fusion.SvIds[switchStream(fs)] == gene2.id()
                && fusion.Chromosomes[fs].equals(gene1.chromosome()) && fusion.Chromosomes[switchStream(fs)].equals(gene2.chromosome())
                && fusion.Positions[fs] == gene1.position() && fusion.Positions[switchStream(fs)] == gene2.position())
                {
                    return true;
                }
            }
        }

        return false;
    }

    private void writeData(final NeoEpitopeFusion fusion)
    {
        try
        {
            if(mFileWriter == null)
            {
                if(mIsMultiSample)
                {
                    String outputFileName = mOutputDir + "LNX_NEO_EPITOPES.csv";

                    mFileWriter = createBufferedWriter(outputFileName, false);

                    mFileWriter.write("sampleId");
                    mFileWriter.write(NeoEpitopeFusion.header());
                    mFileWriter.newLine();
                }
                else
                {
                    String outputFileName = NeoEpitopeFusion.generateFilename(mOutputDir, mSampleId);

                    mFileWriter = createBufferedWriter(outputFileName, false);
                    mFileWriter.write(NeoEpitopeFusion.header());
                    mFileWriter.newLine();
                }
            }

            if(mIsMultiSample)
            {
                mFileWriter.write(String.format("%s",mSampleId));
            }

            mFileWriter.write(NeoEpitopeFusion.toString(fusion));
            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing neo-epitope output file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }

}
