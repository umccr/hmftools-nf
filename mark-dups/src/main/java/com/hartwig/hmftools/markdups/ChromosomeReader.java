package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.common.DuplicateGroupBuilder.findDuplicateFragments;
import static com.hartwig.hmftools.markdups.common.FilterReadsType.readOutsideSpecifiedRegions;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.markdups.common.CandidateDuplicates;
import com.hartwig.hmftools.markdups.common.DuplicateGroupBuilder;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.FragmentStatus;
import com.hartwig.hmftools.markdups.common.FragmentUtils;
import com.hartwig.hmftools.markdups.common.PartitionData;
import com.hartwig.hmftools.markdups.common.PartitionResults;
import com.hartwig.hmftools.markdups.common.Statistics;
import com.hartwig.hmftools.markdups.common.UnmapStats;
import com.hartwig.hmftools.markdups.consensus.ConsensusReads;
import com.hartwig.hmftools.markdups.common.DuplicateGroup;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ChromosomeReader implements Consumer<List<Fragment>>, Callable
{
    private final MarkDupsConfig mConfig;
    private final ChrBaseRegion mRegion;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final PartitionDataStore mPartitionDataStore;
    private final BamWriter mBamWriter;
    private final ReadPositionsCache mReadPositions;
    private final DuplicateGroupBuilder mDuplicateGroupBuilder;
    private final ConsensusReads mConsensusReads;

    private BaseRegion mCurrentPartition;
    private String mCurrentStrPartition;
    private PartitionData mCurrentPartitionData;
    private ChrBaseRegion mExcludedRegion;
    private final List<BaseRegion> mUnmapRegions;

    private final Map<String,List<SAMRecord>> mPendingIncompleteReads;

    private final boolean mLogReadIds;
    private int mPartitionRecordCount;
    private final Statistics mStats;
    private final PerformanceCounter mPcTotal;
    private final PerformanceCounter mPcAcceptPositions;
    private final PerformanceCounter mPcPendingIncompletes;

    public ChromosomeReader(
            final ChrBaseRegion region, final MarkDupsConfig config, FileWriterCache fileWriterCache,
            final PartitionDataStore partitionDataStore)
    {
        mConfig = config;
        mRegion = region;
        mPartitionDataStore = partitionDataStore;
        mBamWriter = fileWriterCache.getBamWriter(mRegion.Chromosome);

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, false);
        mBamSlicer.setKeepUnmapped();

        mReadPositions = new ReadPositionsCache(region.Chromosome, config.BufferSize, !config.NoMateCigar, this);
        mDuplicateGroupBuilder = new DuplicateGroupBuilder(config);
        mConsensusReads = new ConsensusReads(config.RefGenome);
        mConsensusReads.setDebugOptions(config.RunChecks);

        mUnmapRegions = Lists.newArrayList();

        if(!mConfig.SpecificChrRegions.Regions.isEmpty())
        {
            // NOTE: doesn't currently handle multiple regions on the same chromosome
            ChrBaseRegion firstRegion = mConfig.SpecificChrRegions.Regions.stream()
                    .filter(x -> x.Chromosome.equals(mRegion.Chromosome)).findFirst().orElse(mRegion);

            int partitionStart = (firstRegion.start() / mConfig.PartitionSize) * mConfig.PartitionSize;
            mCurrentPartition = new BaseRegion(partitionStart, partitionStart + mConfig.PartitionSize - 1);
        }
        else
        {
            mCurrentPartition = new BaseRegion(1, mConfig.PartitionSize - 1);
        }

        setUnmappedRegions();

        mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);
        mCurrentPartitionData = mPartitionDataStore.getOrCreatePartitionData(mCurrentStrPartition);
        setExcludedRegion(ExcludedRegions.getPolyGRegion(mConfig.RefGenVersion));

        mPendingIncompleteReads = Maps.newHashMap();

        mPartitionRecordCount = 0;

        mStats = mDuplicateGroupBuilder.statistics();

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
        mPcTotal = new PerformanceCounter("Total");
        mPcAcceptPositions = new PerformanceCounter("AcceptPositions");
        mPcPendingIncompletes = new PerformanceCounter("PendingIncompletes");
    }

    // public PerformanceCounter perfCounter() { return mPcTotal; }
    public List<PerformanceCounter> perfCounters()
    {
        return List.of(mPcTotal, mPcAcceptPositions, mPcPendingIncompletes);
    }

    public Statistics statistics() { return mStats; }
    public BamWriter recordWriter() { return mBamWriter; }

    @Override
    public Long call()
    {
        run();
        return (long)1;
    }

    public void run()
    {
        perfCountersStart();

        if(!mConfig.SpecificChrRegions.Regions.isEmpty())
        {
            for(ChrBaseRegion region : mConfig.SpecificChrRegions.Regions)
            {
                if(!region.Chromosome.equals(mRegion.Chromosome))
                    continue;

                MD_LOGGER.debug("processing specific region({})", region);
                mBamSlicer.slice(mSamReader, region, this::processSamRecord);
            }
        }
        else
        {
            MD_LOGGER.info("processing chromosome({})", mRegion.Chromosome);
            mBamSlicer.slice(mSamReader, mRegion, this::processSamRecord);
        }

        onPartitionComplete(false);

        MD_LOGGER.info("chromosome({}) complete, reads({})", mRegion.Chromosome, mStats.TotalReads);

        mConsensusReads.logStats(mRegion.Chromosome);
    }

    private void onPartitionComplete(boolean setupNext)
    {
        mReadPositions.evictAll();

        processPendingIncompletes();

        perfCountersStop();

        MD_LOGGER.debug("partition({}:{}) complete, reads({})", mRegion.Chromosome, mCurrentPartition, mPartitionRecordCount);

        if(mConfig.PerfDebug)
            mCurrentPartitionData.logCacheCounts();

        mPartitionRecordCount = 0;

        if(setupNext)
        {
            // move ahead to the next partition, until the end of the chromosome is reached
            int regionStart = mCurrentPartition.end() + 1;

            if(regionStart > mRegion.end())
            {
                mCurrentPartition = null;
                return;
            }

            mCurrentPartition.setStart(regionStart);
            mCurrentPartition.setEnd(regionStart + mConfig.PartitionSize - 1);

            mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);
            mCurrentPartitionData = mPartitionDataStore.getOrCreatePartitionData(mCurrentStrPartition);
            setExcludedRegion(ExcludedRegions.getPolyGRegion(mConfig.RefGenVersion));
            setUnmappedRegions();

            perfCountersStart();
        }
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(readOutsideSpecifiedRegions(
                read, mConfig.SpecificChrRegions.Regions, mConfig.SpecificChrRegions.Chromosomes, mConfig.SpecificRegionsFilterType))
        {
            return;
        }

        ++mStats.TotalReads;
        ++mPartitionRecordCount;

        if(mConfig.RunChecks)
            mBamWriter.registerRead(read);

        int readStart = read.getAlignmentStart();

        while(mCurrentPartition != null && readStart > mCurrentPartition.end())
        {
            onPartitionComplete(true);

            if(mCurrentPartition == null)
            {
                mBamWriter.writeFragment(new Fragment(read));
                return;
            }
        }

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName())) // debugging only
        {
            MD_LOGGER.debug("specific read: {}", readToString(read));
        }

        if(mConfig.UnmapRegions.enabled())
        {
            mConfig.UnmapRegions.checkTransformRead(read, mUnmapRegions);

            if(read.getSupplementaryAlignmentFlag() && read.getReadUnmappedFlag())
                return; // drop unmapped supplementaries

            if(read.getReadUnmappedFlag() && read.getMateUnmappedFlag())
            {
                mBamWriter.writeFragment(new Fragment(read));
                return;
            }
        }

        try
        {
            if(!mReadPositions.processRead(read))
            {
                ++mStats.Incomplete;

                String basePartition = Fragment.getBasePartition(read, mConfig.PartitionSize);

                if(basePartition == null)
                {
                    // mate or supp is on a non-human chromsome, meaning it won't be retrieved - so write this immediately
                    mBamWriter.writeRead(read, FragmentStatus.UNSET);
                    return;
                }

                processIncompleteRead(read, basePartition);
            }
        }
        catch(Exception e)
        {
            MD_LOGGER.error("read({}) exception: {}", readToString(read), e.toString());
            e.printStackTrace();
            System.exit(1);
        }

        if(!mConfig.NoMateCigar && read.getReadPairedFlag() && !read.getMateUnmappedFlag() && !read.getSupplementaryAlignmentFlag()
                && !read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            ++mStats.MissingMateCigar;
        }
    }

    private void processIncompleteRead(final SAMRecord read, final String basePartition)
    {
        if(basePartition.equals(mCurrentStrPartition))
        {
            PartitionResults partitionResults = mCurrentPartitionData.processIncompleteFragment(read);

            if(partitionResults != null)
            {
                if(partitionResults.umiGroups() != null || partitionResults.resolvedFragments() != null)
                {
                    if(partitionResults.umiGroups() != null)
                        partitionResults.umiGroups().forEach(x -> processDuplicateGroup(x));

                    if(partitionResults.resolvedFragments() != null)
                        mBamWriter.writeFragments(partitionResults.resolvedFragments(), true);
                }
                else if(partitionResults.fragmentStatus() != null && partitionResults.fragmentStatus().isResolved())
                {
                    mBamWriter.writeRead(read, partitionResults.fragmentStatus());
                }
            }
        }
        else
        {
            ++mStats.InterPartition;

            // cache this read and send through as groups when the partition is complete
            List<SAMRecord> pendingFragments = mPendingIncompleteReads.get(basePartition);

            if(pendingFragments == null)
            {
                pendingFragments = Lists.newArrayList();
                mPendingIncompleteReads.put(basePartition, pendingFragments);
            }

            pendingFragments.add(read);
        }
    }

    private void processPendingIncompletes()
    {
        if(mPendingIncompleteReads.isEmpty())
            return;

        if(mPendingIncompleteReads.size() > 100)
        {
            MD_LOGGER.debug("partition({}:{}) processing {} pending incomplete fragments",
                    mRegion.Chromosome, mCurrentPartition, mPendingIncompleteReads.values().stream().mapToInt(x -> x.size()).sum());
        }

        mPcPendingIncompletes.resume();

        for(Map.Entry<String,List<SAMRecord>> entry : mPendingIncompleteReads.entrySet())
        {
            String basePartition = entry.getKey();
            List<SAMRecord> reads = entry.getValue();

            PartitionData partitionData = mPartitionDataStore.getOrCreatePartitionData(basePartition);

            PartitionResults partitionResults = partitionData.processIncompleteFragments(reads);

            if(partitionResults.umiGroups() != null)
                partitionResults.umiGroups().forEach(x -> processDuplicateGroup(x));

            if(partitionResults.resolvedFragments() != null)
                mBamWriter.writeFragments(partitionResults.resolvedFragments(), true);
        }

        mPendingIncompleteReads.clear();

        mPcPendingIncompletes.pause();
    }

    private void processDuplicateGroup(final DuplicateGroup duplicateGroup)
    {
        // form consensus reads for any complete read leg groups and write reads
        List<SAMRecord> completeReads = duplicateGroup.popCompletedReads(mConsensusReads, false);
        mBamWriter.writeDuplicateGroup(duplicateGroup, completeReads);
    }

    public void accept(final List<Fragment> positionFragments)
    {
        if(positionFragments.isEmpty())
            return;

        mPcAcceptPositions.resume();
        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<CandidateDuplicates> candidateDuplicatesList = Lists.newArrayList();
        List<List<Fragment>> positionDuplicateGroups = Lists.newArrayList();

        int posFragmentCount = positionFragments.size();
        boolean logDetails = mConfig.PerfDebug && posFragmentCount > 10000;
        long startTimeMs = logDetails ? System.currentTimeMillis() : 0;

        int position = positionFragments.get(0).initialPosition();

        boolean inExcludedRegion = mExcludedRegion != null && positionFragments.stream()
                .anyMatch(x -> x.reads().stream().anyMatch(y -> FragmentUtils.overlapsExcludedRegion(mExcludedRegion, y)));

        findDuplicateFragments(positionFragments, resolvedFragments, positionDuplicateGroups, candidateDuplicatesList, mConfig.UMIs.Enabled);

        List<Fragment> singleFragments = mConfig.UMIs.Enabled && !inExcludedRegion ?
                resolvedFragments.stream().filter(x -> x.status() == FragmentStatus.NONE).collect(Collectors.toList()) : Collections.EMPTY_LIST;

        List<DuplicateGroup> duplicateGroups = mDuplicateGroupBuilder.processDuplicateGroups(
                positionDuplicateGroups, true, singleFragments, inExcludedRegion);

        if(logDetails)
        {
            double timeTakenSec = (System.currentTimeMillis() - startTimeMs) / 1000.0;

            if(timeTakenSec >= 1.0)
            {
                MD_LOGGER.debug("position({}:{}) fragments({}) resolved({}) candidates({}) processing time({})",
                        mRegion.Chromosome, position, posFragmentCount, resolvedFragments.size(),
                        candidateDuplicatesList.stream().mapToInt(x -> x.fragmentCount()).sum(),
                        format("%.1fs", timeTakenSec));
            }
        }

        startTimeMs = logDetails ? System.currentTimeMillis() : 0;

        mCurrentPartitionData.processPrimaryFragments(resolvedFragments, candidateDuplicatesList, duplicateGroups);

        if(logDetails)
        {
            double timeTakenSec = (System.currentTimeMillis() - startTimeMs) / 1000.0;

            if(timeTakenSec >= 1.0)
            {
                MD_LOGGER.debug("position({}:{}) fragments({}) partition processing time({})",
                        mRegion.Chromosome, position, posFragmentCount, format("%.1fs", timeTakenSec));
            }
        }

        if(duplicateGroups != null)
            duplicateGroups.forEach(x -> processDuplicateGroup(x));

        if(!resolvedFragments.isEmpty())
        {
            mBamWriter.writeFragments(resolvedFragments, true);

            mStats.LocalComplete += (int)resolvedFragments.stream().filter(x -> x.allReadsPresent()).count();
        }

        mPcAcceptPositions.pause();
    }

    @VisibleForTesting
    public void setExcludedRegion(final ChrBaseRegion excludedRegion)
    {
        if(excludedRegion.overlaps(mRegion) && mCurrentPartition.overlaps(excludedRegion))
        {
            mExcludedRegion = excludedRegion;
            mCurrentPartitionData.setExcludedRegion(excludedRegion);
        }
        else
        {
            mExcludedRegion = null;
        }
    }

    private void setUnmappedRegions()
    {
        mUnmapRegions.clear();

        List<BaseRegion> chrRegions = mConfig.UnmapRegions.getRegions(mRegion.Chromosome);

        if(chrRegions != null)
            chrRegions.stream().filter(x -> x.overlaps(mCurrentPartition)).forEach(x -> mUnmapRegions.add(x));
    }

    private void perfCountersStart()
    {
        if(mConfig.PerfDebug)
            mPcTotal.start(format("%s:%s", mRegion.Chromosome, mCurrentPartition));
        else
            mPcTotal.start();

        mPcAcceptPositions.startPaused();
        mPcPendingIncompletes.startPaused();
    }

    private void perfCountersStop()
    {
        mPcTotal.stop();
        mPcAcceptPositions.stop();
        mPcPendingIncompletes.stop();
    }

    @VisibleForTesting
    public void processRead(final SAMRecord read) { processSamRecord(read); }

    @VisibleForTesting
    public void flushReadPositions() { mReadPositions.evictAll(); }

    @VisibleForTesting
    public void flushPendingIncompletes() { processPendingIncompletes(); }

    @VisibleForTesting
    public void onChromosomeComplete() { onPartitionComplete(false); }

    @VisibleForTesting
    public PartitionDataStore partitionDataStore() { return mPartitionDataStore; }

    @VisibleForTesting
    public ConsensusReads consensusReads() { return mConsensusReads; }
}
