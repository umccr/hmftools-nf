package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.LOG_COUNT;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.ReadGroup.mergeChimericReadMaps;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.TaskExecutor;

public class FusionFinder
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final List<ReadGroup> mChimericReadGroups;
    private final Map<String,ReadGroup> mChimericPartialReadGroups;
    private final Set<String> mDuplicateReadIds;

    private final List<FusionTask> mFusionTasks;
    private final FusionWriter mFusionWriter;
    private final FusionGeneFilters mGeneFilters;

    private final PerformanceCounter mPerfCounter;

    public FusionFinder(final IsofoxConfig config, final EnsemblDataCache geneTransCache, final FusionWriter fusionWriter)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mChimericPartialReadGroups = Maps.newHashMap();
        mChimericReadGroups = Lists.newArrayList();
        mDuplicateReadIds = Sets.newHashSet();
        mFusionTasks = Lists.newArrayList();

        mGeneFilters = new FusionGeneFilters(config, geneTransCache);

        mPerfCounter = new PerformanceCounter("Fusions");
        mFusionWriter = fusionWriter;
    }

    public FusionGeneFilters getGeneFilters() { return mGeneFilters; }

    public void addChimericReads(final Map<String,ReadGroup> partialReadGroups)
    {
        mergeChimericReadMaps(mChimericPartialReadGroups, mChimericReadGroups, partialReadGroups);
    }

    public void addDuplicateReadIds(final Set<String> readIds)
    {
        mergeDuplicateReadIds(mDuplicateReadIds, readIds);
    }

    public static void addChimericReads(final Map<String,ReadGroup> chimericReadMap, final ReadRecord read)
    {
        ReadGroup chimericReads = chimericReadMap.get(read.Id);
        if (chimericReads == null)
            chimericReadMap.put(read.Id, new ReadGroup(read));
        else
            chimericReads.Reads.add(read);
    }

    private static final String LOG_READ_ID = "";
    // private static final String LOG_READ_ID = "NB500901:18:HTYNHBGX2:2:23308:18394:18413";

    public void findFusions()
    {
        // convert any set of valid reads into a fragment, and then process these in groups by chromosomal pair
        ISF_LOGGER.info("processing {} chimeric read groups", mChimericReadGroups.size());

        mPerfCounter.start();

        int invalidFragments = 0;
        int hasMissingReads = 0;
        int hasExcesssReads = 0;
        int duplicates = 0;
        int skipped = 0;
        int fragments = 0;
        int missingSuppReads = 0;

        int partialDups = 0;
        int partialSkipped = 0;

        int readGroupCount = 0;
        int nextLog = LOG_COUNT;

        final Map<String,List<FusionFragment>> chrPairFragments = Maps.newHashMap();

        mChimericReadGroups.addAll(mChimericPartialReadGroups.values());

        for(ReadGroup readGroup : mChimericReadGroups)
        {
            ++readGroupCount;

            if(readGroupCount >= nextLog)
            {
                nextLog += LOG_COUNT;
                ISF_LOGGER.info("processed {} chimeric read groups", readGroupCount);
            }

            boolean isComplete = readGroup.isComplete();
            final List<ReadRecord> reads = readGroup.Reads;

            if(reads.get(0).Id.equals(LOG_READ_ID))
            {
                ISF_LOGGER.debug("specific read: {}", reads.get(0));
            }

            if(mDuplicateReadIds.contains(reads.get(0).Id) || reads.stream().anyMatch(x -> x.isDuplicate()))
            {
                ++duplicates;

                if(!isComplete)
                    ++partialDups;

                continue;
            }

            if(reads.stream().anyMatch(x -> mGeneFilters.skipRead(x.mateChromosome(), x.mateStartPosition())))
            {
                ++skipped;

                if(!isComplete)
                    ++partialSkipped;

                continue;
            }

            if(!isComplete)
            {
                if(readGroup.hasSuppAlignment())
                {
                    if(skipMissingReads(reads))
                    {
                        ++skipped;
                        continue;
                    }

                    ++missingSuppReads;
                }

                ++invalidFragments;

                if(reads.size() > 3 || (!readGroup.hasSuppAlignment() && reads.size() > 2))
                    ++hasExcesssReads;
                else
                    ++hasMissingReads;

                mFusionWriter.writeReadData(reads, "INVALID_READ_COUNT");
            }
            else
            {
                if(mChimericPartialReadGroups.containsKey(readGroup.id()))
                {
                    ISF_LOGGER.error("partial read({}) group marked as complete", readGroup.id());

                    if(mChimericReadGroups.stream().anyMatch(x -> x.id().equals(readGroup.id())))
                    {
                        ISF_LOGGER.error("partial read({}) group also in complete list", readGroup.id());
                    }
                }

                FusionFragment fragment = new FusionFragment(readGroup);

                if(fragment.type() == FusionFragmentType.UNKNOWN)
                {
                    ++invalidFragments;
                    mFusionWriter.writeReadData(reads, "INVALID_FRAG");
                    continue;
                }

                ++fragments;

                final String chrPair = formChromosomePair(fragment.chromosomes()[SE_START], fragment.chromosomes()[SE_END]);
                List<FusionFragment> fragmentList = chrPairFragments.get(chrPair);

                if(fragmentList == null)
                {
                    chrPairFragments.put(chrPair, Lists.newArrayList(fragment));
                }
                else
                {
                    fragmentList.add(fragment);
                }
            }
        }

        // allocate fusion pairs evenly amongst threads (if multi-thread)
        mFusionTasks.clear();

        for(int taskId = 0; taskId < max(mConfig.Threads, 1); ++taskId)
        {
            mFusionTasks.add(new FusionTask(String.valueOf(taskId), mConfig, mGeneTransCache, mGeneFilters, mFusionWriter));
        }

        for(List<FusionFragment> chrPairFrags : chrPairFragments.values())
        {
            // allocate the next chr-pair fragment batch to the task with the least
            FusionTask leastAllocated = null;
            int minAllocated = 0;

            for(FusionTask fusionTask : mFusionTasks)
            {
                if(minAllocated == 0 || fusionTask.getFragments().size() < minAllocated)
                {
                    leastAllocated = fusionTask;
                    minAllocated = fusionTask.getFragments().size();
                    if(minAllocated == 0)
                        break;
                }
            }

            leastAllocated.getFragments().addAll(chrPairFrags);
        }

        ISF_LOGGER.info("chimeric groups({} skipped={} dups=({} existing={}) invalid={} miss={} candidates={}) chrPairs({}) tasks({})",
                mChimericReadGroups.size(), skipped, duplicates, mDuplicateReadIds.size(), invalidFragments, missingSuppReads, fragments,
                chrPairFragments.size(), mFusionTasks.size());

        if(!mChimericPartialReadGroups.isEmpty())
        {
            int complete = (int)mChimericPartialReadGroups.values().stream().filter(x -> x.isComplete()).count();

            ISF_LOGGER.info("partial groups({} complete={} dups={} skip={} excess={} miss={})",
                    mChimericPartialReadGroups.size(), complete, partialDups, partialSkipped, hasExcesssReads, hasMissingReads);
        }

        mChimericPartialReadGroups.clear();
        mChimericReadGroups.clear();
        mDuplicateReadIds.clear();

        if(mFusionTasks.isEmpty())
        {
            ISF_LOGGER.warn("no fusion tasks created");
            return;
        }
        else
        {
            final List<Callable> callableList = mFusionTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeChromosomeTask(callableList, mConfig.Threads);
            logPerformanceStats();
        }

        mFusionWriter.close();

        mPerfCounter.stop();

        if(mConfig.Fusions.PerformanceStats)
            mPerfCounter.logStats();

        ISF_LOGGER.info("fusion calling complete");
    }

    private boolean skipMissingReads(final List<ReadRecord> reads)
    {
        for(final ReadRecord read : reads)
        {
            if(read.hasSuppAlignment())
            {
                SupplementaryReadData suppData = SupplementaryReadData.from(read.getSuppAlignment());

                if(mGeneFilters.skipRead(suppData.Chromosome, suppData.Position))
                    return true;

                ISF_LOGGER.debug("read({}) missing supp({})", read, suppData);
            }
        }

        return false;
    }

    public static void mergeDuplicateReadIds(final Set<String> destSet, final Set<String> sourceSet)
    {
        sourceSet.forEach(x -> destSet.add(x));
    }

    private void logPerformanceStats()
    {
        if(!mConfig.Fusions.PerformanceStats)
            return;

        if(!ISF_LOGGER.isDebugEnabled() && mConfig.Functions.size() > 1)
            return;

        final PerformanceCounter perfCounter = mFusionTasks.get(0).getPerfCounter();

        for (int i = 1; i < mFusionTasks.size(); ++i)
        {
            perfCounter.merge(mFusionTasks.get(i).getPerfCounter());
        }

        perfCounter.logStats();
    }

    @VisibleForTesting
    public final Map<String,List<FusionReadData>> getFusionCandidates()
    {
        return mFusionTasks.isEmpty() ? Maps.newHashMap() : mFusionTasks.get(0).getFusionCandidates();
    }

    public final Map<String,List<FusionFragment>> getUnfusedFragments()
    {
        return mFusionTasks.isEmpty() ? Maps.newHashMap() : mFusionTasks.get(0).getUnfusedFragments();
    }

    public void addChimericReads(final List<ReadGroup> readGroups)
    {
        mChimericReadGroups.addAll(readGroups);
    }

    public void clearState()
    {
        mChimericPartialReadGroups.clear();
        mChimericReadGroups.clear();
        mFusionTasks.get(0).clearState();
    }
}
