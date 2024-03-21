package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.PartitionTask;

import htsjdk.io.HtsPath;
import htsjdk.samtools.*;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;


public class BqrThread extends Thread
{
    private final IndexedFastaSequenceFile mRefGenome;
    private final SageConfig mConfig;
    private final SamReader mBamReader;
    private final Queue<PartitionTask> mRegions;
    private final BaseQualityResults mResults;

    private final BqrRegionReader mRegionCounter; // will be reused for each region
    //private final htsjdk.samtools.InputResource.Type InputResource = htsjdk.samtools.InputResource.Type.HTSGET;

    public BqrThread(
            final SageConfig config, final IndexedFastaSequenceFile refGenome, final String bamFile,
            final Queue<PartitionTask> regions, final BaseQualityResults results, final BqrRecordWriter recordWriter)  {
        mRefGenome = refGenome;
        mConfig = config;
        mRegions = regions;
        mResults = results;

        HtsPath path = new HtsPath(bamFile);

        mBamReader = SamReaderFactory.makeDefault()
                .validationStringency(mConfig.BamStringency)
                .referenceSource(new ReferenceSource(mRefGenome))
                .open(SamInputResource.of(path.getURI()));

        mRegionCounter = new BqrRegionReader(mConfig, mBamReader, mRefGenome, mResults, recordWriter);

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                PartitionTask partition = mRegions.remove();

                mRegionCounter.initialise(partition.Partition);

                if(partition.TaskId > 0 && (partition.TaskId % 100) == 0)
                {
                    SG_LOGGER.debug("base-qual regions assigned({}) remaining({})",
                            partition.TaskId, mRegions.size());
                }

                mRegionCounter.run();
            }
            catch(NoSuchElementException e)
            {
                SG_LOGGER.trace("all tasks complete");
                break;
            }
        }
    }
}
