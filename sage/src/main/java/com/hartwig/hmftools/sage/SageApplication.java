package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.coverage.GeneCoverage.populateCoverageBuckets;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.Sets;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.coverage.GeneDepthFile;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.ChromosomePipeline;
import com.hartwig.hmftools.sage.quality.BaseQualityRecalibration;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.vcf.VariantContextFactory;
import com.hartwig.hmftools.sage.vcf.VariantFile;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class SageApplication implements AutoCloseable
{
    private final SageConfig mConfig;
    private final ReferenceData mRefData;

    private final ExecutorService mExecutorService;
    private final PhaseSetCounter mPhaseSetCounter;

    private final VariantVCF mVcfFile;
    private final VariantFile mVariantFile;

    private SageApplication(final CommandLine cmd)
    {
        final VersionInfo version = new VersionInfo("sage.version");
        SG_LOGGER.info("Sage version: {}", version.version());

        mConfig = new SageConfig(false, version.version(), cmd);

        if(!mConfig.isValid())
        {
            System.exit(1);
            SG_LOGGER.error("invalid config, exiting");
        }

        mRefData = new ReferenceData(mConfig, cmd);

        if(!mRefData.load())
        {
            System.exit(1);
            SG_LOGGER.error("invalid reference data, exiting");
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("SAGE-%d").build();
        mExecutorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        mPhaseSetCounter = new PhaseSetCounter();

        mVcfFile = new VariantVCF(mRefData.RefGenome, mConfig);

        if(mConfig.WriteCsv && !mConfig.TumorIds.isEmpty())
        {
            mVariantFile = new VariantFile(mConfig.TumorIds.get(0), mConfig.SampleDataDir);
        }
        else
        {
            mVariantFile = null;
        }

        SG_LOGGER.info("writing to file: {}", mConfig.OutputFile);
    }

    private void run() throws IOException
    {
        long startTime = System.currentTimeMillis();
        final Coverage coverage = createCoverage();

        BaseQualityRecalibration baseQualityRecalibration = new BaseQualityRecalibration(mConfig, mExecutorService, mRefData.RefGenome);
        baseQualityRecalibration.produceRecalibrationMap();
        final Map<String,QualityRecalibrationMap> recalibrationMap = baseQualityRecalibration.getSampleRecalibrationMap();

        final SAMSequenceDictionary dictionary = dictionary();
        for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            final String chromosome = samSequenceRecord.getSequenceName();

            if(mConfig.SpecificChromosomes.isEmpty() || mConfig.SpecificChromosomes.contains(chromosome))
            {
                if(!HumanChromosome.contains(chromosome) && !MitochondrialChromosome.contains(chromosome))
                    continue;

                try
                {
                    final ChromosomePipeline pipeline = new ChromosomePipeline(
                            chromosome, mConfig, mExecutorService, mRefData, recalibrationMap, coverage, mPhaseSetCounter, this::writeVariant);

                    pipeline.process();
                }
                catch(Exception e)
                {
                    SG_LOGGER.error("chromosome({}) failed to execute pipeline tasks: {}", chromosome, e.toString());
                    e.printStackTrace();
                }

                System.gc();
            }
        }

        for(String sample : coverage.samples())
        {
            String filename = mConfig.geneCoverageFile(sample);
            GeneDepthFile.write(filename, coverage.depth(sample));
        }

        long endTime = System.currentTimeMillis();
        double runTime = (endTime - startTime) / 1000.0;

        SG_LOGGER.info("Sage complete, run time({}s)", String.format("%.2f", runTime));
    }

    private SAMSequenceDictionary dictionary() throws IOException
    {
        final String bam = mConfig.ReferenceBams.isEmpty() ? mConfig.TumorBams.get(0) : mConfig.ReferenceBams.get(0);

        SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(mConfig.Stringency)
                .referenceSource(new ReferenceSource(mRefData.RefGenome))
                .open(new File(bam));

        SAMSequenceDictionary dictionary = tumorReader.getFileHeader().getSequenceDictionary();

        tumorReader.close();

        return dictionary;
    }

    public void writeVariant(final SageVariant variant)
    {
        mVcfFile.write(VariantContextFactory.create(variant));

        if(mVariantFile != null)
            mVariantFile.writeToFile(variant);
    }

    private Coverage createCoverage()
    {
        populateCoverageBuckets();

        Set<String> samples = Sets.newHashSet();
        if(!mConfig.CoverageBed.isEmpty())
        {
            samples.addAll(mConfig.TumorIds);
        }

        return new Coverage(samples, mRefData.CoveragePanel.values());
    }

    @Override
    public void close() throws IOException
    {
        mVcfFile.close();

        if(mVariantFile != null)
            mVariantFile.close();

        mRefData.RefGenome.close();
        mExecutorService.shutdown();
    }

    public static void main(final String... args) throws IOException, InterruptedException, ExecutionException
    {
        final Options options = SageConfig.createSageOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            final SageApplication application = new SageApplication(cmd);
            application.run();
            application.close();
        }
        catch(ParseException e)
        {
            SG_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SageApplication", options);
            System.exit(1);
        }
    }

    public static CommandLine createCommandLine(final String[] args, final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}