package com.hartwig.hmftools.svprep.depth;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.ALLELE_FRACTION_DESC;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REF_READPAIR_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REF_READPAIR_COVERAGE_DESC;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REF_READ_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REF_READ_COVERAGE_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.ExcludedRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class DepthAnnotator
{
    private final DepthConfig mConfig;
    private final Map<String,Integer> mSampleVcfGenotypeIds;

    private final Map<String,List<VariantContext>> mChrVariantMap;
    private final List<ChrBaseRegion> mExcludedRegions;

    public DepthAnnotator(final CommandLine cmd)
    {
        mConfig = new DepthConfig(cmd);
        mChrVariantMap = Maps.newHashMap();
        mSampleVcfGenotypeIds = Maps.newHashMap();
        mExcludedRegions = ExcludedRegions.getPolyGRegions(mConfig.RefGenVersion);
    }

    public void run()
    {
        if(mConfig.InputVcf == null || mConfig.OutputVcf == null)
        {
            SV_LOGGER.error("missing VCF");
            System.exit(1);
        }

        if(mConfig.BamFiles.size() != mConfig.Samples.size())
        {
            SV_LOGGER.error("inconsistent samples and BAM files");
            System.exit(1);
        }

        SV_LOGGER.info("SvPrep depth annotation for samples: {}", mConfig.Samples);

        long startTimeMs = System.currentTimeMillis();

        final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                mConfig.InputVcf, new VCFCodec(), false);

        VCFHeader vcfHeader = (VCFHeader)reader.getHeader();

        if(!establishGenotypeIds(vcfHeader))
        {
            System.exit(1);
        }

        int vcfCount = 0;

        try
        {
            String currentChromosome = "";
            List<VariantContext> variantsList = null;

            for(VariantContext variantContext : reader.iterator())
            {
                ++vcfCount;

                if(excludeVariant(variantContext))
                    continue;

                String chromosome = variantContext.getContig();

                if(!chromosome.equals(currentChromosome))
                {
                    variantsList = Lists.newArrayList();
                    mChrVariantMap.put(chromosome, variantsList);
                    currentChromosome = chromosome;
                }

                // create new instances of each variant to set depth values against
                VariantContext newVariant = new VariantContextBuilder(variantContext)
                        .genotypes(variantContext.getGenotypes())
                        .filters(variantContext.getFilters())
                        .make();

                variantsList.add(newVariant);
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("error reading vcf({}): {}", mConfig.InputVcf, e.toString());
            System.exit(1);
        }

        SV_LOGGER.info("loaded {} variants from vcf({})", vcfCount, mConfig.InputVcf);

        if(mChrVariantMap.isEmpty())
        {
            SV_LOGGER.warn("all variants filtered from vcf({})", vcfCount, mConfig.InputVcf);
            return;
        }

        if(mConfig.PerfLogTime > 0)
            analyseVariantDistribution();

        List<DepthTask> depthTasks = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            List<VariantContext> variantsList = mChrVariantMap.get(chrStr);

            if(variantsList == null)
                continue;

            DepthTask depthTask = new DepthTask(chrStr, mConfig, mSampleVcfGenotypeIds);
            depthTask.addVariants(variantsList);
            depthTasks.add(depthTask);
        }

        final List<Callable> callableList = depthTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        // write output VCF
        writeVcf(vcfHeader, depthTasks);

        double timeTakeMins = (System.currentTimeMillis() - startTimeMs) / 60000.0;

        SV_LOGGER.info("SvPrep depth annotation complete, mins({})", format("%.3f", timeTakeMins));

        PerformanceCounter perfCounter = depthTasks.get(0).getPerfCounter();
        for(int i = 1; i < depthTasks.size(); ++i)
        {
            perfCounter.merge(depthTasks.get(i).getPerfCounter());
        }

        perfCounter.logStats();
    }

    private void writeVcf(final VCFHeader header, final List<DepthTask> depthTasks)
    {
        SV_LOGGER.info("writing VCF: {}", mConfig.OutputVcf);

        VariantContextWriter writer = new VariantContextWriterBuilder()
                .setReferenceDictionary(header.getSequenceDictionary())
                .setOutputFile(mConfig.OutputVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();

        if(!header.hasFormatLine(ALLELE_FRACTION))
            header.addMetaDataLine(new VCFFormatHeaderLine(ALLELE_FRACTION, 1, VCFHeaderLineType.Float, ALLELE_FRACTION_DESC));

        if(!header.hasFormatLine(REF_READ_COVERAGE))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(REF_READ_COVERAGE, 1, VCFHeaderLineType.Integer, REF_READ_COVERAGE_DESC));
            header.addMetaDataLine(new VCFInfoHeaderLine(REF_READ_COVERAGE, 1, VCFHeaderLineType.Integer, REF_READ_COVERAGE_DESC));
        }
        else if(mConfig.VcfTagPrefix != null)
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    mConfig.getVcfTag(REF_READ_COVERAGE), 1, VCFHeaderLineType.Integer, REF_READ_COVERAGE_DESC));
            header.addMetaDataLine(new VCFInfoHeaderLine(
                    mConfig.getVcfTag(REF_READ_COVERAGE), 1, VCFHeaderLineType.Integer, REF_READ_COVERAGE_DESC));
        }

        if(!header.hasFormatLine(REF_READPAIR_COVERAGE))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(REF_READPAIR_COVERAGE, 1, VCFHeaderLineType.Integer, REF_READPAIR_COVERAGE_DESC));
            header.addMetaDataLine(new VCFInfoHeaderLine(REF_READPAIR_COVERAGE, 1, VCFHeaderLineType.Integer, REF_READPAIR_COVERAGE_DESC));
        }
        else if(mConfig.VcfTagPrefix != null)
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    mConfig.getVcfTag(REF_READPAIR_COVERAGE), 1, VCFHeaderLineType.Integer, REF_READPAIR_COVERAGE_DESC));

            header.addMetaDataLine(new VCFInfoHeaderLine(
                    mConfig.getVcfTag(REF_READPAIR_COVERAGE), 1, VCFHeaderLineType.Integer, REF_READPAIR_COVERAGE_DESC));
        }

        writer.writeHeader(header);

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());
            DepthTask depthTask = depthTasks.stream().filter(x -> x.chromosome().equals(chrStr)).findFirst().orElse(null);

            if(depthTask == null)
                continue;

            depthTask.variants().forEach(x -> writer.add(x));
        }

        writer.close();
    }

    private boolean establishGenotypeIds(final VCFHeader header)
    {
        List<String> vcfSampleNames = header.getGenotypeSamples();

        for(int s = 0; s < mConfig.Samples.size(); ++s)
        {
            String sampleId = mConfig.Samples.get(s);
            boolean found = false;

            for(int i = 0; i < vcfSampleNames.size(); ++i)
            {
                String vcfSampleName = vcfSampleNames.get(i);

                if(sampleId.equals(vcfSampleName))
                {
                    mSampleVcfGenotypeIds.put(sampleId, i);
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                SV_LOGGER.error("sample({}) not found in VCF genotype names", sampleId);
                return false;
            }
        }

        return true;
    }

    private void analyseVariantDistribution()
    {
        Map<Integer,Integer> groupFrequencies = Maps.newHashMap();

        int totalGroups = 0;
        int totalVariants = 0;
        long estimatedReads = 0;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            List<VariantContext> variants = mChrVariantMap.get(chrStr);

            if(variants == null)
                continue;

            int soloVariants = 0;
            int groupCount = 0;

            int index = 0;
            while(index < variants.size())
            {
                VariantContext variant = variants.get(index);

                int variantCount = 1;

                int posStart = variant.getStart();
                int posEnd = posStart;
                int nextIndex = index + 1;
                while(nextIndex < variants.size())
                {
                    VariantContext nextVariant = variants.get(nextIndex);

                    if(nextVariant.getStart() - posEnd > mConfig.ProximityDistance)
                        break;

                    posEnd = nextVariant.getStart();
                    ++variantCount;
                    ++nextIndex;
                }

                Integer countFrequency = groupFrequencies.get(variantCount);
                groupFrequencies.put(variantCount, countFrequency != null ? countFrequency + 1 : 1);
                ++groupCount;
                ++totalGroups;
                totalVariants += variantCount;

                if(variantCount == 1)
                    ++soloVariants;

                estimatedReads += (posEnd - posStart + 2000);

                index += variantCount;
            }

            SV_LOGGER.debug("chr({}) variants({}) group({}) soloVariants({} pct={})",
                    chrStr, variants.size(), groupCount, soloVariants, format("%.3f", soloVariants / (double)variants.size()));
        }

        int largeGroupCount = 0;
        int largeVariantsCount = 0;

        for(Map.Entry<Integer,Integer> entry : groupFrequencies.entrySet())
        {
            if(entry.getKey() >= 25)
            {
                ++largeGroupCount;
                largeVariantsCount += entry.getKey() * entry.getValue();
            }
            else
            {
                SV_LOGGER.debug("group count({}) frequency({}) total variants({})",
                        entry.getKey(), entry.getValue(), entry.getKey() * entry.getValue());
            }
        }

        SV_LOGGER.debug("large group count({}) total variants({})", largeGroupCount, largeVariantsCount);

        SV_LOGGER.debug("total variants({}) groups({}) estimated reads({})", totalVariants, totalGroups, estimatedReads);
    }

    private boolean excludeVariant(final VariantContext variant)
    {
        if(mExcludedRegions.stream().anyMatch(x -> x.containsPosition(variant.getStart())))
            return true;

        if(mConfig.SpecificRegions.isEmpty())
            return false;

        return mConfig.SpecificRegions.stream().noneMatch(x -> x.containsPosition(variant.getContig(), variant.getStart()));
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        DepthConfig.addOptions(options);
        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        DepthAnnotator bamSvSlicer = new DepthAnnotator(cmd);
        bamSvSlicer.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
