package com.hartwig.hmftools.isofox.unmapped;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaExpressionMatrix.EXPRESSION_SCOPE_GENE;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.START_STR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.startEndStr;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.UNMAPPED_READS;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.SPLICE_TYPE_ACCEPTOR;
import static com.hartwig.hmftools.isofox.unmapped.UnmappedRead.UMR_NO_MATE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.rna.RnaExpressionMatrix;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class UmrCohortAnalyser
{
    private final CohortConfig mConfig;

    // other config
    private final RnaExpressionMatrix mGeneExpression;

    // map of chromosomes to a map of splice-boundary keys to a map of samples to a list of unmapped reads
    private final Map<String,Map<String,Map<String,List<UnmappedRead>>>> mUnmappedReads;

    private final LineElementMatcher mLineElementMatcher;
    private final boolean mCombineFrequencies;

    private final BufferedWriter mWriter;

    private static final String LINX_DIRECTORY = "linx_dir";
    private static final String GENE_EXPRESSION_FILE = "gene_expression_file";
    private static final String COMBINE_FREQUENCIES = "combine_frequencies";
    private static final String LINX_DIR = "linx_dir";
    private static final String SV_VCF_FILE = "sv_vcf";

    public UmrCohortAnalyser(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mUnmappedReads = Maps.newHashMap();

        mGeneExpression = cmd.hasOption(GENE_EXPRESSION_FILE) ?
                new RnaExpressionMatrix(cmd.getOptionValue(GENE_EXPRESSION_FILE), EXPRESSION_SCOPE_GENE) : null;

        mLineElementMatcher = cmd.hasOption(LINX_DIR) && cmd.hasOption(SV_VCF_FILE) ?
                new LineElementMatcher(cmd.getOptionValue(LINX_DIR), cmd.getOptionValue(SV_VCF_FILE)) : null;

        mCombineFrequencies = cmd.hasOption(COMBINE_FREQUENCIES);

        mWriter = initialiseWriter();
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(LINX_DIRECTORY, true, "Path to Linx files");
        options.addOption(GENE_EXPRESSION_FILE, true, "Gene expression file for cohort");
        options.addOption(COMBINE_FREQUENCIES, false, "Determine cohort frequencies for slice candidates");
        options.addOption(LINX_DIR, true, "Linx data directory");
        options.addOption(SV_VCF_FILE, true, "Structural variant VCF");
    }

    public void processSampleFiles()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, UNMAPPED_READS, filenames))
            return;

        int nextLog = 100000;

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path umrFile = filenames.get(i);

            loadFile(sampleId, umrFile);

            if(mCombineFrequencies)
            {
                int totalUmrCount = mUnmappedReads.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();

                if(totalUmrCount >= nextLog)
                {
                    ISF_LOGGER.debug("cached unmapped-read count({})", totalUmrCount);
                    nextLog += 100000;
                }
            }
            else
            {
                mLineElementMatcher.findMatches(sampleId, mUnmappedReads);
                writeUnmappedReads();
                mUnmappedReads.clear();
            }

            if(i > 0 && (i % 100) == 0)
            {
                ISF_LOGGER.debug("processed {} samples", i);
            }
        }

        int totalUmrCount = mUnmappedReads.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();
        ISF_LOGGER.info("loaded {} unmapped-read records", totalUmrCount);

        ISF_LOGGER.info("writing cohort unmapped-read");

        if(mCombineFrequencies)
            writeUnmappedReads();

        closeBufferedWriter(mWriter);
    }

    private void loadFile(final String sampleId, final Path filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);

            lines.remove(0);

            int geneIdIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int geneNameIndex = fieldsIndexMap.get(FLD_GENE_NAME);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int readIndex = fieldsIndexMap.get("ReadId");
            int posStartIndex = fieldsIndexMap.get("ReadStart");
            int posEndIndex = fieldsIndexMap.get("ReadEnd");
            int spliceTypeIndex = fieldsIndexMap.get("SpliceType");
            int scLengthIndex = fieldsIndexMap.get("SoftClipLength");
            int scSideIndex = fieldsIndexMap.get("SoftClipSide");
            int abqIndex = fieldsIndexMap.get("AvgBaseQual");
            int transIndex = fieldsIndexMap.get("TransName");
            int exonRankIndex = fieldsIndexMap.get("ExonRank");
            int exonBoundaryIndex = fieldsIndexMap.get("ExonBoundary");
            int exonDistIndex = fieldsIndexMap.get("ExonDistance");
            int scBasesIndex = fieldsIndexMap.get("SoftClipBases");
            Integer mateIndex = fieldsIndexMap.get("MateCoords");
            Integer cohortFreqIndex = fieldsIndexMap.get("CohortFreq");
            Integer matchesSuppIndex = fieldsIndexMap.get("MatchesChimeric");

            for(String data : lines)
            {
                final String[] values = data.split(DELIMITER);

                UnmappedRead umRead = new UnmappedRead(
                        values[readIndex],
                        new ChrBaseRegion(values[chrIndex], Integer.parseInt(values[posStartIndex]), Integer.parseInt(values[posEndIndex])),
                        Integer.parseInt(values[scLengthIndex]), Integer.parseInt(values[scSideIndex]),
                        Double.parseDouble(values[abqIndex]), values[geneIdIndex], values[geneNameIndex], values[transIndex],
                        Integer.parseInt(values[exonRankIndex]), Integer.parseInt(values[exonBoundaryIndex]),
                        Integer.parseInt(values[exonDistIndex]), values[spliceTypeIndex], values[scBasesIndex],
                        mateIndex != null ? values[mateIndex] : "",
                        cohortFreqIndex != null ? Integer.parseInt(values[cohortFreqIndex]) : 0,
                        matchesSuppIndex != null ? Boolean.parseBoolean(values[matchesSuppIndex]) : false);

                // could filter these in the BAM reading process
                if(umRead.ExonRank == 1 && umRead.SpliceType.equals(SPLICE_TYPE_ACCEPTOR))
                    continue;

                addUnmappedRead(sampleId, umRead);
            }

            ISF_LOGGER.debug("sample({}) loaded {} unmapped-read records", sampleId, lines.size());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load unmapped-read file({}): {}", filename.toString(), e.toString());
        }
    }

    private void addUnmappedRead(final String sampleId, final UnmappedRead umRead)
    {
        Map<String,Map<String,List<UnmappedRead>>> chrUmrs = mUnmappedReads.get(umRead.ReadRegion.Chromosome);

        if(chrUmrs == null)
        {
            chrUmrs = Maps.newHashMap();
            mUnmappedReads.put(umRead.ReadRegion.Chromosome, chrUmrs);
        }

        String umrKey = umRead.positionKey();
        Map<String,List<UnmappedRead>> umrKeyList = chrUmrs.get(umrKey);

        if(umrKeyList == null)
        {
            umrKeyList = Maps.newHashMap();
            chrUmrs.put(umrKey, umrKeyList);
        }

        List<UnmappedRead> sampleReads = umrKeyList.get(sampleId);

        if(sampleReads == null)
        {
            sampleReads = Lists.newArrayList();
            umrKeyList.put(sampleId, sampleReads);
        }

        sampleReads.add(umRead);
    }

    private BufferedWriter initialiseWriter()
    {
        final String outputFile = mConfig.formCohortFilename("combined_unmapped_reads.csv");

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);
            writer.write("SampleId,FragmentCount,UnpairedCount,Chromosome,GeneName,TransName,ExonRank,SpliceType");
            writer.write(",ExonBoundary,ExonDistance,SoftClipSide,AvgBaseQual");
            writer.write(",SoftClipBases,GeneTPM,HasChimericMatch");

            if(mCombineFrequencies)
                writer.write(",CohortFrequency");

            if(mLineElementMatcher != null)
                writer.write(",SvLinxMatches");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to initialise cohort unmapped reads file({}): {}", outputFile, e.toString());
            return null;
        }
    }

    private void writeUnmappedReads()
    {
        if(mWriter == null)
            return;

        try
        {
            for(Map.Entry<String,Map<String,Map<String,List<UnmappedRead>>>> chrEntry : mUnmappedReads.entrySet())
            {
                String chromosome = chrEntry.getKey();

                for(Map<String,List<UnmappedRead>> umrSampleMap : chrEntry.getValue().values())
                {
                    int sampleCount = umrSampleMap.size();

                    for(Map.Entry<String,List<UnmappedRead>> sampleEntry : umrSampleMap.entrySet())
                    {
                        String sampleId = sampleEntry.getKey();

                        List<UnmappedRead> umReads = sampleEntry.getValue();
                        UnmappedRead firstRead = umReads.get(0);

                        String matchedSVsInfo = "";
                        if(mLineElementMatcher != null)
                        {
                            List<StructuralVariant> matchedSVs = mLineElementMatcher.getUmrMatch(firstRead.positionKey());

                            if(matchedSVs == null)
                                continue;

                            StringJoiner sj = new StringJoiner(ITEM_DELIM);

                            // VcfId:Type:Position:OtherPos
                            matchedSVs.forEach(x -> sj.add(String.format("%s:%s:%d:%d",
                                    x.id(), x.type(), x.position(true), x.end() != null ? x.position(false) : 0)));

                            matchedSVsInfo = sj.toString();
                        }

                        // de-dup reads by fragment and across genes sharing the same exon boundary using readId
                        Set<String> readIds = Sets.newHashSet();
                        Set<String> genes = Sets.newHashSet();

                        double geneTpm = 0;
                        boolean hasSuppMatch = false;
                        int unpairedCount = 0;

                        for(UnmappedRead umRead : umReads)
                        {
                            if(!genes.contains(umRead.GeneName))
                            {
                                genes.add(umRead.GeneName);

                                if(mGeneExpression != null)
                                    geneTpm += mGeneExpression.getExpression(umRead.GeneId, sampleId);
                            }

                            readIds.add(umRead.ReadId);
                            hasSuppMatch |= umRead.MatchesChimeric;

                            if(umRead.MateCoords.equals(UMR_NO_MATE))
                                ++unpairedCount;
                        }

                        StringJoiner genesStr = new StringJoiner(ITEM_DELIM);
                        genes.forEach(x -> genesStr.add(x));

                        int fragmentCount = readIds.size();

                        double avgBaseQual = umReads.stream().mapToDouble(x -> x.AvgBaseQual).sum() / umReads.size();

                        mWriter.write(String.format("%s,%d,%d,%s,%s,%s,%d,%s",
                                sampleId, fragmentCount, unpairedCount, chromosome, genesStr.toString(), firstRead.TransName,
                                firstRead.ExonRank, firstRead.SpliceType));

                        mWriter.write(String.format(",%d,%d,%s,%.1f,%s,%4.3e,%s",
                                firstRead.ExonBoundary, firstRead.ExonDistance, startEndStr(firstRead.ScSide),
                                avgBaseQual, firstRead.ScBases, geneTpm, hasSuppMatch));

                        if(mCombineFrequencies)
                        {
                            mWriter.write(String.format(",%d", sampleCount));
                        }

                        if(mLineElementMatcher != null)
                        {
                            mWriter.write(String.format(",%s", matchedSVsInfo));
                        }

                        mWriter.newLine();
                    }
                }
            }

        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write cohort unmapped reads file: {}", e.toString());
        }
    }
}