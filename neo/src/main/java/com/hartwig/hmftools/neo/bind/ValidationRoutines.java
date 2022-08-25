package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.stats.AucCalc;
import com.hartwig.hmftools.common.stats.AucData;
import com.hartwig.hmftools.common.utils.MatrixUtils;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.NotNull;

public class ValidationRoutines
{
    private final TrainConfig mConfig;

    private final Map<String, Map<Integer, List<BindData>>> mAlleleTrainingData;
    private final List<String> mValidationAlleles;

    private final Map<String, Map<Integer, BindCountData>> mAlleleBindCounts; // counts data by peptide length>

    private final PosWeightModel mPosWeightModel;
    private final HlaSequences mHlaSequences;
    private final int mThreads;

    private final BufferedWriter mSummaryWriter;
    private final BufferedWriter mPeptideWriter;

    private static final String TEST_ALLELES_FILE = "test_alleles_file";

    public ValidationRoutines(final CommandLine cmd)
    {
        mConfig = new TrainConfig(cmd);

        mValidationAlleles = Lists.newArrayList();

        if(cmd.hasOption(TEST_ALLELES_FILE))
            mValidationAlleles.addAll(loadDelimitedIdFile(cmd.getOptionValue(TEST_ALLELES_FILE), FLD_ALLELE, DELIM));

        mAlleleTrainingData = Maps.newHashMap();
        mAlleleBindCounts = Maps.newHashMap();

        mHlaSequences = new HlaSequences();
        mHlaSequences.load(cmd.getOptionValue(HLA_DEFINITIONS_FILE));
        mPosWeightModel = new PosWeightModel(mConfig.Constants, mHlaSequences);
        mThreads = parseThreads(cmd);

        mSummaryWriter = initialiseSummaryWriter();
        mPeptideWriter = initialisePeptideWriter();
    }

    public void run()
    {
        if(!loadBindData(mConfig.TrainingDataFile, mConfig.RequiredPeptideLengths, mAlleleTrainingData))
        {
            NE_LOGGER.error("failed to load training data");
            System.exit(1);
        }

        if(mValidationAlleles.isEmpty())
            mAlleleTrainingData.keySet().forEach(x -> mValidationAlleles.add(x));

        NE_LOGGER.info("running validation for {} alleles", mValidationAlleles.size());

        // build bind counts and blend across peptide-lengths - this data isn't affected by leaving out each allele
        buildBindCounts();

        if(mThreads > 1)
        {
            List<LeaveOutAlleleTask> alleleTasks = Lists.newArrayList();

            int threads = min(mThreads, mValidationAlleles.size());

            for(int i = 0; i < threads; ++i)
            {
                alleleTasks.add(new LeaveOutAlleleTask(
                        i, mConfig, mPosWeightModel, mAlleleTrainingData, mAlleleBindCounts, mSummaryWriter, mPeptideWriter));
            }

            int taskIndex = 0;
            for(String allele : mValidationAlleles)
            {
                // if allele isn't in the training set then no need for evaluation
                if(!mAlleleTrainingData.containsKey(allele))
                    continue;

                if(taskIndex >= alleleTasks.size())
                    taskIndex = 0;

                alleleTasks.get(taskIndex).getAlleles().add(allele);
                ++taskIndex;
            }

            final List<Callable> callableList = alleleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, callableList.size());
        }
        else
        {
            LeaveOutAlleleTask alleleTask = new LeaveOutAlleleTask(
                    0, mConfig, mPosWeightModel, mAlleleTrainingData, mAlleleBindCounts, mSummaryWriter, mPeptideWriter);
            alleleTask.getAlleles().addAll(mValidationAlleles);
            alleleTask.call();
        }

        closeBufferedWriter(mSummaryWriter);
        closeBufferedWriter(mPeptideWriter);

        NE_LOGGER.info("validation routine complete");
    }

    private void buildBindCounts()
    {
        // translate binds to counts and weighted counts - ie the step prior to allele-blending
        NE_LOGGER.info("building bind counts for {} alleles", mValidationAlleles.size());

        for(Map.Entry<String, Map<Integer, List<BindData>>> alleleEntry : mAlleleTrainingData.entrySet())
        {
            final String allele = alleleEntry.getKey();

            final Map<Integer, List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            Map<Integer, BindCountData> pepLenCountsMap = Maps.newHashMap(); // counts data by peptide length
            mAlleleBindCounts.put(allele, pepLenCountsMap);

            for(Map.Entry<Integer, List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                List<BindData> bindDataList = pepLenEntry.getValue();

                for(BindData bindData : bindDataList)
                {
                    int peptideLength = bindData.peptideLength();

                    BindCountData bindCounts = pepLenCountsMap.get(bindData.peptideLength());

                    if(bindCounts == null)
                    {
                        bindCounts = new BindCountData(allele, peptideLength);
                        pepLenCountsMap.put(peptideLength, bindCounts);
                    }

                    bindCounts.processBindData(bindData, false);
                }
            }
        }

        // fill in any gaps in alleles or peptide lengths if they are required
        if(!mValidationAlleles.isEmpty())
        {
            mValidationAlleles.stream()
                    .filter(x -> !mAlleleBindCounts.containsKey(x))
                    .forEach(x -> mAlleleBindCounts.put(x, Maps.newHashMap()));
        }

        for(Map.Entry<String, Map<Integer, BindCountData>> alleleEntry : mAlleleBindCounts.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer, BindCountData> pepLenBindCountsMap = alleleEntry.getValue();

            if(!mConfig.RequiredPeptideLengths.isEmpty())
            {
                mConfig.RequiredPeptideLengths.stream()
                        .filter(x -> !pepLenBindCountsMap.containsKey(x))
                        .forEach(x -> pepLenBindCountsMap.put(x, new BindCountData(allele, x)));
            }
        }

        for(Map.Entry<String, Map<Integer, BindCountData>> alleleEntry : mAlleleBindCounts.entrySet())
        {
            final Map<Integer, BindCountData> pepLenBindCountsMap = alleleEntry.getValue();

            final List<BindCountData> pepLenBindCounts = pepLenBindCountsMap.values().stream().collect(Collectors.toList());

            for(BindCountData bindCounts : pepLenBindCounts)
            {
                mPosWeightModel.buildWeightedCounts(bindCounts, pepLenBindCounts);
            }

            mPosWeightModel.buildPositionAdjustedTotals(pepLenBindCounts);
        }
    }

    private class LeaveOutAlleleTask implements Callable
    {
        private final int mTaskId;
        private final TrainConfig mConfig;
        private final List<String> mAlleles;

        private final Map<String,Map<Integer,List<BindData>>> mAlleleTrainingData;

        private final Map<String,Map<Integer,BindCountData>> mAlleleBindCounts; // counts data by peptide length>

        private final PosWeightModel mPosWeightModel;
        private final RandomPeptideDistribution mRandomDistribution;

        private final BufferedWriter mSummaryWriter;
        private final BufferedWriter mPeptideWriter;

        public LeaveOutAlleleTask(
                int taskId, final TrainConfig config, final PosWeightModel posWeightModel,
                final Map<String,Map<Integer,List<BindData>>> alleleTrainingData,
                final Map<String, Map<Integer, BindCountData>> alleleBindCounts,
                final BufferedWriter summaryWriter, final BufferedWriter peptideWriter)
        {
            mTaskId = taskId;
            mConfig = config;
            mPosWeightModel = posWeightModel;
            mAlleleTrainingData = alleleTrainingData;
            mAlleleBindCounts = alleleBindCounts;

            mAlleles = Lists.newArrayList();
            mRandomDistribution = new RandomPeptideDistribution(mConfig.RandomPeptides);
            mSummaryWriter = summaryWriter;
            mPeptideWriter = peptideWriter;
        }

        public final List<String> getAlleles() { return mAlleles; }

        @Override
        public Long call()
        {
            for(int i = 0; i < mAlleles.size(); ++i)
            {
                String allele = mAlleles.get(i);
                int trainingDataCount = mAlleleTrainingData.get(allele).values().stream().mapToInt(x -> x.size()).sum();

                NE_LOGGER.info("{}: scoring without allele({}) and {} training data items", mTaskId, allele, trainingDataCount);

                testWithExcludedAllele(allele);

                if(i > 0 && (i % 5) == 0)
                {
                    NE_LOGGER.info("{}: processed {} of {} allele", mTaskId, i, mAlleles.size());
                }
            }

            return (long)0;
        }

        private void testWithExcludedAllele(final String targetAllele)
        {
            final Map<String, Map<Integer, BindScoreMatrix>> alleleBindMatrices = Maps.newHashMap();

            buildPositionWeightMatrixData(targetAllele, alleleBindMatrices);

            FlankScores flankScores = new FlankScores();

            if(mConfig.ApplyFlanks)
            {
                FlankCounts flankCounts = new FlankCounts();

                for(Map.Entry<String, Map<Integer, List<BindData>>> alleleEntry : mAlleleTrainingData.entrySet())
                {
                    String allele = alleleEntry.getKey();

                    if(allele.equals(targetAllele))
                        continue;

                    for(List<BindData> pepLenEntry : alleleEntry.getValue().values())
                    {
                        for(BindData bindData : pepLenEntry)
                        {
                            flankCounts.processBindData(bindData);
                        }
                    }
                }

                flankScores.createMatrix(flankCounts.getBindCounts());
            }

            // build distributions
            mRandomDistribution.buildDistribution(alleleBindMatrices, flankScores);

            final Map<String, Map<Integer, List<BindData>>> allelePeptideData = Maps.newHashMap();
            allelePeptideData.put(targetAllele, mAlleleTrainingData.get(targetAllele));

            BindScorer scorer = new BindScorer(allelePeptideData, alleleBindMatrices, mRandomDistribution, flankScores, null);
            scorer.runScoring();

            BindingLikelihood bindingLikelihood = new BindingLikelihood();
            bindingLikelihood.buildAllelePeptideLikelihoods(allelePeptideData, null);

            mRandomDistribution.buildLikelihoodDistribution(alleleBindMatrices, flankScores, bindingLikelihood, null);

            runScoring(targetAllele, alleleBindMatrices, flankScores, bindingLikelihood);
        }

        private void buildPositionWeightMatrixData(final String targetAllele,
                final Map<String, Map<Integer, BindScoreMatrix>> alleleBindMatrices)
        {
            Map<Integer, List<BindCountData>> countsByLength = Maps.newHashMap();
            mConfig.RequiredPeptideLengths.forEach(x -> countsByLength.put(x, Lists.newArrayList()));

            Map<Integer, BindCountData> targetAllelePepLenMap = mAlleleBindCounts.get(targetAllele);

            // now factor each allele's weighted counts into all the others
            for(Integer peptideLength : mConfig.RequiredPeptideLengths)
            {
                List<BindCountData> allBindCounts = Lists.newArrayList();
                countsByLength.put(peptideLength, allBindCounts);

                for(Map<Integer, BindCountData> pepLenMap : mAlleleBindCounts.values())
                {
                    BindCountData bindCounts = pepLenMap.get(peptideLength);

                    // exclude the target allele from contributing its counts
                    if(bindCounts.Allele.equals(targetAllele))
                        continue;

                    allBindCounts.add(bindCounts);
                }

                BindCountData bindCounts = targetAllelePepLenMap.get(peptideLength);

                // clear any previously set values
                MatrixUtils.clear(bindCounts.getFinalWeightedCounts());

                mPosWeightModel.buildFinalWeightedCounts(bindCounts, allBindCounts);
            }

            Map<Integer, BindCountData> pepLenBindCounts = mAlleleBindCounts.get(targetAllele);
            Map<Integer, BindScoreMatrix> peptideLengthMatrixMap = Maps.newHashMap();
            alleleBindMatrices.put(targetAllele, peptideLengthMatrixMap);

            for(BindCountData bindCounts : pepLenBindCounts.values())
            {
                BindScoreMatrix matrix = mPosWeightModel.createMatrix(bindCounts);
                peptideLengthMatrixMap.put(matrix.PeptideLength, matrix);
            }
        }

        private void runScoring(
                final String targetAllele, final Map<String, Map<Integer, BindScoreMatrix>> alleleBindMatrices,
                final FlankScores flankScores, BindingLikelihood bindingLikelihood)
        {
            NE_LOGGER.debug("{}: scoring targeted allele({})", mTaskId, targetAllele);

            // score the target allele having left it's training data out
            final Map<Integer, List<BindData>> pepLenBindDataMap = mAlleleTrainingData.get(targetAllele);

            Map<Integer, BindScoreMatrix> pepLenMatrixMap = alleleBindMatrices.get(targetAllele);

            if(pepLenMatrixMap == null)
            {
                NE_LOGGER.warn("allele({}) has no validation data", targetAllele);
                return;
            }

            List<AucData> alleleAucData = Lists.newArrayList();
            TprCalc alleleTprCalc = new TprCalc();

            for(Map.Entry<Integer, List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                int peptideLength = pepLenEntry.getKey();
                final List<BindData> bindDataList = pepLenEntry.getValue();

                TprCalc pepLenTprCalc = new TprCalc();

                for(BindData bindData : bindDataList)
                {
                    BindScoreMatrix matrix = pepLenMatrixMap.get(peptideLength);

                    BindScorer.calcScoreData(
                            bindData, matrix, flankScores, mRandomDistribution, bindingLikelihood, null, null);

                    writePeptideResults(mPeptideWriter, targetAllele, bindData);

                    alleleTprCalc.addRank(bindData.likelihoodRank());
                    pepLenTprCalc.addRank(bindData.likelihoodRank());

                    alleleAucData.add(new AucData(true, bindData.likelihoodRank(), true));
                    // writeResults(targetAllele, bindData);
                }

                writeSummaryResults(mSummaryWriter, targetAllele, String.valueOf(peptideLength), pepLenTprCalc.calc(), 0);
            }

            double aucPerc = AucCalc.calcPercentilesAuc(alleleAucData, Level.TRACE);
            writeSummaryResults(mSummaryWriter, targetAllele, "ALL", alleleTprCalc.calc(), aucPerc);
        }
    }

    private BufferedWriter initialiseSummaryWriter()
    {
        try
        {
            String outputFile = BindCommon.formFilename(mConfig.OutputDir, "validation_summary", mConfig.OutputId);
            BufferedWriter writer = createBufferedWriter(outputFile, false);
            writer.write("ExcludedAllele,PeptideLength,TPR,AUC");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise validation results file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeSummaryResults(
            final BufferedWriter writer, final String excludedAllele, final String peptideLength, double tpr, double auc)
    {
        try
        {
            writer.write(String.format("%s,%s,%.4f,%.4f",
                    excludedAllele, peptideLength, tpr, auc));

            writer.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write validation results file: {}", e.toString());
        }
    }

    private BufferedWriter initialisePeptideWriter()
    {
        try
        {
            String outputFile = BindCommon.formFilename(mConfig.OutputDir, "validation_scores", mConfig.OutputId);
            BufferedWriter writer = createBufferedWriter(outputFile, false);
            writer.write("Allele,Peptide,Score,Rank,Likelihood,LikelihoodRank");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise peptide scores file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writePeptideResults(final BufferedWriter writer, final String excludedAllele, final BindData bindData)
    {
        try
        {
            writer.write(String.format("%s,%s,%.4f,%.6f,%.6f,%.6f",
                    excludedAllele, bindData.Peptide, bindData.score(),
                    bindData.rankPercentile(), bindData.likelihood(), bindData.likelihoodRank()));

            writer.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide scores file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        TrainConfig.addCmdLineArgs(options);
        options.addOption(TEST_ALLELES_FILE, true, "List of alleles to leave out, otherwise will do all in training set");

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        ValidationRoutines validation = new ValidationRoutines(cmd);
        validation.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
