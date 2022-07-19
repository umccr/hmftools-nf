package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.hla.HlaFiles;
import com.hartwig.hmftools.common.hla.HlaType;
import com.hartwig.hmftools.common.hla.HlaTypeDetails;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadLilacData {

    private static final Logger LOGGER = LogManager.getLogger(LoadLilacData.class);

    private static final String SAMPLE = "sample";
    private static final String LILAC_OUTPUT_FILE = "lilac_output_file";
    private static final String LILAC_QC_METRICS_FILE = "lilac_qc_metrics_file";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        final String sample = cmd.getOptionValue(SAMPLE);
        final String lilacOutputFile = cmd.getOptionValue(LILAC_OUTPUT_FILE);
        final String lilacQcMetricsFile = cmd.getOptionValue(LILAC_QC_METRICS_FILE);

        LOGGER.info("Reading output from {} and using QC metrics from {}", lilacOutputFile, lilacQcMetricsFile);
        final HlaType type = HlaFiles.type(lilacOutputFile, lilacQcMetricsFile);
        final List<HlaTypeDetails> details = HlaFiles.typeDetails(lilacOutputFile);

        LOGGER.info("Persisting lilac data to db for {}", sample);
        dbAccess.writeHla(sample, type, details);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(LILAC_OUTPUT_FILE, true, "Lilac output file.");
        options.addOption(LILAC_QC_METRICS_FILE, true, "Lilac QC metrics file.");
        addDatabaseCmdLineArgs(options);
        return options;
    }

}