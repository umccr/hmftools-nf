package com.hartwig.hmftools.cobalt.norm;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.addGcProfilePath;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class NormalisationConfig
{
    public final List<String> SampleIds;
    public final String AmberDir;
    public final String CobaltDir;
    public final String TargetRegionsBed;
    public final String GcProfile;
    public final String OutputFile;
    public final String DetailedFile;
    public final RefGenomeVersion RefGenVersion;

    private static final String COBALT_DIR = "cobalt_dir";
    private static final String AMBER_DIR = "amber_dir";
    private static final String TARGET_REGIONS_BED = "target_regions_bed";
    private static final String OUTPUT_FILE = "output_file";
    private static final String DETAILED_OUTPUT = "detailed_file";

    public NormalisationConfig(final CommandLine cmd)
    {
        SampleIds = loadSampleIdsFile(cmd);
        CobaltDir = cmd.getOptionValue(COBALT_DIR);
        AmberDir = cmd.getOptionValue(AMBER_DIR);
        GcProfile = cmd.getOptionValue(GC_PROFILE);
        TargetRegionsBed = cmd.getOptionValue(TARGET_REGIONS_BED);
        RefGenVersion = RefGenomeVersion.from(cmd);
        OutputFile = cmd.getOptionValue(OUTPUT_FILE);
        DetailedFile = cmd.getOptionValue(DETAILED_OUTPUT);
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, "Sample IDs file");
        options.addOption(AMBER_DIR, true, "Path to amber files");
        options.addOption(COBALT_DIR, true, "Path to cobalt files");
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        options.addOption(TARGET_REGIONS_BED, true, "Target regions BED file");
        options.addOption(OUTPUT_FILE, true, "Output normalisation file");
        options.addOption(DETAILED_OUTPUT, true, "Output normalisation file");
        addGcProfilePath(options);
        addSampleIdFile(options);
        addLoggingOptions(options);
    }
}
