package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.CommandLine;
import org.jetbrains.annotations.NotNull;

public final class DriverGenePanelConfig
{
    public static final String DRIVER_GENE_PANEL_OPTION = "driver_gene_panel";
    public static final String DRIVER_GENE_PANEL_OPTION_DESC = "Path to driver gene panel";

    private DriverGenePanelConfig() {}

    public static boolean isConfigured(final ConfigBuilder configBuilder)
    {
        return configBuilder.hasValue(DRIVER_GENE_PANEL_OPTION);
    }

    public static List<DriverGene> loadDriverGenes(final ConfigBuilder configBuilder)
    {
        if(!isConfigured(configBuilder))
            return Lists.newArrayList();

        try
        {
            return DriverGeneFile.read(configBuilder.getValue(DRIVER_GENE_PANEL_OPTION));
        }
        catch(IOException e)
        {
            return Lists.newArrayList();
        }
    }

    public static boolean isConfigured(@NotNull final CommandLine cmd)
    {
        return cmd.hasOption(DRIVER_GENE_PANEL_OPTION);
    }

    public static List<DriverGene> driverGenes(final ConfigBuilder configBuilder) throws IOException
    {
        return DriverGeneFile.read(configBuilder.getValue(DRIVER_GENE_PANEL_OPTION));
    }

    public static List<DriverGene> driverGenes(final CommandLine cmd) throws IOException
    {
        return DriverGeneFile.read(cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));
    }

    public static List<DriverGene> loadDriverGenes(final CommandLine cmd)
    {
        if(!isConfigured(cmd))
            return Lists.newArrayList();

        try
        {
            return DriverGeneFile.read(cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));
        }
        catch(IOException e)
        {
            return Lists.newArrayList();
        }
    }
}
