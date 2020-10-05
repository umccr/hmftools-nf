package com.hartwig.hmftools.serve.actionability.range;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class RangeEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private RangeEvidenceAnalyzerFactory() {
    }

    @NotNull
    public static RangeEvidenceAnalyzer loadFromActionableRangeTsv(@NotNull String actionableRangeTsv) throws IOException {
        List<ActionableRange> actionableRanges = Lists.newArrayList();
        List<String> lines = Files.readAllLines(new File(actionableRangeTsv).toPath());

        // Skip header line for range
        for (String line : lines.subList(1, lines.size())) {
            actionableRanges.add(fromLine(line));
        }
        return new RangeEvidenceAnalyzer(actionableRanges);
    }

    @NotNull
    private static ActionableRange fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableActionableRange.builder()
                .gene(values[0])
                .chromosome(values[1])
                .start(values[2])
                .end(values[3])
                .mutationType(values[4])
                .source(values[5])
                .drug(values[6])
                .cancerType(values[7])
                .doid(values[8])
                .level(values[9])
                .direction(values[10])
                .build();
    }
}
