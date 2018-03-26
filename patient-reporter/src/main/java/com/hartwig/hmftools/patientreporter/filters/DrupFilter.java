package com.hartwig.hmftools.patientreporter.filters;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.io.reader.LineReader;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;

public class DrupFilter implements Predicate<VariantReport> {

    @NotNull
    private final Set<String> includedGenes;

    public DrupFilter(@NotNull final String drupGenesCsv) throws IOException {
        final List<String> geneLines = LineReader.build().readLines(new File(drupGenesCsv).toPath(),
                line -> line.length() > 0);
        this.includedGenes = Sets.newHashSet(geneLines);
    }

    @Override
    public boolean test(final VariantReport variantReport) {
        return includedGenes.contains(variantReport.gene());
    }
}
