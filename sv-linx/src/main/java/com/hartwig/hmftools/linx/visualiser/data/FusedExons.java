package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.linx.visualiser.data.Exons.downstreamExons;
import static com.hartwig.hmftools.linx.visualiser.data.Exons.upstreamExons;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class FusedExons
{
    private static final String DELIMITER = "\t";

    public static void write(@NotNull final String fileName, @NotNull final List<FusedExon> fusedExons) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(fusedExons));
    }

    @NotNull
    public static List<FusedExon> fusedExons(@NotNull final Fusion fusion, @NotNull final List<Exon> exons)
    {
        final List<FusedExon> result = Lists.newArrayList();

        final List<Exon> upStreamExons = upstreamExons(fusion, exons);
        final List<Exon> downStreamExons = downstreamExons(fusion, exons);
        if (upStreamExons.isEmpty() || downStreamExons.isEmpty()) {
            return result;
        }

        final Exon firstUpstreamExon = upStreamExons.get(0);

        final long upGeneOffset = offset(fusion.strandUp(), firstUpstreamExon);
        final long upGeneStart = start(fusion.strandUp(), upGeneOffset, firstUpstreamExon);
        final long upGeneEnd = convert(fusion.strandUp(), upGeneOffset, fusion.positionUp());

        final ImmutableFusedExon.Builder upFusedExonBuilder = ImmutableFusedExon.builder()
                .sampleId(fusion.sampleId())
                .clusterId(fusion.clusterId())
                .fusion(fusion.name())
                .chromosome(fusion.chromosomeUp())
                .unadjustedGeneStart(fusion.strandUp() < 0 ? firstUpstreamExon.end() : firstUpstreamExon.start())
                .gene(fusion.geneUp())
                .geneStart(upGeneStart)
                .geneEnd(upGeneEnd);

        for (final Exon exon : upStreamExons)
        {
            final long exonStart = start(fusion.strandUp(), upGeneOffset, exon);
            final long exonEnd = end(fusion.strandUp(), upGeneOffset, exon);

            if (exonStart <= upGeneEnd)
            {
                final FusedExon fusedExon = upFusedExonBuilder
                        .start(exonStart)
                        .end(Math.min(exonEnd, upGeneEnd))
                        .rank(exon.rank())
                        .skipped(false) // TODO: Check exon skipped field
                        .build();
                result.add(fusedExon);
            }
        }

        final long downGeneOffset = fusion.positionDown();
        final long downGeneStart = convert(fusion.strandDown(), downGeneOffset, fusion.positionDown()) + upGeneEnd;
        final long downGeneEnd = end(fusion.strandDown(), downGeneOffset, downStreamExons.get(downStreamExons.size() - 1)) + upGeneEnd;

        final ImmutableFusedExon.Builder downFusedExonBuilder = ImmutableFusedExon.builder().from(upFusedExonBuilder.build())
                .chromosome(fusion.chromosomeDown())
                .unadjustedGeneStart(fusion.positionDown())
                .gene(fusion.geneDown())
                .geneStart(downGeneStart)
                .geneEnd(downGeneEnd);

        boolean intronicToExonicFusion = fusion.regionTypeUp().equals("Intronic") && fusion.regionTypeDown().equals("Exonic");

        for (int i = 0; i < downStreamExons.size(); i++)
        {
            final Exon exon = downStreamExons.get(i);

            final long exonStart = start(fusion.strandDown(), downGeneOffset, exon) + upGeneEnd;
            final long exonEnd = end(fusion.strandDown(), downGeneOffset, exon) + upGeneEnd;

            if (exonEnd > downGeneStart)
            {
                final FusedExon fusedExon = downFusedExonBuilder
                        .start(Math.max(exonStart, downGeneStart))
                        .end(exonEnd)
                        .rank(exon.rank())
                        .skipped(exon.rank() == 1 || (i == 0 && intronicToExonicFusion)) // TODO: Check exon skipped field
                        .build();
                result.add(fusedExon);
            }

        }

        return result;
    }


    @NotNull
    static List<String> toLines(@NotNull final List<FusedExon> exons)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        exons.stream().map(FusedExons::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("sampleId")
                .add("clusterId")
                .add("fusion")
                .add("gene")
                .add("geneStart")
                .add("geneEnd")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("rank")
                .add("truncated")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final FusedExon exon)
    {
        return new StringJoiner(DELIMITER)
                .add(exon.sampleId())
                .add(String.valueOf(exon.clusterId()))
                .add(String.valueOf(exon.fusion()))
                .add(String.valueOf(exon.gene()))
                .add(String.valueOf(exon.geneStart()))
                .add(String.valueOf(exon.geneEnd()))
                .add(String.valueOf(exon.chromosome()))
                .add(String.valueOf(exon.start()))
                .add(String.valueOf(exon.end()))
                .add(String.valueOf(exon.rank()))
                .add(String.valueOf(exon.skipped()))
                .toString();
    }


    private static long offset(int strand, @NotNull final Exon exon) {
        return strand < 0 ? exon.end() : exon.start();
    }

    private static long start(int strand, long offset, GenomeRegion region) {
        return strand < 0 ? offset - region.end(): region.start() - offset;
    }

    private static long end(int strand, long offset, GenomeRegion region) {
        return strand < 0 ? offset - region.start() : region.end() - offset;
    }

    private static long convert(int strand, long offset, long position) {
        return strand < 0 ? offset - position : position - offset;
    }
}
