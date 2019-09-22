package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.REF_GENOME;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.common.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidence;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidenceFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BamCountReader
{
    private IndexedFastaSequenceFile mIndexedFastaSequenceFile;
    private File mRefGenomeFile;
    private SamReader mTumorReader;

    private static final int DEFAULT_MIN_BASE_QUALITY = 13;
    private static final int DEFAULT_MIN_MAPPING_QUALITY = 1;

    private static final Logger LOGGER = LogManager.getLogger(BamCountReader.class);

    BamCountReader()
    {
        mIndexedFastaSequenceFile = null;
        mTumorReader = null;
        mRefGenomeFile = null;
    }

    public void initialise(final String refGenomeFile, IndexedFastaSequenceFile ifSeqFile)
    {
        mIndexedFastaSequenceFile = ifSeqFile;

        mRefGenomeFile = new File(refGenomeFile);
    }

    public void readBamCounts(final String bamFile, List<BachelorGermlineVariant> bachRecords)
    {
        LOGGER.debug("reading BAM file: {}", bamFile);

        mTumorReader = SamReaderFactory.makeDefault().referenceSequence(mRefGenomeFile).open(new File(bamFile));

        final Set<VariantHotspot> allHotspots = Sets.newHashSet();

        for(BachelorGermlineVariant variant : bachRecords)
        {
            VariantHotspot variantHotspot = ImmutableVariantHotspotImpl.builder()
                    .chromosome(variant.Chromosome)
                    .position(variant.Position)
                    .ref(variant.Ref)
                    .alt(variant.Alts)
                    .build();

            allHotspots.add(variantHotspot);
        }

        final VariantHotspotEvidenceFactory hotspotEvidenceFactory = new VariantHotspotEvidenceFactory(DEFAULT_MIN_MAPPING_QUALITY, DEFAULT_MIN_BASE_QUALITY, allHotspots);
        final List<VariantHotspotEvidence> tumorEvidence = hotspotEvidenceFactory.evidence(mIndexedFastaSequenceFile, mTumorReader);

        if(tumorEvidence.size() != bachRecords.size())
        {
            LOGGER.error("Incomplete BAM evidence read: evidenceCount({}) vs bachRecords({})", tumorEvidence.size(), bachRecords.size());
            return;
        }

        for(BachelorGermlineVariant variant : bachRecords)
        {
            for(VariantHotspotEvidence evidence : tumorEvidence)
            {
                if(evidence.chromosome().equals(variant.Chromosome) && evidence.position() == variant.Position
                && evidence.ref().equals(variant.Ref) && evidence.alt().equals(variant.Alts))
                {
                    variant.setTumorData(evidence.altSupport(), evidence.readDepth());

                    LOGGER.debug("chr({}) position({}) matched, counts(ref={} alt={} depth={})",
                            variant.Chromosome, variant.Position,
                            variant.getTumorRefCount(), variant.getTumorAltCount(), variant.getTumorReadDepth());

                    break;
                }
            }
        }
    }
}
