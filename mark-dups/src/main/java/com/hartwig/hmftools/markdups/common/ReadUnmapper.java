package com.hartwig.hmftools.markdups.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getChromosomeFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionEndFieldIndex;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getPositionStartFieldIndex;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_INDEX;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.Constants.UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX;
import static com.hartwig.hmftools.markdups.common.Constants.UNMAP_MAX_NON_OVERLAPPING_BASES;
import static com.hartwig.hmftools.markdups.common.Constants.UNMAP_MIN_HIGH_DEPTH;
import static com.hartwig.hmftools.markdups.common.Constants.UNMAP_MIN_SOFT_CLIP;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

import htsjdk.samtools.SAMRecord;

public class ReadUnmapper
{
    private final Map<String, List<HighDepthRegion>> mChrLocationsMap; // keyed by chromosome start
    private boolean mEnabled;
    private final UnmapStats mStats;

    public ReadUnmapper(final String filename)
    {
        this(loadUnmapRegions(filename));

        if(mEnabled)
        {
            MD_LOGGER.info("loaded {} unmapping regions from {}",
                    mChrLocationsMap.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
    }

    public ReadUnmapper(final Map<String, List<HighDepthRegion>> chrLocationsMap)
    {
        mChrLocationsMap = chrLocationsMap;
        mEnabled = mChrLocationsMap != null && !mChrLocationsMap.isEmpty();
        mStats = new UnmapStats();
    }

    public List<HighDepthRegion> getRegions(final String chromosome)
    {
        return mChrLocationsMap.get(chromosome);
    }

    public boolean enabled() { return mEnabled; }

    public UnmapStats stats() { return mStats; }

    public boolean checkTransformRead(final SAMRecord read, final List<HighDepthRegion> readRegions)
    {
        if(!mEnabled || readRegions == null)
            return false;

        /* Criteria for unmapping a read:
            - falls within a region of high depth
            - discordant - INV, BND, one read unmapped or fragment length > 1000
            - soft-clip bases > 20

           Scenarios & logic:
           - both read and mate are already unmapped, then nothing to do
           - check the read's coords, its mate's coords and any supplementary alignment coords vs the loaded unmapping regions
           - for any overlap with a non-high-depth region, additionally check discordant and soft-lip bases as above
           -
           - supplementaries - unmap if their primary will be or if they need to be
         */

        // first check the read's alignment itself
        boolean readUnmapped = read.getReadUnmappedFlag();
        boolean mateUnmapped = mateUnmapped(read);

        if(readUnmapped && mateUnmapped)
            return false; // nothing to do and there won't be a supplementary

        boolean isSupplementary = read.getSupplementaryAlignmentFlag();

        boolean unmapRead = false;

        if(!readUnmapped)
        {
            UnmapReason unmapReason = checkUnmapRead(read, readRegions);

            if(unmapReason != UnmapReason.NONE)
            {
                unmapRead = true;

                switch(unmapReason)
                {
                    case HIGH_DEPTH:
                        mStats.HighDepthCount.incrementAndGet();
                        break;

                    case SOFT_CLIP:
                        mStats.LongSoftClipCount.incrementAndGet();
                        break;

                    case CHIMERIC:
                        mStats.ChimericCount.incrementAndGet();
                        break;

                    default:
                        break;
                }
            }
            else if(isSupplementary)
            {
                // supplementaries are unmapped if their primary will be unmapped
                unmapRead = isSupplementary && checkUnmapSupplementaryRead(read);
            }
        }

        if(unmapRead && isSupplementary)
        {
            // these will be dropped from the BAM
            setUnmappedAttributes(read);
            mStats.SupplementaryCount.incrementAndGet();
            return true;
        }

        boolean unmapMate = !mateUnmapped && checkUnmapMate(read);

        if(unmapRead)
        {
            unmapReadAlignment(read, mateUnmapped, unmapMate);
        }

        if(unmapMate)
        {
            // Note that for supplementary reads, at this point readUnmapped will be false, since a supplementary must be mapped. Also,
            // unmapRead is false, else we would have returned already. Therefore, this will only modify the mate properties and attributes.
            unmapMateAlignment(read, readUnmapped, unmapRead);
        }

        if((readUnmapped || unmapRead) && (mateUnmapped || unmapMate))
        {
            mStats.UnmappedCount.incrementAndGet();
        }
        else if(unmapRead)
        {
            mStats.ReadCount.incrementAndGet();
        }
        else if(unmapMate)
        {
            mStats.MateCount.incrementAndGet();
        }

        // nothing more to do for a supplementary whose primary or mate has been unmapped, or if it was unmapped
        if(isSupplementary)
        {
            return unmapMate;
        }

        boolean unmapSuppAlignment = false;

        if(!unmapRead && checkUnmapSupplementaryAlignments(read))
        {
            clearSupplementaryAlignment(read);
            mStats.SuppAlignmentCount.incrementAndGet();
            unmapSuppAlignment = true;
        }

        return unmapRead || unmapMate || unmapSuppAlignment;
    }

    private enum UnmapReason
    {
        HIGH_DEPTH,
        SOFT_CLIP,
        CHIMERIC,
        NONE;
    }

    private static UnmapReason checkUnmapRead(final SAMRecord read, final List<HighDepthRegion> readRegions)
    {
        RegionMatchType matchType = readMaxDepthRegionOverlap(read.getAlignmentStart(), read.getAlignmentEnd(), readRegions);

        if(matchType == RegionMatchType.NONE)
            return UnmapReason.NONE;

        if(matchType == RegionMatchType.HIGH_DEPTH)
            return UnmapReason.HIGH_DEPTH;

        if(getSoftClipCountFromCigarStr(read.getCigarString()) > UNMAP_MIN_SOFT_CLIP)
            return UnmapReason.SOFT_CLIP;

        if(isChimericRead(read, false))
            return UnmapReason.CHIMERIC;

        return UnmapReason.NONE;
    }

    private boolean checkUnmapMate(final SAMRecord read)
    {
        RegionMatchType matchType = mateMaxDepthRegionOverlap(read);

        if(matchType == RegionMatchType.NONE)
            return false;

        if(matchType == RegionMatchType.HIGH_DEPTH)
            return true;

        if(read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            final String mateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);

            if(getSoftClipCountFromCigarStr(mateCigar) > UNMAP_MIN_SOFT_CLIP)
                return true;
        }

        if(isChimericRead(read, true))
            return true;

        return false;
    }

    private boolean checkUnmapSupplementaryRead(final SAMRecord read)
    {
        final SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);

        if(suppData == null)
            return false;

        if(checkUnmapSupplementaryAlignment(suppData, read.getReadBases().length))
            return true;

        // We don't want to drop supplementaries based on chimeric read criteria, since this is a criteria for a primary read and its
        // primary mate. However, when checking if we are unmapping a supplementary based on whether its primary is unmapped, we do check
        // the chimeric read criteria.
        if(isSupplementaryChimericRead(read))
            return true;

        return false;
    }

    private boolean checkUnmapSupplementaryAlignments(final SAMRecord read)
    {
        // checks if any supplementary alignment qualifies for unmapping
        List<SupplementaryReadData> alignments = SupplementaryReadData.extractAlignments(read);

        if(alignments == null)
            return false;

        int readBaseLength = read.getReadBases().length;

        for(SupplementaryReadData suppData : alignments)
        {
            if(checkUnmapSupplementaryAlignment(suppData, readBaseLength))
                return true;
        }

        return false;
    }

    private boolean checkUnmapSupplementaryAlignment(final SupplementaryReadData suppData, int readBaseLength)
    {
        if(suppData == null)
            return false;

        RegionMatchType matchType = supplementaryMaxDepthRegionOverlap(suppData, readBaseLength);

        if(matchType == RegionMatchType.NONE)
            return false;

        if(matchType == RegionMatchType.HIGH_DEPTH)
            return true;

        if(getSoftClipCountFromCigarStr(suppData.Cigar) >= UNMAP_MIN_SOFT_CLIP)
            return true;

        return false;
    }

    private enum RegionMatchType
    {
        NONE,
        HIGH_DEPTH,
        OTHER;
    }

    private final static int READ_END_POS_BUFFER = 40; // to account for long deletes in a read that push back the alignment end

    @VisibleForTesting
    public RegionMatchType mateMaxDepthRegionOverlap(final SAMRecord read)
    {
        // Returns -1 if there is no overlap.
        final List<HighDepthRegion> mateRegions = mChrLocationsMap.get(read.getMateReferenceName());

        if(mateRegions == null)
            return RegionMatchType.NONE;

        // coordinate check must be identical to how the mate checks itself, which requires knowledge of aligned bases
        if(read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            int approxMateEnd = read.getMateAlignmentStart() + read.getReadBases().length + READ_END_POS_BUFFER;
            return readMaxDepthRegionOverlap(
                    read.getMateAlignmentStart(), approxMateEnd, read.getStringAttribute(MATE_CIGAR_ATTRIBUTE), mateRegions);
        }
        else
        {
            // TODO: this won't work for split reads
            int readBaseLength = read.getReadBases().length;
            return readMaxDepthRegionOverlap(
                    read.getMateAlignmentStart(), read.getMateAlignmentStart() + readBaseLength - 1, mateRegions);
        }
    }

    private RegionMatchType supplementaryMaxDepthRegionOverlap(final SupplementaryReadData suppReadData, final int readLength)
    {
        // Returns -1 if there is no overlap.
        final List<HighDepthRegion> suppRegions = mChrLocationsMap.get(suppReadData.Chromosome);

        if(suppRegions == null)
            return RegionMatchType.NONE;

        int approxSuppEnd = suppReadData.Position + readLength + READ_END_POS_BUFFER;
        return readMaxDepthRegionOverlap(suppReadData.Position, approxSuppEnd, suppReadData.Cigar, suppRegions);
    }

    private static RegionMatchType readMaxDepthRegionOverlap(final int readStart, final int readEnd, final List<HighDepthRegion> regions)
    {
        RegionMatchType matchType = RegionMatchType.NONE;

        for(HighDepthRegion region : regions)
        {
            if(region.start() > readEnd)
                break;

            if(!positionsOverlap(readStart, readEnd, region.start(), region.end()))
                continue;

            if(positionsWithin(readStart, readEnd, region.start(), region.end()))
            {
                return region.maxDepth() >= UNMAP_MIN_HIGH_DEPTH ? RegionMatchType.HIGH_DEPTH : RegionMatchType.OTHER;
            }

            int overlapBases = min(region.end(), readEnd) - max(region.start(), readStart) + 1;
            int readLength = readEnd - readStart + 1;
            int nonOverlappingBases = readLength - overlapBases;

            if(nonOverlappingBases < UNMAP_MAX_NON_OVERLAPPING_BASES)
            {
                if(region.maxDepth() >= UNMAP_MIN_HIGH_DEPTH)
                    return RegionMatchType.HIGH_DEPTH;
                else
                    matchType = RegionMatchType.OTHER;
            }
        }

        return matchType;
    }

    private static RegionMatchType readMaxDepthRegionOverlap(
            final int readStart, int approxReadEnd, final String cigar, final List<HighDepthRegion> regions)
    {
        RegionMatchType matchType = RegionMatchType.NONE;

        int readEnd = -1;
        for(HighDepthRegion region : regions)
        {
            if(region.start() > approxReadEnd)
                break;

            if(!positionsOverlap(readStart, approxReadEnd, region.start(), region.end()))
                continue;

            // now test with precise read end alignment from the Cigar
            if(readEnd < 0)
            {
                readEnd = getReadEndFromCigarStr(readStart, cigar);
            }

            int overlapBases = min(region.end(), readEnd) - max(region.start(), readStart) + 1;
            int readLength = readEnd - readStart + 1;
            int nonOverlappingBases = readLength - overlapBases;

            if(nonOverlappingBases < UNMAP_MAX_NON_OVERLAPPING_BASES)
            {
                if(region.maxDepth() >= UNMAP_MIN_HIGH_DEPTH)
                    return RegionMatchType.HIGH_DEPTH;
                else
                    matchType = RegionMatchType.OTHER;
            }
        }

        return matchType;
    }

    private static int getReadEndFromCigarStr(final int readStart, final String cigarStr)
    {
        int currentPosition = readStart;
        int elementLength = 0;

        for(int i = 0; i < cigarStr.length(); ++i)
        {
            final char c = cigarStr.charAt(i);
            final boolean isAddItem = (c == 'D' || c == 'M' || c == 'N');

            if(isAddItem)
            {
                currentPosition += elementLength;
                elementLength = 0;
                continue;
            }

            final int digit = c - '0';
            if(digit >= 0 && digit <= 9)
            {
                elementLength = elementLength * 10 + digit;
            }
            else
            {
                elementLength = 0;
            }
        }

        // always pointing to the start of the next element, so need to move back a base
        return currentPosition - 1;
    }

    @VisibleForTesting
    public static int getSoftClipCountFromCigarStr(final String cigarStr)
    {
        int softClipCount = 0;
        int elementLength = 0;
        for(int i = 0; i < cigarStr.length(); ++i)
        {
            final char c = cigarStr.charAt(i);
            if(c == 'S')
            {
                softClipCount += elementLength;
                elementLength = 0;
                continue;
            }

            int digit = c - '0';
            if(digit >= 0 && digit <= 9)
            {
                elementLength = elementLength * 10 + digit;
            }
            else
            {
                elementLength = 0;
            }
        }

        return softClipCount;
    }

    private static boolean isChimericRead(final SAMRecord record, boolean checkForMate)
    {
        boolean readUnmapped = record.getReadUnmappedFlag();
        boolean mateUnmapped = mateUnmapped(record);

        if(!checkForMate && readUnmapped)
            return false;

        if(checkForMate && mateUnmapped)
            return false;

        // or a fragment length outside the expected maximum
        if(abs(record.getInferredInsertSize()) > UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX)
            return true;

        // an unmapped mate
        if(!checkForMate && mateUnmapped)
            return true;

        if(checkForMate && readUnmapped)
            return true;

        if(record.getReadPairedFlag())
        {
            // inter-chromosomal
            if(!record.getReferenceName().equals(record.getMateReferenceName()))
                return true;

            // inversion
            if(record.getReadNegativeStrandFlag() == mateNegativeStrand(record))
                return true;
        }

        return false;
    }

    private static boolean isSupplementaryChimericRead(final SAMRecord record)
    {
        // check whether the primary read (detailed in the supp attribute) is chimeric with the mate - note the primary is always mapped
        SupplementaryReadData suppReadData = SupplementaryReadData.extractAlignment(record);

        // insert size is not populated for supplementaries

        if(mateUnmapped(record))
            return true;

        if(record.getReadPairedFlag())
        {
            // inter-chromosomal
            if(!suppReadData.Chromosome.equals(record.getMateReferenceName()))
                return true;

            // inversion
            boolean primaryReadOnNegativeStrand = suppReadData.Strand == SupplementaryReadData.SUPP_NEG_STRAND;

            if(primaryReadOnNegativeStrand == mateNegativeStrand(record))
                return true;
        }

        return false;
    }

    private static void setUnmappedAttributes(final SAMRecord read)
    {
        read.setReadUnmappedFlag(true);
        read.setMappingQuality(0);

        // store the original mapping
        setUnmapCoordsAttribute(read, read.getReferenceName(), read.getAlignmentStart());
    }

    public static void unmapReadAlignment(final SAMRecord read, final boolean mateUnmapped, final boolean unmapRead)
    {
        setUnmappedAttributes(read);

        read.setProperPairFlag(false);

        // clear insert size
        read.setInferredInsertSize(0);

        // clear reference index, reference name and alignment start
        if(mateUnmapped || unmapRead)
        {
            // set both to unknown
            read.setAlignmentStart(0);
            read.setReferenceIndex(NO_CHROMOSOME_INDEX);
            read.setReferenceName(NO_CHROMOSOME_NAME);
        }
        else
        {
            // set to primary's details
            read.setAlignmentStart(read.getMateAlignmentStart());
            read.setReferenceIndex(read.getMateReferenceIndex());
            read.setReferenceName(read.getMateReferenceName());
        }

        // If mate is unmapped, have to clear its reference name and alignment start here, since no unmapping logic will be run for the
        // mate.
        if(mateUnmapped)
        {
            read.setMateAlignmentStart(0);
            read.setMateReferenceIndex(NO_CHROMOSOME_INDEX);
            read.setMateReferenceName(NO_CHROMOSOME_NAME);
        }

        // clear cigar if present
        read.setCigarString(NO_CIGAR);

        // clear the supplementary data too
        clearSupplementaryAlignment(read);
    }

    private static final String UNMAPP_COORDS_DELIM = ":";

    public static String[] parseUnmappedCoords(final SAMRecord read)
    {
        return read.getStringAttribute(UNMAP_ATTRIBUTE).split(UNMAPP_COORDS_DELIM, 2);
    }

    private static void setUnmapCoordsAttribute(final SAMRecord read, final String chromosome, final int position)
    {
        read.setAttribute(UNMAP_ATTRIBUTE, chromosome + UNMAPP_COORDS_DELIM + position);
    }

    public static void unmapMateAlignment(final SAMRecord read, final boolean readUnmapped, final boolean unmapRead)
    {
        // set flag mate unmapped
        read.setMateUnmappedFlag(true);
        read.setProperPairFlag(false);

        // store the mate's original mapping
        setUnmapCoordsAttribute(read, read.getMateReferenceName(), read.getMateAlignmentStart());

        // clear insert size
        read.setInferredInsertSize(0);

        // set mate reference index, reference name and alignment start
        if(readUnmapped || unmapRead)
        {
            // set both to unknown
            read.setMateAlignmentStart(0);
            read.setMateReferenceIndex(NO_CHROMOSOME_INDEX);
            read.setMateReferenceName(NO_CHROMOSOME_NAME);

            // clear the read's coords too - this is not handled elsewhere
            if(readUnmapped)
            {
                read.setAlignmentStart(0);
                read.setReferenceIndex(NO_CHROMOSOME_INDEX);
                read.setReferenceName(NO_CHROMOSOME_NAME);
            }
        }
        else
        {
            // set to primary's details
            read.setMateAlignmentStart(read.getAlignmentStart());
            read.setMateReferenceIndex(read.getReferenceIndex());
            read.setMateReferenceName(read.getReferenceName());
        }

        // clear mate cigar if present
        if(read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            read.setAttribute(MATE_CIGAR_ATTRIBUTE, null);
        }
    }

    public static void clearSupplementaryAlignment(final SAMRecord read)
    {
        // remove supplementary alignment details from attributes
        if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, null);
        }
    }

    private static Map<String, List<HighDepthRegion>> loadUnmapRegions(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            final String delim = FileDelimiters.inferFileDelimiter(filename);

            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(header, delim);
            int chrIndex = getChromosomeFieldIndex(fieldIndexMap);
            int posStartIndex = getPositionStartFieldIndex(fieldIndexMap);
            int posEndIndex = getPositionEndFieldIndex(fieldIndexMap);
            int depthIndex = fieldIndexMap.get("MaxDepth");

            Map<String, List<HighDepthRegion>> chrLocationsMap = Maps.newHashMap();
            HighDepthRegion lastRegion = null;

            for(String line : lines)
            {
                String[] values = line.split(delim, -1);

                String chromosome = values[chrIndex];

                List<HighDepthRegion> regions = chrLocationsMap.get(chromosome);
                if(regions == null)
                {
                    regions = Lists.newArrayList();
                    chrLocationsMap.put(chromosome, regions);
                    lastRegion = null;
                }

                int posStart = Integer.parseInt(values[posStartIndex]);
                int posEnd = Integer.parseInt(values[posEndIndex]);
                int maxDepth = Integer.parseInt(values[depthIndex]);

                HighDepthRegion region = new HighDepthRegion(posStart, posEnd, maxDepth);

                if(lastRegion != null)
                {
                    if(lastRegion.overlaps(region))
                    {
                        MD_LOGGER.warn("unmap regions overlap: current({}) next({})", lastRegion, region);
                    }
                    else if(region.end() < lastRegion.start())
                    {
                        MD_LOGGER.warn("unmap regions not sorted: current({}) next({})", lastRegion, region);
                    }

                    // to difficult to merge
                }

                lastRegion = region;

                regions.add(region);
            }

            return chrLocationsMap;
        }
        catch(IOException e)
        {
            MD_LOGGER.error("failed to read high-depth regions file {}: {}", filename, e.toString());
            System.exit(1);
            return null;
        }
    }

    @VisibleForTesting
    public void addRegion(final String chromosome, final HighDepthRegion region)
    {
        List<HighDepthRegion> regions = mChrLocationsMap.get(chromosome);
        if(regions == null)
        {
            regions = Lists.newArrayList();
            mChrLocationsMap.put(chromosome, regions);
        }

        regions.add(region);

        mEnabled = true;
    }
}