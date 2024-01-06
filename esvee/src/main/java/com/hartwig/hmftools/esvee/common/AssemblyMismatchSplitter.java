package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.common.AssemblyUtils.buildFromJunctionReads;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.read.Read;

public class AssemblyMismatchSplitter
{
    private final AssemblySequence mSequence;

    public AssemblyMismatchSplitter(final AssemblySequence sequence)
    {
        mSequence = sequence;
    }

    public List<AssemblySequence> splitOnMismatches(int minSequenceLength)
    {
        // every remaining mismatch should have 2+ (or whatever configured) supporting reads
        // build unique collections of mismatches for each long enough read
        Set<Read> validMismatchReads = mSequence.support().stream()
                .filter(x -> x.mismatches() > 0)
                .filter(x -> x.readRangeLength() > minSequenceLength) // since includes the junction base
                .map(x -> x.read()).collect(Collectors.toSet());

        if(validMismatchReads.isEmpty())
            return Collections.emptyList();

        // make a support for each read with mismatches and then add if unique
        Map<Read,SequenceMismatches> readSequenceMismatches = Maps.newHashMap();

        for(Map.Entry<Integer,BaseMismatches> entry : mSequence.mismatches().indexedBaseMismatches().entrySet())
        {
            BaseMismatches baseMismatches = entry.getValue();
            int assemblyIndex = entry.getKey();

            for(int j = 0; j < baseMismatches.Mismatches.length; ++j)
            {
                if(baseMismatches.Mismatches[j] == null)
                    continue;

                BaseMismatch baseMismatch = baseMismatches.Mismatches[j];

                Mismatch mismatch = new Mismatch(assemblyIndex, baseMismatch.base(), baseMismatch.MaxQual, baseMismatch.QualTotal);

                for(Read read : baseMismatch.Reads)
                {
                    if(!validMismatchReads.contains(read))
                        continue;

                    SequenceMismatches readMismatches = readSequenceMismatches.get(read);

                    if(readMismatches == null)
                    {
                        readMismatches = new SequenceMismatches(read, mismatch);
                        readSequenceMismatches.put(read, readMismatches);
                    }
                    else
                    {
                        readMismatches.Mismatches.add(mismatch);
                    }
                }
            }
        }

        // now check for identical or compatible (ie containing all) mismatches across the reads
        List<SequenceMismatches> uniqueSequenceMismatches = Lists.newArrayList();
        Set<Read> matchedReads = Sets.newHashSet();

        List<SequenceMismatches> sortedSequenceMismatches = readSequenceMismatches.values().stream()
                .sorted(Comparator.comparingInt(x -> x.Mismatches.size())).collect(Collectors.toList());

        for(SequenceMismatches readMismatches : sortedSequenceMismatches)
        {
            Read read = readMismatches.Reads.get(0);

            if(matchedReads.contains(read))
                continue;

            matchedReads.add(read);

            uniqueSequenceMismatches.add(readMismatches);

            for(SequenceMismatches otherReadMismatches : sortedSequenceMismatches)
            {
                Read otherRead = otherReadMismatches.Reads.get(0);

                if(matchedReads.contains(otherRead))
                    continue;

                if(readMismatches.matchesOrContains(otherReadMismatches))
                {
                    readMismatches.Reads.add(otherRead);
                    matchedReads.add(otherRead);
                }
            }
        }

        // add the 'initial' sequence from reads without mismatches
        List<Read> noMismatchReads = mSequence.support().stream().filter(x -> x.mismatches() == 0).map(x -> x.read()).collect(Collectors.toList());

        AssemblySequence initialSequence = buildFromJunctionReads(mSequence.initialJunction(), noMismatchReads, false);

        List<AssemblySequence> finalSequences = Lists.newArrayListWithCapacity(1 + uniqueSequenceMismatches.size());
        finalSequences.add(initialSequence);

        // and then add each unique mismatched sequence

        for(SequenceMismatches sequenceMismatches : uniqueSequenceMismatches)
        {
            List<Read> candidateReads = sequenceMismatches.Reads;

            AssemblySequence mismatchSequence = buildFromJunctionReads(mSequence.initialJunction(), candidateReads, false);

            for(AssemblySupport support : mSequence.support())
            {
                if(support.mismatches() > 0 && !candidateReads.contains(support.read()))
                {
                    if(mismatchSequence.checkReadMatches(support.read(), 0))
                    {
                        mismatchSequence.addRead(support.read(), false);
                    }
                }
            }

            finalSequences.add(mismatchSequence);
        }

        return finalSequences;
    }

    private class Mismatch
    {
        public final int Index;
        public final byte Base;
        public final byte MaxQual;
        public final int TotalQual;

        public Mismatch(final int index, final byte base, final byte maxQual, final int totalQual)
        {
            Index = index;
            Base = base;
            MaxQual = maxQual;
            TotalQual = totalQual;
        }

        public boolean matches(final Mismatch other) { return Index == other.Index && Base == other.Base; }

        public String toString() { return format("%d: %c", Index, (char)Base); }
    }

    private class SequenceMismatches
    {
        public final List<Read> Reads;
        public final List<Mismatch> Mismatches;

        public SequenceMismatches(final Read read, final Mismatch mismatch)
        {
            Reads = Lists.newArrayList(read);
            Mismatches = Lists.newArrayList(mismatch);
        }

        public boolean matchesOrContains(final SequenceMismatches other)
        {
            if(other.Mismatches.size() > Mismatches.size())
                return false;

            for(int i = 0; i < Mismatches.size(); ++i)
            {
                if(i >= other.Mismatches.size())
                    return true;

                if(!Mismatches.get(i).matches(other.Mismatches.get(i)))
                    return false;
            }

            return true;
        }
    }

}
