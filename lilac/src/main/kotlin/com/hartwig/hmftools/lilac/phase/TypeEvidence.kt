package com.hartwig.hmftools.lilac.phase

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.read.Fragment

class TypeEvidence(private val minBaseCount: Int, private val minFragmentCount: Int) {

    fun evidence(aminoAcidFragments: List<Fragment>): List<PhasedEvidence> {
        val aminoAcidCounts = SequenceCount.aminoAcids(minBaseCount, aminoAcidFragments)

        val heterozygousIndices = aminoAcidCounts.heterozygousIndices()
        val heterozygousEvidence = ExtendedEvidence(minBaseCount, minFragmentCount, heterozygousIndices, aminoAcidFragments)

        val allEvidence = mutableSetOf<PhasedEvidence>()
        val initialEvidence = heterozygousEvidence.initialEvidence()
        var unprocessedEvidence = initialEvidence

        allEvidence.addAll(initialEvidence)

        while (unprocessedEvidence.isNotEmpty()) {
            val topEvidence = unprocessedEvidence[0]
            allEvidence.add(topEvidence)

            val newEvidence = heterozygousEvidence.extendConsecutive(topEvidence, allEvidence)
            allEvidence.addAll(newEvidence)

            val updatedEvidence = mutableSetOf<PhasedEvidence>()
            updatedEvidence.addAll(unprocessedEvidence.drop(1))
            updatedEvidence.addAll(newEvidence)

            unprocessedEvidence = updatedEvidence.sorted()
        }

        return longestEvidence(allEvidence)

    }

    private fun longestEvidence(evidence: Collection<PhasedEvidence>): List<PhasedEvidence> {
        fun Collection<PhasedEvidence>.otherContains(victim: PhasedEvidence): Boolean {
            return this.any { it != victim && it.contains(victim) }
        }
        return evidence
                .filter { !evidence.otherContains(it) }
                .sortedBy { it.aminoAcidIndices[0] }
    }


    private fun exclude(fragment: NucleotideFragment, gene: String, excludedNucleotides: Collection<Int>): Boolean {
        return fragment.alignedGene == gene && excludedNucleotides.any { fragment.containsNucleotide(it) }
    }

}


