package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FeatureTypeExtractor {

    public static final Set<String> AMPLIFICATIONS = Sets.newHashSet("Amplification",
            "amplification",
            "AMPLIFICATION",
            "amp",
            "overexpression",
            "over exp",
            "amp over exp",
            "OVEREXPRESSION",
            "Overexpression");

    public static final Set<String> DELETIONS = Sets.newHashSet("Deletion",
            "deletion",
            "DELETION",
            "del",
            "undexpression",
            "dec exp",
            "UNDEREXPRESSION",
            "loss",
            "LOSS",
            "Copy Number Loss");

    public static final Set<String> SEARCH_FUSION_PAIRS = Sets.newHashSet("Fusion",
            "Disruptive Inframe Deletion",
            "Gene Fusion",
            "fusion",
            "EGFR-KDD",
            "Transcript Regulatory Region Fusion",
            "FGFR3 - BAIAP2L1 Fusion",
            "nonsense");
    public static final Set<String> SEARCH_FUSION_PROMISCUOUS =
            Sets.newHashSet("REARRANGEMENT", "Fusions", "fusion", "rearrange", "Transcript Fusion", "FUSION", "FUSIONS");

    public static final Set<String> IGNORE = Sets.newHashSet("3' EXON DELETION");

    public static final Set<String> INTERNAL_FUSION =
            Sets.newHashSet("is_deletion", "EGFRvIII", "EGFRvV", "EGFRvII", "ITD");

    public static final Set<String> GENE_LEVEL = Sets.newHashSet("Gain-of-function Mutations",
            "Gain-of-Function",
            "Oncogenic Mutations",
            "MUTATION",
            "act mut",
            "pos",
            "positive",
            "inact mut",
            "biallelic inactivation",
            "negative",
            "Loss Of Function Variant",
            "Loss Of Heterozygosity",
            "TRUNCATING MUTATION",
            "Truncating Mutations",
            "mutant",
            "mut",
            "gene_only",
            "ACTIVATING MUTATION",
            "DELETERIOUS MUTATION",
            "feature_truncation",
            "FRAMESHIFT TRUNCATION",
            "FRAMESHIFT MUTATION",
            "SPLICE VARIANT 7",
            "Splice",
            "DNMT3B7",
            "LCS6-variant",
            "AR-V7",
            "ARv567es");

    public static final Set<String> GENE_EXON = Sets.newHashSet("exon");

    public static final Set<String> SIGNATURES = Sets.newHashSet("Microsatellite Instability-High");

    private FeatureTypeExtractor() {
    }

    @NotNull
    public static FeatureType extractType(@NotNull Feature feature) {
        return extractType(feature.name(),
                feature.biomarkerType(),
                feature.provenanceRule(),
                ProteinAnnotationExtractor.extractProteinAnnotation(feature));
    }

    @NotNull
    public static FeatureType extractType(@NotNull String featureName, @Nullable String biomarkerType, @Nullable String provenanceRule,
            @NotNull String proteinAnnotation) {
        String event = Strings.EMPTY;
        if (featureName.toLowerCase().contains("exon")) {
            event = "exon";
        }

        if (biomarkerType != null && provenanceRule != null) {
            if (featureName.contains("+") && !featureName.toLowerCase().contains("c.") && !featureName.contains(">")) {
                return FeatureType.COMBINED;
            } else if (featureName.contains("insertion")) {
                int countInsertion = featureName.split("insertion").length - 1;
                if (countInsertion > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("deletion")) {
                int countDeletion = featureName.split("deletion").length - 1;
                if (countDeletion > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("frameshift")) {
                int countFrameshift = featureName.split("frameshift").length - 1;
                if (countFrameshift > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("insertions") && featureName.contains("deletion")) {
                int countCombined = (featureName.split("insertion").length - 1) + (featureName.split("deletion").length - 1);
                if (countCombined > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("splice")) {
                int countSplice = featureName.split("splice").length - 1;
                if (countSplice > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (proteinAnnotation.equals("p61BRAF-V600E")) {
                return FeatureType.COMBINED;
            }
        }
        if (DetermineHotspot.isHotspot(proteinAnnotation)) {
            return FeatureType.HOTSPOT;
        } else if (FeatureTypeExtractor.SIGNATURES.contains(featureName)) {
            return FeatureType.SIGNATURE;
        } else if (DetermineCopyNumber.isAmplification(featureName, biomarkerType)) {
            return FeatureType.AMPLIFICATION;
        } else if (DetermineCopyNumber.isDeletion(featureName, biomarkerType) && !featureName.toLowerCase().contains("exon")) {
            return FeatureType.DELETION;
        } else if (DetermineFusion.isFusion(featureName, biomarkerType, provenanceRule, proteinAnnotation) && !featureName.contains(
                "p61BRAF")) {
            return FeatureType.FUSION_PAIR;
        } else if (DetermineFusion.isFusionPromiscuous(featureName, biomarkerType, provenanceRule, proteinAnnotation)) {
            return FeatureType.FUSION_PROMISCUOUS;
        } else if (FeatureTypeExtractor.GENE_EXON.contains(event)) {
            if (featureName.toLowerCase().contains("deletion") || featureName.toLowerCase().contains("insertion")
                    || featureName.toLowerCase().contains("proximal") || featureName.toLowerCase().contains("mutation")
                    || featureName.toLowerCase().contains("splice site insertion") || featureName.toLowerCase().contains("frameshift")) {
                return FeatureType.GENE_RANGE_EXON;
            }
        } else if (proteinAnnotation.length() > 1 && proteinAnnotation.substring(proteinAnnotation.length() - 1).equals("X")) {
            return FeatureType.GENE_RANGE_CODON;
        } else if (proteinAnnotation.length() >= 1 && isValidSingleCodonRange(proteinAnnotation)) {
            return FeatureType.GENE_RANGE_CODON;
        } else if (!DetermineHotspot.isHotspot(proteinAnnotation)) {
            if (FeatureTypeExtractor.GENE_LEVEL.contains(biomarkerType) || FeatureTypeExtractor.GENE_LEVEL.contains(featureName)
                    || FeatureTypeExtractor.GENE_LEVEL.contains(provenanceRule) || FeatureTypeExtractor.GENE_LEVEL.contains(
                    proteinAnnotation)) {
                return FeatureType.GENE_LEVEL;
            }
        } return FeatureType.UNKNOWN;
    }

    private static boolean isValidSingleCodonRange(@NotNull String feature) {
        // Features are expected to look something like V600 (1 char - N digits)
        if (feature.length() < 3) {
            return false;
        }

        if (!Character.isLetter(feature.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(feature.charAt(1))) {
            return false;
        }

        if (feature.contains("*")) {
            return false;
        }

        if (feature.contains("/")) {
            return false;
        }

        if (feature.contains("fs")) {
            return false;
        }

        return Character.isDigit(feature.substring(feature.length() - 1).charAt(0));
    }

}
