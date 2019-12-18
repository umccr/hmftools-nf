package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.nullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalJsonArray;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalJsonObject;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalNullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalStringList;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.string;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.stringList;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.toStringList;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAgreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightLeftRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchBreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchExonBoundaries;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchExonsInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchFusion;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchGRCh37Location;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchMutation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchParent;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchPosition;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTag;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTranscriptConsequencesGRCh37;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchWGSALocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchWGSAMap;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAgreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightLeftRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchBreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchExonBoundaries;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchExonsInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchFusion;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchGRCh37Location;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchMutation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchParent;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchPosition;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTag;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTranscriptConsequencesGRCh37;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchWGSALocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchWGSAMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class MolecularMatchObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(MolecularMatchObjectFactory.class);

    private static final List<Integer> EXPECTED_MOLECULARMATCH_AST_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LEFT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_LEFT_RIGHT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES = Lists.newArrayList(3, 29, 30, 31);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_CRITERIA_UNMET_SIZES = Lists.newArrayList(8, 9, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LOCATIONGRCH37_SIZES = Lists.newArrayList(9);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LOCATIONS_SIZES = Lists.newArrayList(3, 11);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_PREVALENCE_SIZES = Lists.newArrayList(4, 6);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_SOURCE_SIZES = Lists.newArrayList(8, 9, 10, 11, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TAGS_SIZES = Lists.newArrayList(3, 8, 9, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TIER_EXPLANATION_SIZES = Lists.newArrayList(4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES__GRCH37_SIZES = Lists.newArrayList(6);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_VARIANTINFO_SIZES = Lists.newArrayList(10);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_PARENTS_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_FUSIONS_SIZES = Lists.newArrayList(8);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_WGSAMAP_SIZES = Lists.newArrayList(7, 9);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_EXONSINFO_SIZES = Lists.newArrayList(3, 7);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_EXONSBOUNDARIES_SIZES =
            Lists.newArrayList(1, 2, 3, 5, 6, 7, 8, 9, 11, 13, 16, 17, 20, 21, 22, 24, 26, 27, 28, 29, 38, 41);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_POSITIONS_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_FUSIONDATA_SIZES = Lists.newArrayList(1, 2, 15, 17);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_BREG_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_AREG_SIZES = Lists.newArrayList(2);

    private MolecularMatchObjectFactory() {
    }

    @NotNull
    static MolecularMatch create(@NotNull JsonObject molecularMatchObject) {
        ViccDatamodelCheckerFactory.molecularMatchEntryChecker().check(molecularMatchObject);

        return ImmutableMolecularMatch.builder()
                .direction(string(molecularMatchObject, "direction"))
                .biomarkerClass(string(molecularMatchObject, "biomarkerClass"))
                .mutations(createMutations(molecularMatchObject.getAsJsonArray("mutations")))
                .variantInfos(createVariantInfos(molecularMatchObject.getAsJsonArray("variantInfo")))
                .prevalences(createPrevalences(molecularMatchObject.getAsJsonArray("prevalence")))
                .score(string(molecularMatchObject, "_score"))
                .sources(createSources(molecularMatchObject.getAsJsonArray("sources")))
                .clinicalSignificance(string(molecularMatchObject, "clinicalSignificance"))
                .tier(string(molecularMatchObject, "tier"))
                .tierExplanations(createTierExplanations(molecularMatchObject.getAsJsonArray("tierExplanation")))
                .ampcap(string(molecularMatchObject, "ampcap"))
                .civicValue(string(molecularMatchObject, "civic"))
                .regulatoryBody(string(molecularMatchObject, "regulatoryBody"))
                .regulatoryBodyApproved(string(molecularMatchObject, "regulatoryBodyApproved"))
                .guidelineBody(optionalString(molecularMatchObject, "guidelineBody"))
                .guidelineVersion(optionalString(molecularMatchObject, "guidelineVersion"))
                .includeCondition1(stringList(molecularMatchObject, "includeCondition1"))
                .includeMutation1(optionalStringList(molecularMatchObject, "includeMutation1"))
                .includeDrug1(optionalStringList(molecularMatchObject, "includeDrug1"))
                .includeStage0(optionalStringList(molecularMatchObject, "includeStage0"))
                .includeGene0(optionalStringList(molecularMatchObject, "includeDrug0"))
                .includeCondition0(stringList(molecularMatchObject, "includeCondition0"))
                .includeMutation0(optionalStringList(molecularMatchObject, "includeMutation0"))
                .criteriaMets(stringList(molecularMatchObject, "criteriaMet"))
                .criteriaUnmets(createCriteriaUnmets(molecularMatchObject.getAsJsonArray("criteriaUnmet")))
                .ast(createAst(molecularMatchObject.getAsJsonObject("ast")))
                .institutions(optionalStringList(molecularMatchObject, "institution"))
                .tags(createTags(molecularMatchObject.getAsJsonArray("tags")))
                .classifications(createClassifications(molecularMatchObject.getAsJsonArray("classifications")))
                .noTherapyAvailable(optionalString(molecularMatchObject, "noTherapyAvailable"))
                .therapeuticContexts(createTherapeuticContexts(molecularMatchObject.getAsJsonArray("therapeuticContext")))
                .sixtier(string(molecularMatchObject, "sixtier"))
                .mvld(string(molecularMatchObject, "mvld"))
                .autoGenerateNarrative(string(molecularMatchObject, "autoGenerateNarrative"))
                .narrative(string(molecularMatchObject, "narrative"))
                .expression(string(molecularMatchObject, "expression"))
                .customer(string(molecularMatchObject, "customer"))
                .version(string(molecularMatchObject, "version"))
                .id(string(molecularMatchObject, "id"))
                .externalIds(optionalStringList(molecularMatchObject, "external_id"))
                .uniqueKey(string(molecularMatchObject, "uniqueKey"))
                .hashKey(string(molecularMatchObject, "hashKey"))
                .build();
    }

    @NotNull
    private static List<MolecularMatchMutation> createMutations(@NotNull JsonArray mutationArray) {
        List<MolecularMatchMutation> mutationList = Lists.newArrayList();
        ViccDatamodelChecker mutationChecker = ViccDatamodelCheckerFactory.molecularMatchMutationChecker();

        for (JsonElement mutationElement : mutationArray) {
            JsonObject mutationObject = mutationElement.getAsJsonObject();
            mutationChecker.check(mutationObject);

            mutationList.add(ImmutableMolecularMatchMutation.builder()
                    .geneSymbol(string(mutationObject, "geneSymbol"))
                    .name(string(mutationObject, "name"))
                    .transcriptRecognized(optionalString(mutationObject, "transcriptRecognized"))
                    .transcript(string(mutationObject, "transcript"))
                    .longestTranscript(optionalString(mutationObject, "longestTranscript"))
                    .uniprotTranscript(optionalString(mutationObject, "uniprotTranscript"))
                    .transcriptConsequences(createTranscriptConsequences(optionalJsonArray(mutationObject, "transcriptConsequence")))
                    .parents(createParents(mutationObject.getAsJsonArray("parents")))
                    .wgsaLocations(createWGSALocations(optionalJsonObject(mutationObject, "wgsaData")))
                    .wgsaMaps(createWGSAMaps(optionalJsonArray(mutationObject, "wgsaMap")))
                    .exonsInfo(createExonsInfo(optionalJsonObject(mutationObject, "exonsInfo")))
                    .fusionData(createFusionData(optionalJsonArray(mutationObject, "fusionData")))
                    .mutationTypes(stringList(mutationObject, "mutation_type"))
                    .sources(stringList(mutationObject, "sources"))
                    .synonyms(stringList(mutationObject, "synonyms"))
                    .grch37Locations(createGRCh37Locations(mutationObject.getAsJsonArray("GRCh37_location")))
                    .pathology(stringList(mutationObject, "pathology"))
                    .cDNA(stringList(mutationObject, "cdna"))
                    .description(string(mutationObject, "description"))
                    .src(string(mutationObject, "_src"))
                    .id(string(mutationObject, "id"))
                    .build());
        }
        return mutationList;
    }

    @NotNull
    private static List<MolecularMatchTranscriptConsequence> createTranscriptConsequences(@Nullable JsonArray transcriptConsequenceArray) {
        if (transcriptConsequenceArray == null) {
            return Lists.newArrayList();
        }

        List<MolecularMatchTranscriptConsequence> transcriptConsequenceList = Lists.newArrayList();
        ViccDatamodelChecker transcriptConsequenceChecker = ViccDatamodelCheckerFactory.molecularMatchTranscriptConsequenceChecker();

        for (JsonElement transcriptConsequenceElement : transcriptConsequenceArray) {
            JsonObject transcriptConsequenceObject = transcriptConsequenceElement.getAsJsonObject();
            transcriptConsequenceChecker.check(transcriptConsequenceObject);

            transcriptConsequenceList.add(ImmutableMolecularMatchTranscriptConsequence.builder()
                    .chr(optionalString(transcriptConsequenceObject, "chr"))
                    .start(optionalString(transcriptConsequenceObject, "start"))
                    .stop(optionalString(transcriptConsequenceObject, "stop"))
                    .ref(optionalString(transcriptConsequenceObject, "ref"))
                    .alt(optionalString(transcriptConsequenceObject, "alt"))
                    .referenceGenome(string(transcriptConsequenceObject, "referenceGenome"))
                    .transcript(string(transcriptConsequenceObject, "transcript"))
                    .strand(string(transcriptConsequenceObject, "strand"))
                    .cdna(optionalString(transcriptConsequenceObject, "cdna"))
                    .aminoAcidChange(optionalNullableString(transcriptConsequenceObject, "amino_acid_change"))
                    .intronNumber(nullableString(transcriptConsequenceObject, "intronNumber"))
                    .exonNumbers(optionalStringList(transcriptConsequenceObject, "exonNumber"))
                    .suppress(string(transcriptConsequenceObject, "suppress"))
                    .custom(string(transcriptConsequenceObject, "custom"))
                    .validated(string(transcriptConsequenceObject, "validated"))
                    .compositeKey(string(transcriptConsequenceObject, "compositeKey"))
                    .build());
        }
        return transcriptConsequenceList;
    }

    @NotNull
    private static List<MolecularMatchWGSALocation> createWGSALocations(@Nullable JsonObject wgsaDataObject) {
        if (wgsaDataObject == null) {
            return Lists.newArrayList();
        }
        ViccDatamodelCheckerFactory.molecularMatchWGSADataChecker().check(wgsaDataObject);

        List<MolecularMatchWGSALocation> molecularMatchWGSALocationList = Lists.newArrayList();
        ViccDatamodelChecker wgsaLocationChecker = ViccDatamodelCheckerFactory.molecularMatchWGSALocationChecker();

        for (JsonElement wgsaLocationElement : wgsaDataObject.get("locations").getAsJsonArray()) {
            JsonObject wgsaLocationObject = wgsaLocationElement.getAsJsonObject();
            wgsaLocationChecker.check(wgsaLocationObject);

            molecularMatchWGSALocationList.add(ImmutableMolecularMatchWGSALocation.builder()
                    .genes(stringList(wgsaLocationObject, "Gene"))
                    .chr(string(wgsaLocationObject, "Chr"))
                    .start(string(wgsaLocationObject, "Start"))
                    .end(string(wgsaLocationObject, "End"))
                    .ref(string(wgsaLocationObject, "Ref"))
                    .alt(string(wgsaLocationObject, "Alt"))
                    .chrStartRefAlt(string(wgsaLocationObject, "Chr_Start_Ref_Alt"))
                    .transcript(string(wgsaLocationObject, "Transcript"))
                    .nucleotideChange(string(wgsaLocationObject, "NucleotideChange"))
                    .aa(optionalString(wgsaLocationObject, "AA"))
                    .fullAAs(stringList(wgsaLocationObject, "FullAA"))
                    .exonicFunc(optionalString(wgsaLocationObject, "ExonicFunc"))
                    .popFreqMax(string(wgsaLocationObject, "PopFreqMax"))
                    .clinVarDiseases(optionalStringList(wgsaLocationObject, "ClinVar_DIS"))
                    .clinVarSigs(optionalStringList(wgsaLocationObject, "ClinVar_SIG"))
                    .clinVarStates(optionalStringList(wgsaLocationObject, "ClinVar_STATUS"))
                    .clinVarDbIds(optionalStringList(wgsaLocationObject, "ClinVar_DBID"))
                    .exacAFR(optionalString(wgsaLocationObject, "ExAC_AFR"))
                    .exacAMR(optionalString(wgsaLocationObject, "ExAC_AMR"))
                    .exacEAS(optionalString(wgsaLocationObject, "ExAC_EAS"))
                    .exacFIN(optionalString(wgsaLocationObject, "ExAC_FIN"))
                    .exacNFE(optionalString(wgsaLocationObject, "ExAC_NFE"))
                    .exacSAS(optionalString(wgsaLocationObject, "ExAC_SAS"))
                    .exacFreq(optionalString(wgsaLocationObject, "ExAC_Freq"))
                    .g1000AFR(optionalString(wgsaLocationObject, "1000G_AFR"))
                    .g1000AMR(optionalString(wgsaLocationObject, "1000G_AMR"))
                    .g1000EUR(optionalString(wgsaLocationObject, "1000G_EUR"))
                    .g1000EAS(optionalString(wgsaLocationObject, "1000G_EAS"))
                    .g1000SAS(optionalString(wgsaLocationObject, "1000G_SAS"))
                    .g1000ALL(optionalString(wgsaLocationObject, "1000G_ALL"))
                    .fathmm(string(wgsaLocationObject, "FATHMM"))
                    .fathmmPred(string(wgsaLocationObject, "FATHMM_Pred"))
                    .esp6500siAA(optionalString(wgsaLocationObject, "ESP6500si_AA"))
                    .esp6500siEA(optionalString(wgsaLocationObject, "ESP6500si_EA"))
                    .dbSNP(optionalString(wgsaLocationObject, "dbSNP"))
                    .cosmicId(optionalString(wgsaLocationObject, "COSMIC_ID"))
                    .phyloP46wayPlacental(string(wgsaLocationObject, "phyloP46way_placental"))
                    .phyloP100wayVertebrate(string(wgsaLocationObject, "phyloP100way_vertebrate"))
                    .siPhy29wayLogOdds(string(wgsaLocationObject, "SiPhy_29way_logOdds"))
                    .gwasSNP(optionalString(wgsaLocationObject, "GWAS_SNP"))
                    .gwasDIS(optionalString(wgsaLocationObject, "GWAS_DIS"))
                    .gwasPubmed(optionalString(wgsaLocationObject, "GWAS_PUBMED"))
                    .gerpRS(string(wgsaLocationObject, "GERP++_RS"))
                    .func(string(wgsaLocationObject, "Func"))
                    .wgRna(optionalString(wgsaLocationObject, "wgRna"))
                    .targetScanS(optionalString(wgsaLocationObject, "targetScanS"))
                    .key(string(wgsaLocationObject, "_key"))
                    .build());
        }
        return molecularMatchWGSALocationList;
    }

    @NotNull
    private static List<MolecularMatchTherapeuticContext> createTherapeuticContexts(@NotNull JsonArray arrayTherapeuticContext) {
        List<MolecularMatchTherapeuticContext> therapeuticContextList = Lists.newArrayList();
        for (JsonElement therapeuticContext : arrayTherapeuticContext) {
            Set<String> keysTherapeuticContext = therapeuticContext.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES.contains(keysTherapeuticContext.size())) {
                LOGGER.warn("Found {} in molecular match therapeutic context rather than the expected {}",
                        keysTherapeuticContext.size(),
                        EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES);
                LOGGER.warn(keysTherapeuticContext);
            }

            therapeuticContextList.add(ImmutableMolecularMatchTherapeuticContext.builder()
                    .facet(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .suppress(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .valid(!therapeuticContext.getAsJsonObject().has("valid")
                            ? null
                            : therapeuticContext.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .name(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return therapeuticContextList;
    }

    @NotNull
    private static List<MolecularMatchClassification> createClassifications(@NotNull JsonArray objectClassifications) {
        List<MolecularMatchClassification> classificationList = Lists.newArrayList();
        for (JsonElement classification : objectClassifications) {
            Set<String> keysClassification = classification.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES.contains(keysClassification.size())) {
                LOGGER.warn("Found {} in molecular match classification rather than the expected {}",
                        keysClassification.size(),
                        EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES);
                LOGGER.warn(keysClassification);
            }

            classificationList.add(ImmutableMolecularMatchClassification.builder()
                    .end(!classification.getAsJsonObject().has("End")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("End")))
                    .classification(classification.getAsJsonObject().getAsJsonPrimitive("classification").getAsString())
                    .classificationOverride(
                            !classification.getAsJsonObject().has("classificationOverride") || classification.getAsJsonObject()
                                    .get("classificationOverride")
                                    .isJsonNull()
                                    ? null
                                    : classification.getAsJsonObject().getAsJsonPrimitive("classificationOverride").getAsString())
                    .start(!classification.getAsJsonObject().has("Start")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Start")))
                    .chr(!classification.getAsJsonObject().has("Chr")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Chr")))
                    .geneSymbol(!classification.getAsJsonObject().has("geneSymbol")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("geneSymbol").getAsString())
                    .pathology(!classification.getAsJsonObject().has("pathology")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("pathology")))
                    .ref(!classification.getAsJsonObject().has("Ref")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Ref")))
                    .description(!classification.getAsJsonObject().has("description")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .priority(!classification.getAsJsonObject().has("priority")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .nucleotideChange(!classification.getAsJsonObject().has("NucleotideChange")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("NucleotideChange")))
                    .parents(!classification.getAsJsonObject().has("parents")
                            ? null
                            : createParents(classification.getAsJsonObject().getAsJsonArray("parents")))
                    .expandGeneSearch(!classification.getAsJsonObject().has("expandGeneSearch")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("expandGeneSearch").getAsString())
                    .drugsExperimentalCount(!classification.getAsJsonObject().has("drugsExperimentalCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("drugsExperimentalCount").getAsString())
                    .exon(!classification.getAsJsonObject().has("Exon")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Exon")))
                    .drugsApprovedOffLabelCount(!classification.getAsJsonObject().has("drugsApprovedOffLabelCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("drugsApprovedOffLabelCount").getAsString())
                    .exonicFunc(!classification.getAsJsonObject().has("ExonicFunc")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("ExonicFunc")))
                    .popFreqMax(!classification.getAsJsonObject().has("PopFreqMax")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("PopFreqMax")))
                    .copyNumberType(!classification.getAsJsonObject().has("copyNumberType") || classification.getAsJsonObject()
                            .get("copyNumberType")
                            .isJsonNull() ? null : classification.getAsJsonObject().getAsJsonPrimitive("copyNumberType").getAsString())
                    .publicationCount(!classification.getAsJsonObject().has("publicationCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("publicationCount").getAsString())
                    .transcript(!classification.getAsJsonObject().has("transcript") || classification.getAsJsonObject()
                            .get("transcript")
                            .isJsonNull() ? null : classification.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .dbSNP(!classification.getAsJsonObject().has("dbSNP")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("dbSNP")))
                    .alt(!classification.getAsJsonObject().has("Alt")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Alt")))
                    .name(!classification.getAsJsonObject().has("name")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .rootTerm(!classification.getAsJsonObject().has("rootTerm")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("rootTerm").getAsString())
                    .sources(!classification.getAsJsonObject().has("sources")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("sources")))
                    .drugsApprovedOnLabelCount(!classification.getAsJsonObject().has("drugsApprovedOnLabelCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("drugsApprovedOnLabelCount").getAsString())
                    .trialCount(!classification.getAsJsonObject().has("trialCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("trialCount").getAsString())
                    .alias(!classification.getAsJsonObject().has("alias")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("alias").getAsString())
                    .cosmicId(!classification.getAsJsonObject().has("COSMIC_ID")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("COSMIC_ID")))
                    .transcripts(!classification.getAsJsonObject().has("transcripts")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("transcripts")))
                    .build());
        }
        return classificationList;
    }

    @NotNull
    private static List<MolecularMatchParent> createParents(@NotNull JsonArray arrayParents) {
        List<MolecularMatchParent> parentsList = Lists.newArrayList();
        for (JsonElement parents : arrayParents) {
            Set<String> keysParents = parents.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_PARENTS_SIZES.contains(keysParents.size())) {
                LOGGER.warn("Found {} in molecular match parents rather than the expected {}",
                        keysParents.size(),
                        EXPECTED_MOLECULARMATCH_PARENTS_SIZES);
                LOGGER.warn(keysParents);
            }

            parentsList.add(ImmutableMolecularMatchParent.builder()
                    .transcripts(toStringList(parents.getAsJsonObject().getAsJsonArray("transcripts")))
                    .type(!parents.getAsJsonObject().has("type") || parents.getAsJsonObject().get("type").isJsonNull()
                            ? null
                            : parents.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .name(parents.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .actionableParent(!parents.getAsJsonObject().has("actionableParent")
                            ? null
                            : parents.getAsJsonObject().getAsJsonPrimitive("actionableParent").getAsString())
                    .build());
        }
        return parentsList;
    }

    @NotNull
    private static List<MolecularMatchTag> createTags(@NotNull JsonArray arrayTags) {
        List<MolecularMatchTag> tagsList = Lists.newArrayList();
        for (JsonElement tags : arrayTags) {
            Set<String> keysTags = tags.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TAGS_SIZES.contains(keysTags.size())) {
                LOGGER.warn("Found {} in molecular match tags rather than the expected {}",
                        keysTags.size(),
                        EXPECTED_MOLECULARMATCH_TAGS_SIZES);
                LOGGER.warn(keysTags);
            }

            tagsList.add(ImmutableMolecularMatchTag.builder()
                    .priority(tags.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .compositeKey(!tags.getAsJsonObject().has("compositeKey")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .suppress(!tags.getAsJsonObject().has("suppress")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(!tags.getAsJsonObject().has("filterType")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(tags.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .primary(!tags.getAsJsonObject().has("primary")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("primary").getAsString())
                    .facet(tags.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .valid(!tags.getAsJsonObject().has("valid") ? null : tags.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .custom(!tags.getAsJsonObject().has("custom")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .isNew(!tags.getAsJsonObject().has("isNew") ? null : tags.getAsJsonObject().getAsJsonPrimitive("isNew").getAsString())
                    .generatedBy(!tags.getAsJsonObject().has("generatedBy")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("generatedBy").getAsString())
                    .manualSuppress(!tags.getAsJsonObject().has("manualSuppress")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("manualSuppress").getAsString())
                    .generatedByTerm(!tags.getAsJsonObject().has("generatedByTerm")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("generatedByTerm").getAsString())
                    .transcript(!tags.getAsJsonObject().has("transcript")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .build());
        }
        return tagsList;
    }

    @NotNull
    private static List<MolecularMatchTierExplanation> createTierExplanations(@NotNull JsonArray arrarTierExplanation) {
        List<MolecularMatchTierExplanation> tierExplanationList = Lists.newArrayList();
        for (JsonElement tierExplanation : arrarTierExplanation) {
            Set<String> keysTierExplanation = tierExplanation.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TIER_EXPLANATION_SIZES.contains(keysTierExplanation.size())) {
                LOGGER.warn("Found {} in molecular match tier explanation rather than the expected {}",
                        keysTierExplanation.size(),
                        EXPECTED_MOLECULARMATCH_TIER_EXPLANATION_SIZES);
                LOGGER.warn(keysTierExplanation);
            }

            tierExplanationList.add(ImmutableMolecularMatchTierExplanation.builder()
                    .tier(tierExplanation.getAsJsonObject().getAsJsonPrimitive("tier").getAsString())
                    .step(tierExplanation.getAsJsonObject().getAsJsonPrimitive("step").getAsString())
                    .message(tierExplanation.getAsJsonObject().getAsJsonPrimitive("message").getAsString())
                    .success(tierExplanation.getAsJsonObject().getAsJsonPrimitive("success").getAsString())
                    .build());
        }
        return tierExplanationList;
    }

    @NotNull
    private static List<MolecularMatchVariantInfo> createVariantInfos(@NotNull JsonArray arrayVariantInfo) {
        List<MolecularMatchVariantInfo> variantInfoList = Lists.newArrayList();

        for (JsonElement variantInfo : arrayVariantInfo) {
            Set<String> keysVariantInfo = variantInfo.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_VARIANTINFO_SIZES.contains(keysVariantInfo.size())) {
                LOGGER.warn("Found {} in molecular match variant info rather than the expected {}",
                        keysVariantInfo.size(),
                        EXPECTED_MOLECULARMATCH_VARIANTINFO_SIZES);
                LOGGER.warn(keysVariantInfo);
            }

            variantInfoList.add(ImmutableMolecularMatchVariantInfo.builder()
                    .classification(variantInfo.getAsJsonObject().getAsJsonPrimitive("classification").getAsString())
                    .name(variantInfo.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .consequences(toStringList(variantInfo.getAsJsonObject().getAsJsonArray("consequences")))
                    .fusions(createFusions(variantInfo.getAsJsonObject().getAsJsonArray("fusions")))
                    .locations(createLocations(variantInfo.getAsJsonObject().getAsJsonArray("locations")))
                    .geneFusionPartner(variantInfo.getAsJsonObject().getAsJsonPrimitive("geneFusionPartner").getAsString())
                    .cosmicId(variantInfo.getAsJsonObject().get("COSMIC_ID").isJsonNull()
                            ? null
                            : variantInfo.getAsJsonObject().getAsJsonPrimitive("COSMIC_ID").getAsString())
                    .gene(variantInfo.getAsJsonObject().getAsJsonPrimitive("gene").getAsString())
                    .transcript(variantInfo.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .popFreqMax(variantInfo.getAsJsonObject().getAsJsonPrimitive("popFreqMax").getAsString())
                    .build());
        }
        return variantInfoList;
    }

    @NotNull
    private static List<MolecularMatchFusion> createFusions(@NotNull JsonArray arrayFusions) {
        List<MolecularMatchFusion> fusionsList = Lists.newArrayList();

        for (JsonElement fusions : arrayFusions) {
            Set<String> keysFusions = fusions.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_FUSIONS_SIZES.contains(keysFusions.size())) {
                LOGGER.warn("Found {} in molecular match fusions rather than the expected {}",
                        keysFusions.size(),
                        EXPECTED_MOLECULARMATCH_FUSIONS_SIZES);
                LOGGER.warn(keysFusions);
            }
            fusionsList.add(ImmutableMolecularMatchFusion.builder()
                    .referenceGenome(fusions.getAsJsonObject().get("referenceGenome").getAsString())
                    .LBPWREP(fusions.getAsJsonObject().get("LBPWREP").getAsString())
                    .RBPWREP(fusions.getAsJsonObject().get("RBPWREP").getAsString())
                    .exonNumber(fusions.getAsJsonObject().get("exonNumber").getAsString())
                    .chr(fusions.getAsJsonObject().get("chr").getAsString())
                    .RBPWLEP(fusions.getAsJsonObject().get("RBPWLEP").getAsString())
                    .intronNumber(fusions.getAsJsonObject().get("intronNumber").getAsString())
                    .LBPWLEP(fusions.getAsJsonObject().get("LBPWLEP").getAsString())
                    .build());
        }
        return fusionsList;
    }

    @NotNull
    private static List<MolecularMatchLocation> createLocations(@NotNull JsonArray arrayLocations) {
        List<MolecularMatchLocation> locationsList = Lists.newArrayList();
        for (JsonElement locations : arrayLocations) {
            Set<String> keysLocations = locations.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_LOCATIONS_SIZES.contains(keysLocations.size())) {
                LOGGER.warn("Found {} in molecular match locations rather than the expected {}",
                        keysLocations.size(),
                        EXPECTED_MOLECULARMATCH_LOCATIONS_SIZES);
                LOGGER.warn(keysLocations);
            }

            locationsList.add(ImmutableMolecularMatchLocation.builder()
                    .aminoAcidChange(!locations.getAsJsonObject().has("amino_acid_change")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("amino_acid_change").getAsString())
                    .intronNumber(!locations.getAsJsonObject().has("intronNumber")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .exonNumbers(!locations.getAsJsonObject().has("exonNumber") ? null : createArrayExonNumber(locations))
                    .stop(locations.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .start(locations.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(locations.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .strand(!locations.getAsJsonObject().has("strand")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .alt(!locations.getAsJsonObject().has("alt")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .referenceGenome(!locations.getAsJsonObject().has("referenceGenome")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("referenceGenome").getAsString())
                    .ref(!locations.getAsJsonObject().has("ref")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .cdna(!locations.getAsJsonObject().has("cdna")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .build());
        }
        return locationsList;
    }

    @NotNull
    private static List<String> createArrayExonNumber(@NotNull JsonElement elementExonNumber) {
        if (elementExonNumber.getAsJsonObject().get("exonNumber").isJsonArray()) {
            return toStringList(elementExonNumber.getAsJsonObject().getAsJsonArray("exonNumber"));
        } else {
            return Collections.singletonList(elementExonNumber.getAsJsonObject().getAsJsonPrimitive("exonNumber").getAsString());
        }
    }

    @NotNull
    private static MolecularMatchAst createAst(@NotNull JsonObject objectAst) {
        Set<String> keysAst = objectAst.keySet();
        if (!EXPECTED_MOLECULARMATCH_AST_SIZES.contains(keysAst.size())) {
            LOGGER.warn("Found {} in molecular match ast rather than the expected {}", keysAst.size(), EXPECTED_MOLECULARMATCH_AST_SIZES);
            LOGGER.warn(keysAst);
        }
        return ImmutableMolecularMatchAst.builder()
                .raw(objectAst.get("raw") == null ? null : objectAst.getAsJsonPrimitive("raw").getAsString())
                .value(!objectAst.has("value") ? null : objectAst.getAsJsonPrimitive("value").getAsString())
                .operator(!objectAst.has("operator") ? null : objectAst.getAsJsonPrimitive("operator").getAsString())
                .right(!objectAst.has("right") ? null : createRight(objectAst.getAsJsonObject("right")))
                .type(objectAst.getAsJsonPrimitive("type").getAsString())
                .left(!objectAst.has("left") ? null : createLeft(objectAst.getAsJsonObject("left")))
                .build();
    }

    @NotNull
    private static MolecularMatchAstLeft createLeft(@NotNull JsonObject objectLeft) {
        Set<String> keysLeft = objectLeft.keySet();
        if (!EXPECTED_MOLECULARMATCH_LEFT_SIZES.contains(keysLeft.size())) {
            LOGGER.warn("Found {} in molecular match ast left rather than the expected {}",
                    keysLeft.size(),
                    EXPECTED_MOLECULARMATCH_LEFT_SIZES);
            LOGGER.warn(keysLeft);
        }

        return ImmutableMolecularMatchAstLeft.builder()
                .operator(!objectLeft.has("operator") ? null : objectLeft.getAsJsonPrimitive("operator").getAsString())
                .raw(!objectLeft.has("raw") ? null : objectLeft.getAsJsonPrimitive("raw").getAsString())
                .type(objectLeft.getAsJsonPrimitive("type").getAsString())
                .value(!objectLeft.has("value") ? null : objectLeft.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRight createRight(@NotNull JsonObject objectRight) {
        Set<String> keysRight = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_SIZES.contains(keysRight.size())) {
            LOGGER.warn("Found {} in molecular match ast right rather than the expected {}",
                    keysRight.size(),
                    EXPECTED_MOLECULARMATCH_RIGHT_SIZES);
            LOGGER.warn(keysRight);
        }

        return ImmutableMolecularMatchAstRight.builder()
                .operator(!objectRight.has("operator") ? null : objectRight.getAsJsonPrimitive("operator").getAsString())
                .left(!objectRight.has("left") ? null : createRightLeft(objectRight.getAsJsonObject("left")))
                .right(!objectRight.has("right") ? null : createRightRight(objectRight.getAsJsonObject("right")))
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRightRight createRightRight(@NotNull JsonObject objectRight) {
        Set<String> keysRightRight = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES.contains(keysRightRight.size())) {
            LOGGER.warn("Found {} in molecular match ast right right rather than the expected {}",
                    keysRightRight.size(),
                    EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES);
            LOGGER.warn(keysRightRight);
        }

        return ImmutableMolecularMatchAstRightRight.builder()
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRightLeft createRightLeft(@NotNull JsonObject objectRight) {
        Set<String> keysRightLeft = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES.contains(keysRightLeft.size())) {
            LOGGER.warn("Found {} in molecular match ast right left rather than the expected {}",
                    keysRightLeft.size(),
                    EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES);
            LOGGER.warn(keysRightLeft);
        }

        return ImmutableMolecularMatchAstRightLeft.builder()
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .right(!objectRight.has("right") ? null : createRightLeftRight(objectRight.getAsJsonObject("right")))
                .left(!objectRight.has("left") ? null : createLeft(objectRight.getAsJsonObject("left")))
                .build();
    }

    @NotNull
    private static MolecularMatchAstRightLeftRight createRightLeftRight(@NotNull JsonObject objectRight) {
        Set<String> keysLeftRight = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_LEFT_RIGHT_SIZES.contains(keysLeftRight.size())) {
            LOGGER.warn("Found {} in molecular match ast right left right rather than the expected {}",
                    keysLeftRight.size(),
                    EXPECTED_MOLECULARMATCH_RIGHT_LEFT_RIGHT_SIZES);
            LOGGER.warn(keysLeftRight);
        }

        return ImmutableMolecularMatchAstRightLeftRight.builder()
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static List<MolecularMatchSource> createSources(@NotNull JsonArray arraySources) {
        List<MolecularMatchSource> sourcesList = Lists.newArrayList();
        for (JsonElement source : arraySources) {
            Set<String> keysSource = source.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_SOURCE_SIZES.contains(keysSource.size())) {
                LOGGER.warn("Found {} in molecular match source rather than the expected {}",
                        keysSource.size(),
                        EXPECTED_MOLECULARMATCH_SOURCE_SIZES);
                LOGGER.warn(keysSource);
            }

            sourcesList.add(ImmutableMolecularMatchSource.builder()
                    .name(source.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .suppress(source.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .pubId(source.getAsJsonObject().getAsJsonPrimitive("pubId").getAsString())
                    .subType(!source.getAsJsonObject().has("subType")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("subType").getAsString())
                    .valid(source.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .link(source.getAsJsonObject().getAsJsonPrimitive("link").getAsString())
                    .year(source.getAsJsonObject().getAsJsonPrimitive("year").getAsString())
                    .trialId(!source.getAsJsonObject().has("trialId")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("trialId").getAsString())
                    .type(source.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .id(source.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .institution(!source.getAsJsonObject().has("institution")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("institution").getAsString())
                    .trialPhase(!source.getAsJsonObject().has("trialPhase")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("trialPhase").getAsString())
                    .functionalConsequence(!source.getAsJsonObject().has("functionalConsequence")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("functionalConsequence").getAsString())
                    .trustRating(!source.getAsJsonObject().has("trustRating")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("trustRating").getAsString())
                    .build());
        }
        return sourcesList;
    }

    @NotNull
    private static List<MolecularMatchCriteriaUnmet> createCriteriaUnmets(@NotNull JsonArray arrayCriteriaUnmet) {
        List<MolecularMatchCriteriaUnmet> criteriaUnmetList = Lists.newArrayList();

        for (JsonElement criteriaUnmet : arrayCriteriaUnmet) {
            Set<String> keysCriteriaUnmet = criteriaUnmet.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_CRITERIA_UNMET_SIZES.contains(keysCriteriaUnmet.size())) {
                LOGGER.warn("Found {} in molecular match criteria unmet rather than the expected {}",
                        keysCriteriaUnmet.size(),
                        EXPECTED_MOLECULARMATCH_CRITERIA_UNMET_SIZES);
                LOGGER.warn(keysCriteriaUnmet);
            }

            criteriaUnmetList.add(ImmutableMolecularMatchCriteriaUnmet.builder()
                    .priority(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .compositeKey(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .isNew(!criteriaUnmet.getAsJsonObject().has("isNew")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("isNew").getAsString())
                    .generatedBy(!criteriaUnmet.getAsJsonObject().has("generatedBy")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("generatedBy").getAsString())
                    .manualSuppress(!criteriaUnmet.getAsJsonObject().has("manualSuppress")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("manualSuppress").getAsString())
                    .generatedByTerm(!criteriaUnmet.getAsJsonObject().has("generatedByTerm")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("generatedByTerm").getAsString())
                    .suppress(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .primary(!criteriaUnmet.getAsJsonObject().has("primary")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("primary").getAsString())
                    .facet(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .valid(!criteriaUnmet.getAsJsonObject().has("valid")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .custom(!criteriaUnmet.getAsJsonObject().has("custom")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .transcript(!criteriaUnmet.getAsJsonObject().has("transcript")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .build());
        }
        return criteriaUnmetList;
    }

    @NotNull
    private static List<MolecularMatchPrevalence> createPrevalences(@NotNull JsonArray arrayPrevelance) {
        List<MolecularMatchPrevalence> prevalenceList = Lists.newArrayList();

        for (JsonElement prevalence : arrayPrevelance) {
            Set<String> keysPrevalence = prevalence.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_PREVALENCE_SIZES.contains(keysPrevalence.size())) {
                LOGGER.warn("Found {} in molecular match prevalence rather than the expected {}",
                        keysPrevalence.size(),
                        EXPECTED_MOLECULARMATCH_PREVALENCE_SIZES);
                LOGGER.warn(keysPrevalence);
            }

            prevalenceList.add(ImmutableMolecularMatchPrevalence.builder()
                    .count(prevalence.getAsJsonObject().getAsJsonPrimitive("count").getAsString())
                    .percent(prevalence.getAsJsonObject().getAsJsonPrimitive("percent").getAsString())
                    .studyId(prevalence.getAsJsonObject().getAsJsonPrimitive("studyId").getAsString())
                    .samples(prevalence.getAsJsonObject().getAsJsonPrimitive("samples").getAsString())
                    .molecular(!prevalence.getAsJsonObject().has("molecular")
                            ? null
                            : prevalence.getAsJsonObject().getAsJsonPrimitive("molecular").getAsString())
                    .condition(!prevalence.getAsJsonObject().has("condition")
                            ? null
                            : prevalence.getAsJsonObject().getAsJsonPrimitive("condition").getAsString())
                    .build());
        }
        return prevalenceList;
    }

    @NotNull
    private static List<MolecularMatchFusionData> createFusionData(@Nullable JsonArray arrayFusionData) {
        if (arrayFusionData == null) {
            return Lists.newArrayList();
        }

        List<MolecularMatchFusionData> fusionDataList = Lists.newArrayList();
        for (JsonElement fusionData : arrayFusionData.getAsJsonArray()) {
            Set<String> keysFusionData = fusionData.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_FUSIONDATA_SIZES.contains(keysFusionData.size())) {
                LOGGER.warn("Found {} in molecular match fusion data rather than the expected {}",
                        keysFusionData.size(),
                        EXPECTED_MOLECULARMATCH_FUSIONDATA_SIZES);
                LOGGER.warn(keysFusionData);
            }

            fusionDataList.add(ImmutableMolecularMatchFusionData.builder()
                    .Bgreg(!fusionData.getAsJsonObject().has("Bgreg")
                            ? null
                            : createBreg(fusionData.getAsJsonObject().getAsJsonArray("Bgreg")))
                    .Bchr(!fusionData.getAsJsonObject().has("Bchr")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bchr")))
                    .synonym(!fusionData.getAsJsonObject().has("synonym")
                            ? null
                            : fusionData.getAsJsonObject().getAsJsonPrimitive("synonym").getAsString())
                    .Agene(!fusionData.getAsJsonObject().has("Agene")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Agene")))
                    .Btx(!fusionData.getAsJsonObject().has("Btx") ? null : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Btx")))
                    .Achr(!fusionData.getAsJsonObject().has("Achr")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Achr")))
                    .ins(!fusionData.getAsJsonObject().has("ins") ? null : toStringList(fusionData.getAsJsonObject().getAsJsonArray("ins")))
                    .source(!fusionData.getAsJsonObject().has("source")
                            ? null
                            : fusionData.getAsJsonObject().getAsJsonPrimitive("source").getAsString())
                    .Agreg(!fusionData.getAsJsonObject().has("Agreg")
                            ? null
                            : createAreg(fusionData.getAsJsonObject().getAsJsonArray("Agreg")))
                    .Bgene(!fusionData.getAsJsonObject().has("Bgene")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bgene")))
                    .Acoord(!fusionData.getAsJsonObject().has("Acoord")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Acoord")))
                    .Bori(!fusionData.getAsJsonObject().has("Bori")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bori")))
                    .Aband(!fusionData.getAsJsonObject().has("Aband")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Aband")))
                    .Bband(!fusionData.getAsJsonObject().has("Bband")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bband")))
                    .Aori(!fusionData.getAsJsonObject().has("Aori")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Aori")))
                    .Atx(!fusionData.getAsJsonObject().has("Atx") ? null : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Atx")))
                    .Bcoord(!fusionData.getAsJsonObject().has("Bcoord")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bcoord")))
                    .Paper(!fusionData.getAsJsonObject().has("Paper")
                            ? null
                            : fusionData.getAsJsonObject().getAsJsonPrimitive("Paper").getAsString())
                    .build());
        }
        return fusionDataList;
    }

    @NotNull
    private static List<MolecularMatchAgreg> createAreg(@NotNull JsonArray arrayAreg) {
        List<MolecularMatchAgreg> aregList = Lists.newArrayList();
        for (JsonElement areg : arrayAreg.getAsJsonArray()) {
            Set<String> keysAreg = areg.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_AREG_SIZES.contains(keysAreg.size())) {
                LOGGER.warn("Found {} in molecular match areg rather than the expected {}",
                        keysAreg.size(),
                        EXPECTED_MOLECULARMATCH_AREG_SIZES);
                LOGGER.warn(keysAreg);
            }

            aregList.add(ImmutableMolecularMatchAgreg.builder()
                    .type(areg.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .num(areg.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .build());
        }
        return aregList;
    }

    @NotNull
    private static List<MolecularMatchBreg> createBreg(@NotNull JsonArray arrayBreg) {
        List<MolecularMatchBreg> bregList = Lists.newArrayList();
        for (JsonElement breg : arrayBreg.getAsJsonArray()) {
            Set<String> keysBreg = breg.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_BREG_SIZES.contains(keysBreg.size())) {
                LOGGER.warn("Found {} in molecular match breg rather than the expected {}",
                        keysBreg.size(),
                        EXPECTED_MOLECULARMATCH_BREG_SIZES);
                LOGGER.warn(keysBreg);
            }

            bregList.add(ImmutableMolecularMatchBreg.builder()
                    .type(breg.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .num(breg.getAsJsonObject().getAsJsonPrimitive("num").getAsString())
                    .build());
        }
        return bregList;
    }

    @Nullable
    private static MolecularMatchExonsInfo createExonsInfo(@Nullable JsonObject objectExonsInfo) {
        if (objectExonsInfo == null) {
            return null;
        }
        Set<String> keysExonsInfo = objectExonsInfo.keySet();
        if (!EXPECTED_MOLECULARMATCH_EXONSINFO_SIZES.contains(keysExonsInfo.size())) {
            LOGGER.warn("Found {} in molecular match exon info rather than the expected {}",
                    keysExonsInfo.size(),
                    EXPECTED_MOLECULARMATCH_EXONSINFO_SIZES);
            LOGGER.warn(keysExonsInfo);
        }

        return ImmutableMolecularMatchExonsInfo.builder()
                .exonBoundaries(createExonBoundaries(objectExonsInfo.getAsJsonObject("exonBoundaries")))
                .txStart(!objectExonsInfo.has("txStart") ? null : objectExonsInfo.getAsJsonPrimitive("txStart").getAsString())
                .cdsEnd(!objectExonsInfo.has("cdsEnd") ? null : objectExonsInfo.getAsJsonPrimitive("cdsEnd").getAsString())
                .chr(objectExonsInfo.getAsJsonPrimitive("chr").getAsString())
                .cdsStart(!objectExonsInfo.has("cdsEnd") ? null : objectExonsInfo.getAsJsonPrimitive("cdsEnd").getAsString())
                .transcript(objectExonsInfo.getAsJsonPrimitive("transcript").getAsString())
                .txEnd(!objectExonsInfo.has("txEnd") ? null : objectExonsInfo.getAsJsonPrimitive("txEnd").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchPosition createMolecularPositions(@NotNull JsonObject objectPositions) {
        Set<String> keysPositions = objectPositions.keySet();
        if (!EXPECTED_MOLECULARMATCH_POSITIONS_SIZES.contains(keysPositions.size())) {
            LOGGER.warn("Found {} in molecular match positions rather than the expected {}",
                    keysPositions.size(),
                    EXPECTED_MOLECULARMATCH_POSITIONS_SIZES);
            LOGGER.warn(keysPositions);
        }

        return ImmutableMolecularMatchPosition.builder()
                .start(objectPositions.getAsJsonPrimitive("start").getAsString())
                .stop(objectPositions.getAsJsonPrimitive("stop").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchExonBoundaries createExonBoundaries(@NotNull JsonObject objectExonsBoundaries) {
        Set<String> keysExonsBoundaries = objectExonsBoundaries.keySet();
        if (!EXPECTED_MOLECULARMATCH_EXONSBOUNDARIES_SIZES.contains(keysExonsBoundaries.size())) {
            LOGGER.warn("Found {} in molecular match exons boundaries rather than the expected {}",
                    keysExonsBoundaries.size(),
                    EXPECTED_MOLECULARMATCH_EXONSBOUNDARIES_SIZES);
            LOGGER.warn(keysExonsBoundaries);
        }

        return ImmutableMolecularMatchExonBoundaries.builder()
                .exon1(!objectExonsBoundaries.has("1") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("1")))
                .exon2(!objectExonsBoundaries.has("2") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("2")))
                .exon3(!objectExonsBoundaries.has("3") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("3")))
                .exon4(!objectExonsBoundaries.has("4") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("4")))
                .exon5(!objectExonsBoundaries.has("5") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("5")))
                .exon6(!objectExonsBoundaries.has("6") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("6")))
                .exon7(!objectExonsBoundaries.has("7") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("7")))
                .exon8(!objectExonsBoundaries.has("8") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("8")))
                .exon9(!objectExonsBoundaries.has("9") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("9")))
                .exon10(!objectExonsBoundaries.has("10") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("10")))
                .exon11(!objectExonsBoundaries.has("11") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("11")))
                .exon12(!objectExonsBoundaries.has("12") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("12")))
                .exon13(!objectExonsBoundaries.has("13") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("13")))
                .exon14(!objectExonsBoundaries.has("14") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("14")))
                .exon15(!objectExonsBoundaries.has("15") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("15")))
                .exon16(!objectExonsBoundaries.has("16") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("16")))
                .exon17(!objectExonsBoundaries.has("17") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("17")))
                .exon18(!objectExonsBoundaries.has("18") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("18")))
                .exon19(!objectExonsBoundaries.has("19") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("19")))
                .exon20(!objectExonsBoundaries.has("20") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("20")))
                .exon21(!objectExonsBoundaries.has("21") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("21")))
                .exon22(!objectExonsBoundaries.has("22") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("22")))
                .exon23(!objectExonsBoundaries.has("23") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("23")))
                .exon24(!objectExonsBoundaries.has("24") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("24")))
                .exon25(!objectExonsBoundaries.has("25") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("25")))
                .exon26(!objectExonsBoundaries.has("26") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("26")))
                .exon27(!objectExonsBoundaries.has("27") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("27")))
                .exon28(!objectExonsBoundaries.has("28") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("28")))
                .exon29(!objectExonsBoundaries.has("29") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("29")))
                .exon30(!objectExonsBoundaries.has("30") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("30")))
                .exon31(!objectExonsBoundaries.has("31") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("31")))
                .exon32(!objectExonsBoundaries.has("32") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("32")))
                .exon33(!objectExonsBoundaries.has("33") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("33")))
                .exon34(!objectExonsBoundaries.has("34") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("34")))
                .exon35(!objectExonsBoundaries.has("35") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("35")))
                .exon36(!objectExonsBoundaries.has("36") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("36")))
                .exon37(!objectExonsBoundaries.has("37") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("37")))
                .exon38(!objectExonsBoundaries.has("38") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("38")))
                .exon39(!objectExonsBoundaries.has("39") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("39")))
                .exon40(!objectExonsBoundaries.has("40") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("40")))
                .exon41(!objectExonsBoundaries.has("41") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("41")))
                .build();
    }

    @NotNull
    private static List<MolecularMatchWGSAMap> createWGSAMaps(@Nullable JsonArray objectWgsaMap) {
        if (objectWgsaMap == null) {
            return Lists.newArrayList();

        }
        List<MolecularMatchWGSAMap> molecularMatchWGSaMapList = Lists.newArrayList();
        for (JsonElement wgsDataMap : objectWgsaMap.getAsJsonArray()) {
            Set<String> keysWgsaMap = wgsDataMap.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_WGSAMAP_SIZES.contains(keysWgsaMap.size())) {
                LOGGER.warn("Found {} in molecular match wgsa map rather than the expected {}",
                        keysWgsaMap.size(),
                        EXPECTED_MOLECULARMATCH_WGSAMAP_SIZES);
                LOGGER.warn(keysWgsaMap);
            }

            molecularMatchWGSaMapList.add(ImmutableMolecularMatchWGSAMap.builder()
                    .aa(!wgsDataMap.getAsJsonObject().has("AA")
                            ? null
                            : wgsDataMap.getAsJsonObject().getAsJsonPrimitive("AA").getAsString())
                    .name(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .grch37ChrStartRefAlt(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("GRCh37_Chr_Start_Ref_Alt").getAsString())
                    .synonyms(toStringList(wgsDataMap.getAsJsonObject().getAsJsonArray("Synonyms")))
                    .protCoords(toStringList(wgsDataMap.getAsJsonObject().getAsJsonArray("ProtCoords")))
                    .nucleotideChange(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("NucleotideChange").getAsString())
                    .exon(!wgsDataMap.getAsJsonObject().has("Exon")
                            ? null
                            : wgsDataMap.getAsJsonObject().getAsJsonPrimitive("Exon").getAsString())
                    .gene(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("Gene").getAsString())
                    .transcript(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("Transcript").getAsString())
                    .build());
        }
        return molecularMatchWGSaMapList;
    }

    @NotNull
    private static List<MolecularMatchGRCh37Location> createGRCh37Locations(@NotNull JsonArray arrayLocation) {
        List<MolecularMatchGRCh37Location> grch37LocationList = Lists.newArrayList();

        for (JsonElement location : arrayLocation) {
            Set<String> keysLocation = location.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_LOCATIONGRCH37_SIZES.contains(keysLocation.size())) {
                LOGGER.warn("Found {} in molecular match grch37 location rather than the expected {}",
                        keysLocation.size(),
                        EXPECTED_MOLECULARMATCH_LOCATIONGRCH37_SIZES);
                LOGGER.warn(keysLocation);
            }

            grch37LocationList.add(ImmutableMolecularMatchGRCh37Location.builder()
                    .compositeKey(location.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .ref(location.getAsJsonObject().get("ref").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .stop(location.getAsJsonObject().get("stop").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .start(location.getAsJsonObject().get("start").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(location.getAsJsonObject().get("chr").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .alt(location.getAsJsonObject().get("alt").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .validated(location.getAsJsonObject().getAsJsonPrimitive("validated").getAsString())
                    .transcriptConsequences(createConsequencesGRCH37(location.getAsJsonObject().getAsJsonArray("transcript_consequences")))
                    .strand(location.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .build());
        }
        return grch37LocationList;
    }

    @NotNull
    private static List<MolecularMatchTranscriptConsequencesGRCh37> createConsequencesGRCH37(
            @NotNull JsonArray arrayTranscriptConsequence) {
        List<MolecularMatchTranscriptConsequencesGRCh37> transcriptConsequencesGRCH37List = Lists.newArrayList();
        for (JsonElement transcriptConsequences : arrayTranscriptConsequence) {
            Set<String> keysTranscriptConsequences = transcriptConsequences.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES__GRCH37_SIZES.contains(keysTranscriptConsequences.size())) {
                LOGGER.warn("Found {} in molecular match transcript consequences grch37 rather than the expected {}",
                        keysTranscriptConsequences.size(),
                        EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES__GRCH37_SIZES);
                LOGGER.warn(keysTranscriptConsequences);
            }

            transcriptConsequencesGRCH37List.add(ImmutableMolecularMatchTranscriptConsequencesGRCh37.builder()
                    .aminoAcidChange(transcriptConsequences.getAsJsonObject().get("amino_acid_change").isJsonNull()
                            ? null
                            : transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("amino_acid_change").getAsString())
                    .txSites(toStringList(transcriptConsequences.getAsJsonObject().getAsJsonArray("txSites")))
                    .exonNumbers(transcriptConsequences.getAsJsonObject().get("exonNumber").isJsonNull()
                            ? null
                            : createArrayExonNumber(transcriptConsequences))
                    .intronNumber(transcriptConsequences.getAsJsonObject().get("intronNumber").isJsonNull()
                            ? null
                            : transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .transcript(transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .cdna(transcriptConsequences.getAsJsonObject().get("cdna").isJsonNull()
                            ? null
                            : transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .build());
        }
        return transcriptConsequencesGRCH37List;
    }
}
