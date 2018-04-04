package com.hartwig.hmftools.patientdb.readers.cpct;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Map;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableBaselineData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BaselineReader {
    private static final Logger LOGGER = LogManager.getLogger(BaselineReader.class);

    private static final String STUDY_BASELINE = "SE.BASELINE";
    private static final String STUDY_ENDSTUDY = "SE.ENDSTUDY";

    private static final String FORM_DEMOGRAPHY = "FRM.DEMOGRAPHY";
    private static final String FORM_INFORMED_CONSENT = "FRM.INFORMEDCONSENT";
    private static final String FORM_CARCINOMA = "FRM.CARCINOMA";
    private static final String FORM_ELIGIBILITY = "FRM.ELIGIBILITY";
    private static final String FORM_SELCRIT = "FRM.SELCRIT";
    private static final String FORM_DEATH = "FRM.DEATH";

    private static final String ITEMGROUP_DEMOGRAPHY = "GRP.DEMOGRAPHY.DEMOGRAPHY";
    private static final String ITEMGROUP_INFORMED_CONSENT = "GRP.INFORMEDCONSENT.INFORMEDCONSENT";
    private static final String ITEMGROUP_CARCINOMA = "GRP.CARCINOMA.CARCINOMA";
    private static final String ITEMGROUP_ELIGIBILITY = "GRP.ELIGIBILITY.ELIGIBILITY";
    private static final String ITEMGROUP_SELCRIT = "GRP.SELCRIT.SELCRIT";
    private static final String ITEMGROUP_DEATH = "GRP.DEATH.DEATH";

    public static final String FIELD_SEX = "FLD.DEMOGRAPHY.SEX";
    public static final String FIELD_INFORMED_CONSENT_DATE = "FLD.INFORMEDCONSENT.ICDTC";
    public static final String FIELD_REGISTRATION_DATE1 = "FLD.ELIGIBILITY.REGDTC";
    public static final String FIELD_REGISTRATION_DATE2 = "FLD.SELCRIT.NREGDTC";
    public static final String FIELD_BIRTH_YEAR1 = "FLD.SELCRIT.NBIRTHYEAR";
    public static final String FIELD_BIRTH_YEAR2 = "FLD.ELIGIBILITY.BIRTHYEAR";
    public static final String FIELD_BIRTH_YEAR3 = "FLD.ELIGIBILITY.BIRTHDTCES";

    public static final String FIELD_PRIMARY_TUMOR_LOCATION = "FLD.CARCINOMA.PTUMLOC";
    public static final String FIELD_PRIMARY_TUMOR_LOCATION_OTHER = "FLD.CARCINOMA.PTUMLOCS";

    public static final String FIELD_DEATH_DATE = "FLD.DEATH.DDEATHDTC";

    public static final String FIELD_HOSPITAL1 = "FLD.ELIGIBILITY.HOSPITAL";
    public static final String FIELD_HOSPITAL2 = "FLD.SELCRIT.NHOSPITAL";

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final TumorLocationCurator tumorLocationCurator;
    @NotNull
    private final Map<Integer, String> hospitals;

    BaselineReader(@NotNull TumorLocationCurator tumorLocationCurator, @NotNull Map<Integer, String> hospitals) {
        this.tumorLocationCurator = tumorLocationCurator;
        this.hospitals = hospitals;
    }

    @NotNull
    BaselineData read(@NotNull final EcrfPatient patient) {
        final ImmutableBaselineData.Builder baselineBuilder = ImmutableBaselineData.builder().hospital(getHospital(patient, hospitals));

        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BASELINE)) {
            setDemographyData(baselineBuilder, studyEvent);
            setPrimaryTumorData(baselineBuilder, studyEvent);
            setRegistrationAndBirthData(baselineBuilder, studyEvent);
            setInformedConsent(baselineBuilder, studyEvent);
        }
        setDeathData(baselineBuilder, patient);
        return baselineBuilder.build();
    }

    @Nullable
    private static String getHospital(@NotNull final EcrfPatient patient, @NotNull final Map<Integer, String> hospitals) {
        final Integer hospitalCode = Integer.parseInt(patient.patientId().substring(6, 8));
        final String hospital = hospitals.get(hospitalCode);
        if (hospital == null) {
            LOGGER.warn(FIELD_HOSPITAL1 + ", " + FIELD_HOSPITAL2 + " contained no Hospital with code " + hospitalCode);
        }
        return hospital;
    }

    private void setDemographyData(@NotNull final ImmutableBaselineData.Builder builder, @NotNull final EcrfStudyEvent studyEvent) {
        boolean hasDemographyForm = false;
        for (final EcrfForm demographyForm : studyEvent.nonEmptyFormsPerOID(FORM_DEMOGRAPHY)) {
            for (final EcrfItemGroup demographyItemGroup : demographyForm.nonEmptyItemGroupsPerOID(ITEMGROUP_DEMOGRAPHY)) {
                hasDemographyForm = true;
                builder.gender(demographyItemGroup.readItemString(FIELD_SEX, 0, false));
                builder.demographyStatus(demographyForm.status());
                builder.demographyLocked(demographyForm.locked());
            }
        }

        if (!hasDemographyForm) {
            builder.demographyStatus(FormStatusState.UNKNOWN);
            builder.demographyLocked(false);
        }
    }

    private void setPrimaryTumorData(@NotNull final ImmutableBaselineData.Builder builder, @NotNull final EcrfStudyEvent studyEvent) {
        boolean hasPrimaryTumorForm = false;
        for (final EcrfForm carcinomaForm : studyEvent.nonEmptyFormsPerOID(FORM_CARCINOMA)) {
            for (final EcrfItemGroup carcinomaItemGroup : carcinomaForm.nonEmptyItemGroupsPerOID(ITEMGROUP_CARCINOMA)) {
                hasPrimaryTumorForm = true;
                String primaryTumorLocation = carcinomaItemGroup.readItemString(FIELD_PRIMARY_TUMOR_LOCATION, 0);
                if (primaryTumorLocation != null && primaryTumorLocation.trim().toLowerCase().startsWith("other")) {
                    primaryTumorLocation = carcinomaItemGroup.readItemString(FIELD_PRIMARY_TUMOR_LOCATION_OTHER, 0);
                }
                builder.cancerType(tumorLocationCurator.search(primaryTumorLocation));
                builder.primaryTumorStatus(carcinomaForm.status());
                builder.primaryTumorLocked(carcinomaForm.locked());
            }
        }

        if (!hasPrimaryTumorForm) {
            builder.primaryTumorStatus(FormStatusState.UNKNOWN);
            builder.primaryTumorLocked(false);
        }
    }

    private static void setRegistrationAndBirthData(@NotNull final ImmutableBaselineData.Builder builder,
            @NotNull final EcrfStudyEvent studyEvent) {
        LocalDate registrationDate1 = null;
        LocalDate registrationDate2 = null;
        String birthYear1 = null;
        String birthYear2 = null;
        LocalDate birthYear3 = null;

        boolean hasEligibilityForm = false;
        for (final EcrfForm eligibilityForm : studyEvent.nonEmptyFormsPerOID(FORM_ELIGIBILITY)) {
            for (final EcrfItemGroup eligibilityItemGroup : eligibilityForm.nonEmptyItemGroupsPerOID(ITEMGROUP_ELIGIBILITY)) {
                hasEligibilityForm = true;
                registrationDate1 = eligibilityItemGroup.readItemDate(FIELD_REGISTRATION_DATE1, 0, DATE_FORMATTER);
                birthYear2 = eligibilityItemGroup.readItemString(FIELD_BIRTH_YEAR2, 0);
                birthYear3 = eligibilityItemGroup.readItemDate(FIELD_BIRTH_YEAR3, 0, DATE_FORMATTER);
                builder.eligibilityStatus(eligibilityForm.status());
                builder.eligibilityLocked(eligibilityForm.locked());
            }
        }

        if (!hasEligibilityForm) {
            builder.eligibilityStatus(FormStatusState.UNKNOWN);
            builder.eligibilityLocked(false);
        }

        boolean hasSelectionCriteriaForm = false;
        for (final EcrfForm selcritForm : studyEvent.nonEmptyFormsPerOID(FORM_SELCRIT)) {
            for (final EcrfItemGroup selcritItemGroup : selcritForm.nonEmptyItemGroupsPerOID(ITEMGROUP_SELCRIT)) {
                hasSelectionCriteriaForm = true;
                birthYear1 = selcritItemGroup.readItemString(FIELD_BIRTH_YEAR1, 0);
                if (registrationDate1 == null) {
                    registrationDate2 = selcritItemGroup.readItemDate(FIELD_REGISTRATION_DATE2, 0, DATE_FORMATTER);
                    builder.selectionCriteriaStatus(selcritForm.status());
                    builder.selectionCriteriaLocked(selcritForm.locked());
                }
            }
        }

        if (!hasSelectionCriteriaForm) {
            builder.selectionCriteriaStatus(FormStatusState.UNKNOWN);
            builder.selectionCriteriaLocked(false);
        }

        final LocalDate registrationDate = registrationDate2 == null ? registrationDate1 : registrationDate2;
        final Integer birthYear = determineBirthYear(birthYear1, birthYear2, birthYear3);
        builder.registrationDate(registrationDate);
        builder.birthYear(birthYear);
    }

    @Nullable
    private static Integer determineBirthYear(@Nullable final String birthYear1, @Nullable final String birthYear2,
            @Nullable final LocalDate birthYear3) {
        if (birthYear1 != null) {
            return Integer.parseInt(birthYear1);
        }
        if (birthYear2 != null) {
            return Integer.parseInt(birthYear2);
        }
        if (birthYear3 != null) {
            return birthYear3.getYear();
        }
        return null;
    }

    private void setInformedConsent(@NotNull final ImmutableBaselineData.Builder builder, @NotNull final EcrfStudyEvent studyEvent) {
        boolean hasInformedConsentForm = false;
        for (final EcrfForm informedConsentForm : studyEvent.nonEmptyFormsPerOID(FORM_INFORMED_CONSENT)) {
            for (final EcrfItemGroup informedConsentItemGroup : informedConsentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_INFORMED_CONSENT)) {
                hasInformedConsentForm = true;
                builder.informedConsentDate(informedConsentItemGroup.readItemDate(FIELD_INFORMED_CONSENT_DATE, 0, DATE_FORMATTER));
                builder.informedConsentStatus(informedConsentForm.status());
                builder.informedConsentLocked(informedConsentForm.locked());
            }
        }

        if (!hasInformedConsentForm) {
            builder.informedConsentStatus(FormStatusState.UNKNOWN);
            builder.informedConsentLocked(false);
        }
    }

    private static void setDeathData(@NotNull final ImmutableBaselineData.Builder builder, @NotNull final EcrfPatient patient) {
        boolean hasDeathDateForm = false;
        for (final EcrfStudyEvent endStudyEvent : patient.studyEventsPerOID(STUDY_ENDSTUDY)) {
            for (final EcrfForm deathForm : endStudyEvent.nonEmptyFormsPerOID(FORM_DEATH)) {
                for (final EcrfItemGroup deathItemGroup : deathForm.nonEmptyItemGroupsPerOID(ITEMGROUP_DEATH)) {
                    hasDeathDateForm = true;
                    builder.deathDate(deathItemGroup.readItemDate(FIELD_DEATH_DATE, 0, DATE_FORMATTER));
                    builder.deathStatus(deathForm.status());
                    builder.deathLocked(deathForm.locked());
                }
            }
        }

        if (!hasDeathDateForm) {
            builder.deathStatus(FormStatusState.UNKNOWN);
            builder.deathLocked(false);
        }
    }
}
