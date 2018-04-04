package com.hartwig.hmftools.patientdb.matchers;

import static com.hartwig.hmftools.patientdb.readers.cpct.BiopsyTreatmentResponseReader.FORM_TUMOR_MEASUREMENT;

import java.time.LocalDate;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentResponseData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class TreatmentResponseMatcher {
    private TreatmentResponseMatcher() {
    }

    @NotNull
    public static MatchResult<BiopsyTreatmentResponseData> matchTreatmentResponsesToTreatments(@NotNull String patientId,
            @NotNull List<BiopsyTreatmentData> treatments, @NotNull List<BiopsyTreatmentResponseData> responses) {
        final List<BiopsyTreatmentResponseData> matchedResponses = Lists.newArrayList();
        final List<ValidationFinding> findings = Lists.newArrayList();
        Collections.sort(responses);

        List<BiopsyTreatmentData> sortedTreatments = sortAndFilter(treatments);
        if (hasOverlappingTreatments(sortedTreatments)) {
            if (!responses.isEmpty()) {
                findings.add(responseMatchFinding(patientId, "treatments are overlapping. Cannot match any response.",
                        "treatments: " + sortedTreatments));
            }
            return new MatchResult<>(responses, findings);
        }

        Iterator<BiopsyTreatmentData> treatmentIterator = sortedTreatments.iterator();
        BiopsyTreatmentData currentTreatment = treatmentIterator.hasNext() ? treatmentIterator.next() : null;
        LocalDate firstTreatmentStart = currentTreatment != null ? currentTreatment.startDate() : null;
        BiopsyTreatmentData nextTreatment = treatmentIterator.hasNext() ? treatmentIterator.next() : null;
        boolean hasNewBaselineResponseFound = false;
        for (BiopsyTreatmentResponseData response : responses) {
            LocalDate responseDate = response.date();
            if (responseMatchable(responseDate, firstTreatmentStart)) {
                LocalDate nextTreatmentStart = nextTreatment != null ? nextTreatment.startDate() : null;
                while (nextTreatmentStart != null && responseDate.isAfter(nextTreatmentStart)) {
                    currentTreatment = nextTreatment;
                    nextTreatment = treatmentIterator.hasNext() ? treatmentIterator.next() : null;
                    nextTreatmentStart = nextTreatment != null ? nextTreatment.startDate() : null;
                    hasNewBaselineResponseFound = false;
                }

                if (hasNewBaselineResponseFound) {
                    matchedResponses.add(response);
                    findings.add(responseMatchFinding(patientId,
                            "response after new baseline and before next treatment",
                            "response: " + response));
                } else {
                    String actualResponse = response.response();
                    if (actualResponse != null && actualResponse.equalsIgnoreCase("ne")) {
                        matchedResponses.add(response);
                        hasNewBaselineResponseFound = true;
                    } else {
                        matchedResponses.add(ImmutableBiopsyTreatmentResponseData.builder()
                                .from(response)
                                .treatmentId(currentTreatment.id())
                                .build());
                    }
                }
            } else {
                matchedResponses.add(response);
            }
        }

        return new MatchResult<>(matchedResponses, findings);
    }

    private static boolean hasOverlappingTreatments(@NotNull List<BiopsyTreatmentData> sortedTreatments) {
        Iterator<BiopsyTreatmentData> iterator = sortedTreatments.iterator();

        BiopsyTreatmentData current = iterator.hasNext() ? iterator.next() : null;
        while (iterator.hasNext()) {
            assert current != null;
            BiopsyTreatmentData next = iterator.next();
            LocalDate currentEndDate = current.endDate();
            LocalDate nextStartDate = next.startDate();
            assert nextStartDate != null;
            if (currentEndDate == null || currentEndDate.isAfter(nextStartDate)) {
                return true;
            }
            current = next;
        }

        return false;
    }

    @NotNull
    private static List<BiopsyTreatmentData> sortAndFilter(@NotNull List<BiopsyTreatmentData> treatments) {
        List<BiopsyTreatmentData> sortedTreatments = Lists.newArrayList(treatments);
        Collections.sort(sortedTreatments);
        sortedTreatments.removeIf(treatment -> treatment.startDate() == null);
        return sortedTreatments;
    }

    private static boolean responseMatchable(@Nullable LocalDate responseDate, @Nullable LocalDate firstTreatmentStart) {
        return responseDate != null && firstTreatmentStart != null && responseDate.isAfter(firstTreatmentStart);
    }

    @NotNull
    private static ValidationFinding responseMatchFinding(@NotNull String patientIdentifier, @NotNull String message,
            @NotNull String details) {
        return ValidationFinding.of("match", patientIdentifier, FORM_TUMOR_MEASUREMENT, message, FormStatus.unknown(), details);
    }
}
