package com.hartwig.hmftools.serve.actionability.util;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.sources.ImmutableSources;
import com.hartwig.hmftools.serve.sources.Sources;

import org.jetbrains.annotations.NotNull;

public final class ActionableFileFunctions {

    public static final String FIELD_DELIMITER = "\t";

    private static final String URL_DELIMITER = ",";

    private ActionableFileFunctions() {
    }

    @NotNull
    public static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("sourceEvent")
                .add("source")
                .add("treatment")
                .add("cancerType")
                .add("doid")
                .add("tumorLocationBlacklisting")
                .add("level")
                .add("direction")
                .add("sourceUrls")
                .add("evidenceUrls")
                .toString();
    }

    @NotNull
    public static ActionableEvent fromLine(@NotNull String[] values, int startingPosition) {

        return new ActionableEvent() {

            @NotNull
            @Override
            public Sources source() {
                return ImmutableSources.builder()
                        .sourceEvent(values[startingPosition])
                        .source(Knowledgebase.valueOf(values[startingPosition + 1]))
                        .build();
            }

            @NotNull
            @Override
            public String treatment() {
                return values[startingPosition + 2];
            }

            @NotNull
            @Override
            public String cancerType() {
                return values[startingPosition + 3];
            }

            @NotNull
            @Override
            public String doid() {
                return values[startingPosition + 4];
            }

            @NotNull
            @Override
            public String tumorLocationBlacklisting() {
                return values[startingPosition + 5];
            }

            @NotNull
            @Override
            public EvidenceLevel level() {
                return EvidenceLevel.valueOf(values[startingPosition + 6]);
            }

            @NotNull
            @Override
            public EvidenceDirection direction() {
                return EvidenceDirection.valueOf(values[startingPosition + 7]);
            }

            @NotNull
            @Override
            public Set<String> sourceUrls() {
                int urlPosition = startingPosition + 8;
                return values.length > urlPosition ? stringToUrls(values[urlPosition]) : Sets.newHashSet();
            }

            @NotNull
            @Override
            public Set<String> evidenceUrls() {
                int urlPosition = startingPosition + 9;
                return values.length > urlPosition ? stringToUrls(values[urlPosition]) : Sets.newHashSet();
            }
        };
    }

    @NotNull
    public static String toLine(@NotNull ActionableEvent event) {
        return new StringJoiner(FIELD_DELIMITER).add(event.source().sourceEvent())
                .add(event.source().source().toString())
                .add(event.treatment())
                .add(event.cancerType())
                .add(event.doid())
                .add(event.tumorLocationBlacklisting())
                .add(event.level().toString())
                .add(event.direction().toString())
                .add(urlsToString(event.sourceUrls()))
                .add(urlsToString(event.evidenceUrls()))
                .toString();
    }

    @NotNull
    private static Set<String> stringToUrls(@NotNull String fieldValue) {
        return Sets.newHashSet(fieldValue.split(URL_DELIMITER));
    }

    @NotNull
    private static String urlsToString(@NotNull Set<String> urls) {
        StringJoiner joiner = new StringJoiner(URL_DELIMITER);
        for (String url : urls) {
            joiner.add(url);
        }
        return joiner.toString();
    }
}