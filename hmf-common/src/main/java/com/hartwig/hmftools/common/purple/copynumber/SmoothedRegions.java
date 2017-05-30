package com.hartwig.hmftools.common.purple.copynumber;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

class SmoothedRegions {

    private final List<PurpleCopyNumber> smoothedRegions = Lists.newArrayList();

    private final List<PurpleCopyNumber> broadRegions;
    private final List<FittedRegion> copyNumbers;

    private static final double DIPLOID_MIN_RATIO = 0.75;
    private static final double DIPLOID_MAX_RATIO = 1.25;

    private static final double MIN_COPY_NUMBER_RANGE = 0.3;
    private static final double MAX_COPY_NUMBER_RANGE = 1.3;

    SmoothedRegions(@NotNull final List<PurpleCopyNumber> broadRegions,
            @NotNull final List<FittedRegion> fittedRegions) {
        this.broadRegions = broadRegions;
        this.copyNumbers = fittedRegions;

        run();
    }

    @NotNull
    List<PurpleCopyNumber> getSmoothedRegions() {
        return smoothedRegions;
    }

    private void run() {
        if (!broadRegions.isEmpty()) {
            int largestIncludedIndex = -1;
            PurpleCopyNumberBuilder currentBuilder;

            for (int i = 0; i < broadRegions.size(); i++) {
                final PurpleCopyNumber currentRegion = broadRegions.get(i);
                int startOfRegionIndex = indexOfStart(largestIncludedIndex + 1, currentRegion);
                int endOfRegionIndex = indexOfEnd(startOfRegionIndex, currentRegion);

                // JOBA: Start new builder
                currentBuilder = new PurpleCopyNumberBuilder(copyNumbers.get(startOfRegionIndex));

                // JOBA: Go backwards to previous end
                currentBuilder = backwards(startOfRegionIndex - 1, largestIncludedIndex + 1, currentBuilder);

                // JOBA: Go forwards to end of region
                currentBuilder = forwards(startOfRegionIndex + 1, endOfRegionIndex, currentBuilder);

                boolean isLastBroadRegion = i == broadRegions.size() - 1;
                if (isLastBroadRegion) {
                    currentBuilder = forwards(endOfRegionIndex + 1, copyNumbers.size() - 1, currentBuilder);
                    smoothedRegions.add(currentBuilder.build());
                } else {
                    int nextStartIndex = indexOfStart(endOfRegionIndex + 1, broadRegions.get(i + 1));
                    largestIncludedIndex = forwardsUntilDifferent(endOfRegionIndex + 1, nextStartIndex - 1,
                            currentBuilder);
                    smoothedRegions.add(currentBuilder.build());
                }
            }
        }
    }

    private int forwardsUntilDifferent(int startIndex, int endIndex, @NotNull PurpleCopyNumberBuilder builder) {
        for (int i = startIndex; i <= endIndex; i++) {
            FittedRegion copyNumber = copyNumbers.get(i);
            if (isSimilar(copyNumber, builder)) {
                builder.extendRegion(copyNumber);
            } else {
                return i - 1;
            }
        }

        return endIndex;
    }

    @NotNull
    private PurpleCopyNumberBuilder forwards(int startIndex, int endIndex,
            final @NotNull PurpleCopyNumberBuilder builder) {
        PurpleCopyNumberBuilder current = builder;
        for (int i = startIndex; i <= endIndex; i++) {
            FittedRegion copyNumber = copyNumbers.get(i);
            if (isSimilar(copyNumber, builder)) {
                builder.extendRegion(copyNumber);
            } else {
                smoothedRegions.add(builder.build());
                current = new PurpleCopyNumberBuilder(copyNumber);
            }
        }

        return current;
    }

    @NotNull
    private PurpleCopyNumberBuilder backwards(int startIndex, int endIndex,
            @NotNull final PurpleCopyNumberBuilder forwardBuilder) {
        final Deque<PurpleCopyNumber> preRegions = new ArrayDeque<>();
        PurpleCopyNumberBuilder reverseBuilder = forwardBuilder;

        for (int i = startIndex; i >= endIndex; i--) {
            final FittedRegion copyNumber = copyNumbers.get(i);
            if (isSimilar(copyNumber, reverseBuilder)) {
                reverseBuilder.extendRegion(copyNumber);
            } else {
                if (reverseBuilder != forwardBuilder) {
                    preRegions.addFirst(reverseBuilder.build());
                }
                reverseBuilder = new PurpleCopyNumberBuilder(copyNumber);
            }
        }

        if (reverseBuilder != forwardBuilder) {
            preRegions.addFirst(reverseBuilder.build());
        }

        smoothedRegions.addAll(preRegions);
        return forwardBuilder;
    }

    private static boolean isSimilar(@NotNull final FittedRegion copyNumber,
            @NotNull final PurpleCopyNumberBuilder builder) {
        int bafCount = copyNumber.bafCount();
        if (!isDiploid(copyNumber)) {
            return true;
        }

        double tumorCopyNumberDeviation = Math.abs(copyNumber.tumorCopyNumber() - builder.averageTumorCopyNumber());
        double refNormalisedCopyNumberDeviation = Math.abs(
                copyNumber.refNormalisedCopyNumber() - builder.averageRefNormalisedCopyNumber());
        double copyNumberDeviation = Math.min(tumorCopyNumberDeviation, refNormalisedCopyNumberDeviation);
        if (Doubles.greaterThan(copyNumberDeviation, allowedCopyNumberDeviation(bafCount))) {
            return false;
        }

        if (bafCount > 0 && !Doubles.isZero(builder.averageObservedBAF())) {
            double bafDeviation = Math.abs(copyNumber.observedBAF() - builder.averageObservedBAF());
            if (Doubles.greaterThan(bafDeviation, allowedBAFDeviation(bafCount))) {
                return false;
            }
        }

        return true;
    }

    static double allowedBAFDeviation(int bafCount) {
        return 1d / Math.pow(Math.max(1, bafCount), 0.95) * 0.38 + 0.022;
    }

    static double allowedCopyNumberDeviation(int bafCount) {
        if (bafCount >= 10) {
            return MIN_COPY_NUMBER_RANGE;
        }
        return (MIN_COPY_NUMBER_RANGE - MAX_COPY_NUMBER_RANGE) / 10 * bafCount + MAX_COPY_NUMBER_RANGE;
    }

    private static boolean isDiploid(@NotNull final FittedRegion copyNumber) {
        return Doubles.greaterOrEqual(copyNumber.observedNormalRatio(), DIPLOID_MIN_RATIO) && Doubles.lessOrEqual(
                copyNumber.observedNormalRatio(), DIPLOID_MAX_RATIO);
    }

    private int indexOfEnd(int minIndex, @NotNull PurpleCopyNumber region) {
        return indexOf(minIndex, copyNumber -> copyNumber.end() == region.end());
    }

    private int indexOfStart(int minIndex, @NotNull PurpleCopyNumber region) {
        return indexOf(minIndex, copyNumber -> copyNumber.start() == region.start());
    }

    private int indexOf(int minIndex, @NotNull Predicate<FittedRegion> predicate) {
        for (int i = minIndex; i < copyNumbers.size(); i++) {
            FittedRegion copyNumber = copyNumbers.get(i);
            if (predicate.test(copyNumber)) {
                return i;
            }
        }

        throw new IllegalArgumentException();
    }
}
