package com.hartwig.hmftools.datamodel.purple;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.Set;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleFit {

    @NotNull
    public abstract Set<PurpleQCStatus> qcStatus();

    public abstract boolean hasSufficientQuality();

    public abstract boolean containsTumorCells();

    public abstract double purity();

    public abstract double ploidy();
}
