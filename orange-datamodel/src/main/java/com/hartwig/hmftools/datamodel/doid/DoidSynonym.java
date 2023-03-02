package com.hartwig.hmftools.datamodel.doid;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Value.Immutable
@Value.Style(allParameters = true,
        passAnnotations = {NotNull.class, Nullable.class})
public abstract class DoidSynonym {

    @NotNull
    public abstract String pred();

    @NotNull
    public abstract String val();

    @NotNull
    public abstract List<String> xrefs();
}

