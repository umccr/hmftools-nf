package com.hartwig.hmftools.protect.serve;

import com.hartwig.hmftools.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.sources.Source;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ServeEvidenceItem {

    @NotNull
    public abstract String genomicEvent();

    @NotNull
    public abstract Source source();

    public abstract boolean reported();

    @NotNull
    public abstract String treatment();

    public abstract boolean onLabel();

    @NotNull
    public abstract EvidenceLevel level();

    @NotNull
    public abstract EvidenceDirection direction();

    @NotNull
    public abstract String url();
}
