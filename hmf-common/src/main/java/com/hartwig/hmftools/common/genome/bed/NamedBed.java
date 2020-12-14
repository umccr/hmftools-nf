package com.hartwig.hmftools.common.genome.bed;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface NamedBed extends GenomeRegion {

    @NotNull
    String name();
}
