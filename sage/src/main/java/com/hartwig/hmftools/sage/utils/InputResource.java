package com.hartwig.hmftools.sage.utils;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;

import java.io.File;
import java.net.URI;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

public class InputResource {
    private final String locator;

    public InputResource(String locator) {
        this.locator = locator;
    }
//    public SamReader reader() {
//
//    }
    public boolean isValid() {
        // Determines if we are dealing with htsget/http/etc... resources instead of local files first
        URI uri = URI.create(this.locator);

        if (uri.getScheme() != null) {
            SamInputResource res = SamInputResource.of(uri, null);
            SG_LOGGER.info("Remote URI supplied for tumor bam({})", res.toString());
            return true;
        } else if (!new File(this.locator).exists()) {
            SG_LOGGER.error("Unable to locate tumor bam({})", this.locator);
            return false;
        } else {
            // All checks succeeded
            return true;
        }
    }
}
