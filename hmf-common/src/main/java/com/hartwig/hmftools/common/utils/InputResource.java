package com.hartwig.hmftools.common.utils;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Paths;

public class InputResource {
    private String locator;
    private static SamReader reader = null;
    private static final Logger LOGGER = LogManager.getLogger(InputResource.class);

    public InputResource(String locator) {
        this.locator = locator;
    }

    /***
     * Returns a SamReader irrespective of file residing in local filesystem or
     * whether it's a remote resource.
     *
     * @param genome
     * @param locator
     * @return SamReader reader
     */
    public SamReader reader(RefGenomeSource genome, String locator) {
        this.locator = locator;

        try {
            URI uri = URI.create(locator);

            if (uri.getScheme() != null) {
                if (!Files.exists(Paths.get(locator)))
                    this.reader = SamReaderFactory.makeDefault()
                            .referenceSource(new ReferenceSource(genome.refGenomeFile()))
                            .open(new File(locator));
            } else {
                this.reader = SamReaderFactory.makeDefault()
                        .referenceSource(new ReferenceSource(genome.refGenomeFile()))
                        .open(SamInputResource.of(uri, null));
            }
        } catch (Exception e) {
            LOGGER.error(e.toString());
        }

        return reader;
    }

    /***
     * Determines if the object is valid iif:
     *
     *  1. Exists and is readable in local filesystem.
     *  2. Exists and is readable in remote URI?
     *
     * @return true if valid
     */
    public boolean isValid() {
        // Determines if we are dealing with htsget/http/etc... resources instead of local files first
        URI uri = URI.create(this.locator);

        if (uri.getScheme() != null) {
            SamInputResource res = SamInputResource.of(uri, null);
            LOGGER.info("Remote URI supplied for tumor bam({})", res.toString());
            return true;
        } else if (!new File(this.locator).exists()) {
            LOGGER.error("Unable to locate tumor bam({})", this.locator);
            return false;
        } else {
            // All checks succeeded
            return true;
        }
    }
}
