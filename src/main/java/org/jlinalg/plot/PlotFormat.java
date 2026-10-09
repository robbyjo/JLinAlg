/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.nio.file.Path;
import java.util.Locale;

/** File format used when exporting a plot. */
public enum PlotFormat {
    PNG("png"), JPEG("jpg"), BMP("bmp"), GIF("gif"), TIFF("tiff"),
    SVG("svg"), EPS("eps"), PDF("pdf");

    private final String extension;

    PlotFormat(String extension) {
        this.extension = extension;
    }

    /** Conventional filename extension, without a leading dot. */
    public String extension() {
        return extension;
    }

    /** Infers a format from a destination filename. */
    public static PlotFormat fromPath(Path path) {
        if (path == null) throw new IllegalArgumentException("plot path is required");
        String name = path.getFileName().toString();
        int dot = name.lastIndexOf('.');
        if (dot < 0 || dot == name.length() - 1)
            throw new IllegalArgumentException("plot filename needs an extension: " + name);
        String value = name.substring(dot + 1).toLowerCase(Locale.ROOT);
        if (value.equals("jpeg")) return JPEG;
        if (value.equals("tif")) return TIFF;
        for (PlotFormat format : values())
            if (format.extension.equals(value)) return format;
        throw new IllegalArgumentException("unsupported plot format: " + value);
    }
}
