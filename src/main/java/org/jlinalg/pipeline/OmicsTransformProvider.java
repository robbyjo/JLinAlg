/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.util.Map;

/**
 * Service-provider interface for trusted custom row-wise transforms.
 *
 * <p>Providers are registered through
 * {@code META-INF/services/org.jlinalg.pipeline.OmicsTransformProvider}.</p>
 */
public interface OmicsTransformProvider {
    String name();
    OmicsTransform create(Map<String, String> parameters);
}
