/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

/** One ordered repeated-fit failure. */
public record AssociationFailure(int index, String name,
                                 String exceptionType, String message) { }
