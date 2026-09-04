/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

import org.junit.jupiter.api.Test;

class LdClumpDispatchTest {
    @Test
    void topLevelCliDispatchesClumpHelp() {
        ByteArrayOutputStream standardOutput = new ByteArrayOutputStream();
        ByteArrayOutputStream errorOutput = new ByteArrayOutputStream();

        int status = JLinAlgCli.run(new String[] {"clump", "--help"},
            new PrintStream(standardOutput, true, StandardCharsets.UTF_8),
            new PrintStream(errorOutput, true, StandardCharsets.UTF_8));

        assertEquals(0, status);
        assertTrue(standardOutput.toString(StandardCharsets.UTF_8)
            .contains("--ld-threshold"));
        assertEquals("", errorOutput.toString(StandardCharsets.UTF_8));
    }
}
