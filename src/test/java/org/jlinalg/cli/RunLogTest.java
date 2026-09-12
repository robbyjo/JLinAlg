/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.file.*;
import java.time.Instant;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class RunLogTest {
    @TempDir Path directory;

    @Test void readableRuntimeIncludesDaysAndSubsecondPrecision() {
        assertEquals("0.000 s",RunTiming.readable(0));
        assertEquals("0.007 s",RunTiming.readable(7));
        assertEquals("59.999 s",RunTiming.readable(59999));
        assertEquals("1 min 0.000 s",RunTiming.readable(60000));
        assertEquals("1 h 0 min 0.001 s",RunTiming.readable(3600001));
        assertEquals("2 d 3 h 4 min 5.006 s",RunTiming.readable(183845006));
        assertThrows(IllegalArgumentException.class,()->RunTiming.readable(-1));
    }

    @Test void startIsFlushedBeforeWorkAndCompletionIsWrittenOnlyOnce() throws Exception {
        for(boolean plain:new boolean[]{false,true}) {
            Path path=directory.resolve(plain+".log");
            try(RunLog log=plain?RunLog.openPlain(path,false):RunLog.open(path,false)) {
                String live=Files.readString(path);
                assertTrue(live.contains("started="));assertFalse(live.contains("finished="));
                log.info("work=done");log.complete("complete");log.complete("failed");
            }
            assertTiming(path,"complete");
            assertEquals(1,Files.readString(path).split("finished=",-1).length-1);
        }
    }

    @Test void failedAndAppendedRunsHaveTheirOwnCompleteTimingBlocks() throws Exception {
        Path path=directory.resolve("append.log");
        try(RunLog log=RunLog.openPlain(path,false)) { log.info("work=failed"); }
        assertTiming(path,"failed");
        String before=Files.readString(path);
        try(RunLog log=RunLog.openPlain(path,true)) { log.complete("complete"); }
        String after=Files.readString(path);
        assertTrue(after.startsWith(before));
        assertEquals(2,after.split("started=",-1).length-1);
        assertEquals(2,after.split("finished=",-1).length-1);
        assertEquals(2,after.split("elapsed=",-1).length-1);
    }

    @Test void previouslyUnloggedCommandsReceiveALogWithoutChangingStdout() throws Exception {
        Path log=directory.resolve("catalog.log");
        var stdout=new ByteArrayOutputStream();var stderr=new ByteArrayOutputStream();
        assertEquals(0,JLinAlgCli.run(new String[]{"ld-db","list","--log",log.toString()},
            new PrintStream(stdout),new PrintStream(stderr)));
        assertTrue(stdout.toString().contains("1000g"));
        assertFalse(stdout.toString().contains("started="));
        assertTiming(log,"complete");
        assertEquals(0,JLinAlgCli.run(new String[]{"ld-db","list","--no-log"},
            new PrintStream(stdout),new PrintStream(stderr)));
        Path failure=directory.resolve("failed-command.log");
        assertEquals(2,JLinAlgCli.run(new String[]{"ld-db","unknown","--log",failure.toString()},
            new PrintStream(stdout),new PrintStream(stderr)));
        assertTiming(failure,"failed");
    }

    @Test void explicitLogCannotOverwriteAnInput() throws Exception {
        Path input=directory.resolve("input.tsv");Files.writeString(input,"keep this input");
        var error=new ByteArrayOutputStream();
        assertEquals(2,JLinAlgCli.run(new String[]{"mr-estimate","--input",input.toString(),
            "--log",input.toString()},new PrintStream(new ByteArrayOutputStream()),new PrintStream(error)));
        assertEquals("keep this input",Files.readString(input));
    }

    static void assertTiming(Path path,String status) throws IOException {
        String text=Files.readString(path);
        String started=field(text,"started"),finished=field(text,"finished");
        assertTrue(started.endsWith("Z"));assertTrue(finished.endsWith("Z"));
        assertFalse(Instant.parse(finished).isBefore(Instant.parse(started)));
        long elapsed=Long.parseLong(field(text,"elapsed_ms"));
        assertTrue(elapsed>=0);assertEquals(RunTiming.readable(elapsed),field(text,"elapsed"));
        assertEquals(status,field(text,"status"));
    }
    private static String field(String text,String key) {
        for(String line:text.split("\\R")) {
            int position=line.indexOf(key+"=");
            if(position==0 || position>0 && line.substring(0,position).endsWith("] "))
                return line.substring(position+key.length()+1);
        }
        throw new AssertionError("missing log field: "+key+" in "+text);
    }
}
