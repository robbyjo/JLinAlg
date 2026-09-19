/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class ProjectConfigurationFlagsTest {
    @TempDir Path dir;
    @Test void commandSpecificSwitchesEnableAndDisableWithoutBooleanArguments() throws Exception {
        Map<String,List<String>> switches=Map.of(
            "association",List.of("conditional-gwas-summary","dry-run","explain"),
            "meta-regression",List.of("no-intercept"),
            "iv-regression",List.of("no-intercept"),
            "arima-regression",List.of("no-intercept","smooth","parameter-uncertainty"),
            "penalized-regression",List.of("no-intercept","no-standardize"),
            "confounders",List.of("write-adjusted","center","no-center","scale","add-mean","nonparametric","mean-only"),
            "rare-meta",List.of("leave-variant-out","leave-cohort-out","cohort-results"),
            "coloc",List.of("no-trim"),"mr-estimate",List.of("joint"));
        for(var entry:switches.entrySet())for(String key:entry.getValue())for(boolean enabled:new boolean[]{true,false}) {
            config(entry.getKey(),key,Boolean.toString(enabled));
            String[] args=entry.getKey().equals("association")?new String[0]:new String[]{entry.getKey()};
            var expected=new ArrayList<>(List.of(args));if(enabled)expected.add("--"+key);
            assertEquals(expected,List.of(ProjectConfiguration.resolve(args,dir,dir.resolve("home")).arguments()));
        }
    }
    @Test void valuedBooleansAndSwitchPrecedenceRemainDistinct() throws Exception {
        for(String command:List.of("twas","pwas"))for(boolean value:new boolean[]{true,false}) {
            config(command,"joint",Boolean.toString(value));
            assertArrayEquals(new String[]{command,"--joint",String.valueOf(value)},resolve(command));
        }
        config("meta-regression","no-intercept","false");
        Path local=dir.resolve("home/.jlinalg/config.yaml");Files.createDirectories(local.getParent());
        Files.writeString(local,"schema_version: 1\ncommands:\n  meta-regression:\n    no-intercept: true\n");
        assertArrayEquals(new String[]{"meta-regression"},resolve("meta-regression"));
        assertArrayEquals(new String[]{"meta-regression","--no-intercept"},resolve("meta-regression","--no-intercept"));
        config("meta-regression","no-intercept","\"false\"");
        assertThrows(IllegalArgumentException.class,()->resolve("meta-regression"));
    }
    @Test void configuredMetaRegressionRunsAndChangesTheDesign() throws Exception {
        Path mods=dir.resolve("mods.tsv");Files.writeString(mods,"cohort\tdose\na\t1\nb\t2\nc\t3\n");
        List<String> cohorts=new ArrayList<>();
        for(int i=0;i<3;i++) {
            Path file=dir.resolve("c"+i+".tsv");Files.writeString(file,"feature_id\tbeta\tse\ngene\t"+(2+i)+"\t0.5\n");
            cohorts.addAll(List.of("--cohort",(char)('a'+i)+"="+file));
        }
        for(boolean enabled:new boolean[]{true,false}) {
            config("meta-regression","no-intercept",Boolean.toString(enabled));
            Path out=dir.resolve("out-"+enabled+".tsv");
            List<String> args=new ArrayList<>(List.of("meta-regression","--config",dir.resolve("jlinalg.yaml").toString(),
                "--local-config",dir.resolve("empty.yaml").toString(),"--model","fixed","--moderator-file",mods.toString(),"--moderators","dose","--out",out.toString()));
            Files.writeString(dir.resolve("empty.yaml"),"schema_version: 1\n");args.addAll(cohorts);
            var messages=new ByteArrayOutputStream();var stream=new PrintStream(messages);
            assertEquals(0,JLinAlgCli.run(args.toArray(String[]::new),stream,stream),messages.toString());
            var table=DelimitedData.read(out);assertEquals(enabled?1:2,table.rows().size());
            String[] dose=table.rows().stream().filter(r->r[1].equals("dose")).findFirst().orElseThrow();
            assertEquals(enabled?20./14:1,Double.parseDouble(dose[table.header().indexOf("beta")]),1e-12);
        }
    }
    private void config(String command,String key,String value) throws IOException {
        Files.writeString(dir.resolve("jlinalg.yaml"),"schema_version: 1\ncommands:\n  "+command+":\n    "+key+": "+value+"\n");
    }
    private String[] resolve(String... args) throws IOException {
        return ProjectConfiguration.resolve(args,dir,dir.resolve("home")).arguments();
    }
}
