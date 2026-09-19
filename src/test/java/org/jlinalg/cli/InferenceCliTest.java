/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class InferenceCliTest {
    @TempDir Path directory;
    private Path fixture(String name) throws IOException {
        Path file=directory.resolve(name+".tsv");
        if(!Files.exists(file))try(var in=getClass().getResourceAsStream("/r-reference/regression-inference/"+name+".tsv")){Files.copy(in,file);}
        return file;
    }
    private int run(String... args) {
        var errors=new ByteArrayOutputStream();
        int code=JLinAlgCli.run(args,new PrintStream(new ByteArrayOutputStream()),new PrintStream(errors));
        if(code!=0)lastError=errors.toString();
        return code;
    }
    private String lastError;
    @Test void errorsPreserveExistingFilesAndDoNotPublishPartialResults() throws IOException {
        Path input=fixture("survey"),out=directory.resolve("protected.tsv");Files.writeString(out,"keep\n");
        assertEquals(2,run("rare-events-logit","--input",input.toString(),"--response","rare","--out",out.toString(),"--no-log"));
        assertEquals("keep\n",Files.readString(out));assertFalse(Files.exists(Path.of(out+".metadata.tsv")));
        Path bad=directory.resolve("invalid.tsv");
        assertEquals(2,run("rare-events-logit","--input",input.toString(),"--response","rare","--prevalence","0","--out",bad.toString(),"--no-log"));
        assertFalse(Files.exists(bad));
        assertEquals(2,run("survey-regression","--input",input.toString(),"--response","y","--out",bad.toString(),"--no-log"));
        assertFalse(Files.exists(bad));
        for(String command:List.of("censored-regression","ordinal-regression","rare-events-logit","survey-regression"))assertEquals(0,run(command,"--help"));
    }
    @Test void likelihoodAndSamplingCommandsProducePredictionsAndDesignMetadata() throws IOException {
        Path censored=fixture("censored"),ordinal=fixture("ordinal"),survey=fixture("survey");
        for(String distribution:List.of("gaussian","weibull","lognormal","exponential")) {
            boolean gaussian=distribution.equals("gaussian");Path out=directory.resolve(distribution+".tsv");
            var args=new ArrayList<>(List.of("censored-regression","--input",censored.toString(),"--response",gaussian?"tobit":"time",
                "--censor",gaussian?"cens":"censTime","--predictors","c,z1","--distribution",distribution,"--predict",censored.toString(),"--out",out.toString()));
            args.addAll(gaussian?List.of("--lower","-0.4","--upper","2.5"):List.of("--time","2"));
            assertEquals(0,run(args.toArray(String[]::new)),lastError);
            assertEquals(241,Files.readAllLines(Path.of(out+".predictions.tsv")).size());
        }
        for(String link:List.of("logit","probit")) {
            Path out=directory.resolve("ordered-"+link+".tsv");
            assertEquals(0,run("ordinal-regression","--input",ordinal.toString(),"--response","ordinal","--predictors","c,z1",
                "--categories","4","--link",link,"--predict",ordinal.toString(),"--out",out.toString()),lastError);
            assertTrue(Files.readString(Path.of(out+".predictions.tsv")).contains("category_3"));
        }
        Path rare=directory.resolve("rare.tsv"),weighted=directory.resolve("weighted.tsv");
        assertEquals(0,run("rare-events-logit","--input",survey.toString(),"--response","rare","--predictors","c,z1","--prevalence","0.04","--out",rare.toString()),lastError);
        assertTrue(Files.readString(Path.of(rare+".metadata.tsv")).contains("prior_intercept_shift"));
        assertEquals(0,run("survey-regression","--input",survey.toString(),"--response","y","--predictors","c,z1","--weights","weight",
            "--strata","stratum","--psu","cluster","--out",weighted.toString()),lastError);
        assertTrue(Files.readString(Path.of(weighted+".metadata.tsv")).contains("degrees_of_freedom\t55"));
    }
}
