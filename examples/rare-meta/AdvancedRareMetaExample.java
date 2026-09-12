import java.util.*;
import org.jlinalg.settest.*;
import org.jlinalg.meta.*;
import org.jlinalg.compute.BackendPolicy;

/** Synthetic aligned Gaussian summaries; all covariance values are on the U scale. */
public class AdvancedRareMetaExample {
    public static void main(String[] args) {
        var a=new SetTestScoreState(new double[]{2,-3,4},new double[]{4,1,1,1,9,2,1,2,4},3);
        var b=new SetTestScoreState(new double[]{4,1,2},new double[]{4,1,1,1,9,2,1,2,4},3);
        // Condition each cohort on variant 3 before pooling variants 1 and 2.
        var ca=SummaryScoreModels.condition(a,new int[]{0,1},new int[]{2});
        var cb=SummaryScoreModels.condition(b,new int[]{0,1},new int[]{2});
        var pooled=ScoreMetaAnalysis.pool(new double[][]{ca.scores(),cb.scores()},
            new double[][]{ca.information(),cb.information()},1).state();
        double[] weights={1,1};
        System.out.println("conditional burden: "+SummarySetTests.burden("gene",pooled,weights));
        var defaults=SetTestOptions.defaults();
        var deterministic=new SetTestOptions(defaults.variantFilter(),defaults.missingPolicy(),
            defaults.skatORhoGrid(),0,0,SkatOCalibration.DETERMINISTIC);
        System.out.println("deterministic SKAT-O: "+SummarySetTests.skatO("gene",pooled,weights,deterministic));
        var vt=SummaryScoreModels.variableThreshold("gene",pooled,weights,new double[]{.01,.02},100000,1234);
        System.out.println("VT selected burden (descriptive): "+vt.selectedBurden()+"; adjusted p="+vt.adjustedPValue());
        var hetero=SummaryScoreModels.heterogeneous(List.of(ca,cb));
        System.out.println("heterogeneous SKAT: "+SummarySetTests.skat("gene",hetero,new double[]{1,1,1,1}));
        // Comparable burden definitions in complete cohorts support scalar meta-analysis.
        var ba=SummarySetTests.burden("gene",ca,weights);var bb=SummarySetTests.burden("gene",cb,weights);
        var effects=List.of(new MetaStudy("A",ba.beta(),ba.standardError()),new MetaStudy("B",bb.beta(),bb.standardError()));
        System.out.println("random burden: "+MetaAnalysis.fit(effects,MetaAnalysisOptions.randomEffects(),BackendPolicy.CPU));
        // Influence diagnostics: fixed weights, retained null projections, recalculated tests.
        var omitted=SummaryScoreModels.subset(pooled,new int[]{1});
        System.out.println("omit variant 1: "+SummarySetTests.skatO("gene",omitted,new double[]{weights[1]},deterministic));
        System.out.println("omit cohort B: "+SummarySetTests.skatO("gene",ca,weights,deterministic));
    }
}
