import org.jlinalg.settest.ScoreMetaAnalysis;
import org.jlinalg.settest.SetTestOptions;
import org.jlinalg.settest.SummarySetTests;

/** Aligned, synthetic quantitative-trait summaries from independent cohorts. */
public class RareMetaExample {
    public static void main(String[] args) {
        // Cohort by variant; covariance blocks are row-major on the information scale.
        double[][] scores = {{2, -3}, {4, 1}};
        double[][] covariance = {{4, 1, 1, 9}, {16, 2, 2, 4}};
        var pooled = ScoreMetaAnalysis.pool(scores, covariance, 1);
        if (pooled == null) throw new IllegalStateException("No eligible variants");
        var state = pooled.state();
        double[] u = state.scores(), v = state.information();
        for (int j = 0; j < state.variants(); j++) {
            var single = SummarySetTests.singleVariant("v" + (pooled.indices()[j] + 1),
                u[j], v[j * state.variants() + j]);
            System.out.println("single\tv" + (pooled.indices()[j] + 1)
                + "\t" + pooled.directions()[j] + "\t" + single.pValue());
        }
        double[] equal = {1, 1}, weighted = {1, 2};
        var burden = SummarySetTests.burden("GENE1", state, equal);
        var weightedBurden = SummarySetTests.burden("GENE1", state, weighted);
        var skat = SummarySetTests.skat("GENE1", state, weighted);
        var skato = SummarySetTests.skatO("GENE1", state, weighted, SetTestOptions.defaults());
        System.out.println("equal_burden\tbeta=" + burden.beta() + "\tse=" + burden.standardError()
            + "\tp=" + burden.pValue());
        System.out.println("weighted_burden\tbeta=" + weightedBurden.beta()
            + "\tse=" + weightedBurden.standardError() + "\tp=" + weightedBurden.pValue());
        System.out.println("skat\tQ=" + skat.statistic() + "\tp=" + skat.pValue());
        System.out.println("skat-o\tadjusted_p=" + skato.adjustedPValue()
            + "\tsimulations=" + skato.simulations() + "\tseed=" + skato.randomSeed());
    }
}
