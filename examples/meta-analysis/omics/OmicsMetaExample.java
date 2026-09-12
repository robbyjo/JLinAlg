import org.jlinalg.meta.MetaAnalysis;
import org.jlinalg.meta.MetaAnalysisOptions;

/** Synthetic expression effects: log2 fold change for the same case/control contrast. */
public class OmicsMetaExample {
    public static void main(String[] args) {
        String[] genes = {"gene_A", "gene_B", "gene_C"};
        // Feature-major rows, with cohorts A, B, C in each row.
        double[] beta = {.2, .5, .1, -.3, -.1, -.2, .4, Double.NaN, Double.NaN};
        double[] se = {.1, .2, .15, .2, .1, .15, .1, Double.NaN, Double.NaN};
        var batch = MetaAnalysis.prepareBatch(beta, se, genes.length, 3, 2);
        var fixed = batch.fit(MetaAnalysisOptions.fixedEffect(), 1);
        var random = batch.fit(MetaAnalysisOptions.randomEffects(), 1);
        double[] fixedBeta = fixed.pooledEffectSizes(), fixedSe = fixed.standardErrors();
        double[] randomBeta = random.pooledEffectSizes(), randomSe = random.standardErrors();
        int[] counts = batch.cohortCounts();
        System.out.println("feature_id\tn_cohorts\tdirection\tfixed_beta\tfixed_se\trandom_beta\trandom_se");
        for (int i = 0; i < genes.length; i++) {
            System.out.println(genes[i] + "\t" + counts[i] + "\t" + batch.direction(i)
                + "\t" + fixedBeta[i] + "\t" + fixedSe[i]
                + "\t" + randomBeta[i] + "\t" + randomSe[i]);
        }
    }
}
