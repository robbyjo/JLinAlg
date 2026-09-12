import java.nio.file.Path;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.genetics.GenomicRelationshipMatrix;
import org.jlinalg.genetics.GenomicRelationshipOptions;
import org.jlinalg.pipeline.VariantSources;

/** Reproduce the tutorial's four-sample matrix with the streaming Java API. */
public class GrmExample {
    public static void main(String[] args) throws Exception {
        var source = VariantSources.open(Path.of("examples/grm/dosages.tsv"));
        var grm = GenomicRelationshipMatrix.fromSource(source,
            new GenomicRelationshipOptions(0, 0.5), BackendPolicy.CPU, 2);
        System.out.println("considered=" + grm.variantsConsidered()
            + " used=" + grm.variantsUsed() + " excluded=" + grm.variantsExcluded());
        System.out.println("sample\t" + String.join("\t", grm.sampleIds()));
        for (String row : grm.sampleIds()) {
            System.out.print(row);
            for (String column : grm.sampleIds())
                System.out.print("\t" + grm.relationship(row, column));
            System.out.println();
        }
        System.out.println("K(a,b)=" + grm.relationship("a", "b"));
        System.out.println("kinship(a,b)=" + grm.kinshipCoefficient("a", "b"));
    }
}
