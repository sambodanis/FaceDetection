import java.util.ArrayList;
import java.util.Collections;

/**
 * Computes the singular value decomposition (SVD) of a matrix.
 *
 */
public class SVD {

    public SVD(Matrix A) {
        U = calcU(A);
        S = calcS(A);
        V = calcV(A);
    }

    public Matrix U() {
        return this.U;
    }

    public Matrix S() {
        return this.S;
    }

    public Matrix V() {
        return this.V;
    }

    private Matrix calcV(Matrix A) {
        Matrix ATA = A.transpose().matrixMultiplication(A);
        return ATA.eigenvectors().reverseColumns();
    }

    private Matrix calcS(Matrix A) {
        Matrix AEigenvalues = A.matrixMultiplication(A.transpose()).eigenvalues();
        ArrayList<Double> positiveEigenvalues = new ArrayList<Double>();
        for (int i = 0; i < AEigenvalues.getRows(); i++) {
            if (AEigenvalues.objectAtPoint(i, 0) > 0) positiveEigenvalues.add(AEigenvalues.objectAtPoint(i, 0));
        }

        Matrix S = new Matrix(new double[positiveEigenvalues.size() - 1]
                                        [positiveEigenvalues.size() - 1]);
        Collections.sort(positiveEigenvalues);
        Collections.reverse(positiveEigenvalues);

        int eigIndex = 0;
        for (int i = 0; i < S.getRows(); i++) {
            for (int j = 0; j < S.getColumns(); j++) {
                if (i == j) {
                    S.setObjectAtPoint(i, j, Math.sqrt(positiveEigenvalues.get(eigIndex++)));
                }
            }
        }
        return S;
    }

    private Matrix calcU(Matrix A) {
        Matrix AAT = A.matrixMultiplication(A.transpose());
        return AAT.eigenvectors().reverseColumns();
    }

    private Matrix U;
    private Matrix S;
    private Matrix V;
}
