import java.util.ArrayList;

/**
 * Computes the QR decomposition of a matrix.
 * Uses the Gram Schmidt Process http://en.wikipedia.org/wiki/QR_decomposition
 */
public class QRDecomposition {

    public QRDecomposition(Matrix A) {
        eVectors = new ArrayList<Matrix>();
        aVectors = new ArrayList<Matrix>();
        uVectors = new ArrayList<Matrix>();
        for (int i = 0; i < A.getColumns(); i++) {
            aVectors.add(A.extractColumn(i));
            if (i == 0) {
                uVectors.add(aVectors.get(i));
            } else {
                Matrix projectionSum = new Matrix(new
                        double[aVectors.get(aVectors.size() - 1).getRows()]
                              [aVectors.get(aVectors.size() - 1).getColumns()]);
                for (int j = 0; j < i; j++) {
                    projectionSum = projectionSum.matrixAddition(
                            Matrix.projection(eVectors.get(j), aVectors.get(aVectors.size() - 1)));
//                                  .scalarMultiplication(-1));
                }
                uVectors.add(aVectors.get(i)
                        .matrixAddition(projectionSum.scalarMultiplication(-1)));
            }
            eVectors.add(Matrix.normalise(uVectors.get(uVectors.size() - 1)));
        }
    }

    public Matrix Q() {
        Matrix result = eVectors.get(0);
        if (eVectors.size() == 1) return result;
        for (int i = 1; i < result.getRows(); i++) {
            result = result.stickToMatrix(eVectors.get(i));
        }
        return result.scalarMultiplication(-1);
    }

    public Matrix R() {
        Matrix result = new Matrix(new double[eVectors.size()][aVectors.size()]);
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                if (i <= j) {
                    double innerProduct = eVectors.get(i).transpose()
                            .matrixMultiplication(aVectors.get(j))
                            .objectAtPoint(0,0);
                    result.setObjectAtPoint(i, j, innerProduct);
                }
            }
        }
        return result.scalarMultiplication(-1);
    }

    ArrayList<Matrix> eVectors;
    ArrayList<Matrix> aVectors;
    ArrayList<Matrix> uVectors;
}
