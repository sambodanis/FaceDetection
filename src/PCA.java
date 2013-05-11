/**
 * Created with IntelliJ IDEA.
 * User: sambodanis
 * Date: 04/05/2013
 * Time: 16:02
 * To change this template use File | Settings | File Templates.
 */
public class PCA {

    public PCA(Matrix X) {
        X = meanNormalise(X);

        X = stDevNormalise(X);

        Matrix covariance = X.transpose().matrixMultiplication(X);

        covariance = covariance.scalarMultiplication(1.0/X.getRows());

        SVD svd = new SVD(covariance);

        Matrix U = svd.U();
        Matrix S = svd.S();
        int K = 1;
        double variation = 1.0;

        while (variation >= 0.01) {
            variation = 1.0 - S.extractMatrix(0, K, 0, K).trace() / S.trace();
            K++;
        }

        PrincipalComponents = U.extractMatrix(0, U.getRows(), 0, K);
        Matrix Z = PrincipalComponents.transpose().matrixMultiplication(X.transpose());
        compressedMatrix = PrincipalComponents.matrixMultiplication(Z).transpose();
    }

    public Matrix getPrincipalComponents() {
        return PrincipalComponents;
    }

    public Matrix getCompressedMatrix() {
        return compressedMatrix;
    }

    // Divides every value in the matrix by its standard deviation.
    private Matrix stDevNormalise(Matrix X) {
        double mean = X.mean();
        double diffTotal = 0;
        for (int i = 0; i < X.getRows(); i++) {
            for (int j = 0; j < X.getColumns(); j++) {
                diffTotal += Math.pow(X.objectAtPoint(i, j) - mean, 2);
            }
        }
        double stDev = Math.sqrt(diffTotal / X.numElements());
        return X.scalarDivision(-1 * stDev);
    }

    // Minuses the mean from every element of the matrix.
    private Matrix meanNormalise(Matrix X) {
        return X.scalarAddition(-1 * X.mean());
    }

    private Matrix compressedMatrix;
    private Matrix PrincipalComponents;
}
