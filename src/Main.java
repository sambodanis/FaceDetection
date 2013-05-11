import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class Main implements Runnable {

    public static void main(String[] args) {

        Main program = new Main();
        SwingUtilities.invokeLater((Runnable) program);
    }

    public void run() {
//        testMatrices();
//        testSVD();

        ArrayList<Matrix> faces = getFaceImageMatrices(".");
        String imageName = "people4test.jpg";
        ImagePreprocessor imProcessor2 = new ImagePreprocessor(imageName);
        int numFaces = 4;
        Matrix pictureData = imProcessor2.getImageData();
//        pictureData = pictureData.extractMatrix(0, pictureData.getRows() / 2, 0, pictureData.getColumns() / 2);
        ArrayList<int[]> faceLocations = findFacesInPicture(faces, pictureData);
        KMeansClustering kMeans = new KMeansClustering(numFaces, faceLocations);
        ArrayList<int[]> clusters = kMeans.getCentroids();

        for (int[] c : clusters) {
            System.out.println("Cluster: " + c[0] + ", " + c[1]);
        }
        faceLocationFrame faceDisplay = new faceLocationFrame(clusters, imProcessor2.getBuffImage());

    }

    // Reads in all the face image files,
    // converts them to matrices and stores
    // them in an arrayList of matrices.
    public ArrayList<Matrix> getFaceImageMatrices(String directory) {
        ArrayList<Matrix> result = new ArrayList<Matrix>();
        List<String> trainingImageFiles = filesOfType(directory, ".jpg");
        for (String filename : trainingImageFiles) {
            System.out.println("a " + filename);
            if (filename.contains("yale")) {
                ImagePreprocessor im = new ImagePreprocessor(filename);
                result.add(im.getImageData());
            }
        }
        return result;
    }

    // Returns a list of all files of a specific type within a directory.
    private List<String> filesOfType(String directory, String type) {
        List<String> textFiles = new ArrayList<String>();
        File dir = new File(directory);
        for (File file : dir.listFiles()) {
            if (file.getName().endsWith(type)) {
                textFiles.add(file.getName());
            }
        }
        return textFiles;
    }


    // Uses a sliding window technique to scan through a picture
    // at every step, it gets the principal components from the submatrix
    // within the window and calculates its distance from each
    // face image and takes the average.
    // If this distance is below a certain threshold, there is a
    // good chance that there is a face within the window's submatrix.
    private ArrayList<int[]> findFacesInPicture(ArrayList<Matrix> facesRaw, Matrix pictureData) {
        int scaleEdgeSize = 30;
        ArrayList<Matrix> faces = new ArrayList<Matrix>();
        for (Matrix mat : facesRaw) {
            Matrix compress = mat.scaleDown(30);
            PCA pca = new PCA(compress);
            faces.add(pca.getPrincipalComponents());
        }
        ArrayList<int[]> result = new ArrayList<int[]>();

        int windowLength = 2 * scaleEdgeSize;
        for (int i = 0; i < pictureData.getRows() - windowLength; i += (2.0 / 3) * windowLength) {
            for (int j = 0; j < pictureData.getColumns() - windowLength; j += 10) {
                Matrix testArea = pictureData.extractMatrix(i, i + windowLength, j, j + windowLength);
                testArea = testArea.scaleDown(scaleEdgeSize);
                if (testArea.mean() == 0) continue;
                PCA pcaTemp = new PCA(testArea);
                Matrix testAreaPCs = pcaTemp.getPrincipalComponents();
                double difference = matrixDifference(faces, testAreaPCs);
                if (difference < 25) {
                    int[] newFacePlace = {j, i};
                    result.add(newFacePlace);
//                    System.out.printf("Min difference (%.1f) at %d,%d\n", difference, j, i);
                }
            }
        }
        return result;
    }

    // Calculates the distance between matrices.
    private double matrixDifference(ArrayList<Matrix> faces, Matrix B) {
        double totalDistance = 0;
        for (Matrix A : faces) {
            for (int i = 0; i < A.getRows(); i++) {
                for (int j = 0; j < Math.min(A.getColumns(), B.getColumns()); j++) {
                    totalDistance += Math.abs(A.objectAtPoint(i, j) - B.objectAtPoint(i, j));
                }
            }
        }
        return totalDistance / (double)faces.size();
    }

    private void testSVD() {
        double[][] t1 = {{104, 8,   90,  108, 0},
                {8,   87,  9,   12,  109},
                {90,  9,   90,  111, 0},
                {108, 12,  111, 138, 0},
                {0,   109, 0,   0,   149}};

        SVD svd = new SVD(new Matrix(t1));
        System.out.println("U");
        svd.U().displayMatrix();
        System.out.println("\nS");
        svd.S().displayMatrix();
        System.out.println("\nV");
        svd.V().displayMatrix();
        System.out.println("Tests:\nUtU = I");
        svd.U().matrixMultiplication(svd.U().transpose()).displayMatrix();
        System.out.println("VtV = I");
        svd.V().matrixMultiplication(svd.V().transpose()).displayMatrix();
        System.out.println("A = USVt");
        svd.U().matrixMultiplication(svd.S()).matrixMultiplication(svd.V().transpose()).displayMatrix();


    }

    private void testMatrices() {
        double[][] testMatrixInternal = {{1, -1, 4}, {1, 4, -2}, {1, 4, 2}, {1, -1, 0}};
        Matrix testMatrix = new Matrix(testMatrixInternal);
        testMatrix.displayMatrix();
        System.out.println("maxAbsElem: " + (testMatrix.maxAbsElem() == 4));
        System.out.println("NumElements: " + (testMatrix.numElements() == 12));
        System.out.println("Extract: " + (testMatrix.extractMatrix(0, testMatrix.getRows(),
                0, testMatrix.getColumns()).equals1(testMatrix)));
        double[][] testTemp = {{1}, {1}, {1}, {1}};
        System.out.println("Extract 2: " + testMatrix.extractMatrix(0, 4, 0, 1).equals1(new Matrix(testTemp)));

        double[][] t4 = {{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}};

        Matrix t5 = new Matrix(t4);

        System.out.println("Eigenvalues: ");
        t5.eigenvectors().displayMatrix();

        System.out.println("Eigenvectors");
        t5.eigenvectors().displayMatrix();

        double[][] t9 = {{104, 8,   90,  108, 0},
                         {8,   87,  9,   12,  109},
                         {90,  9,   90,  111, 0},
                         {108, 12,  111, 138, 0},
                         {0,   109, 0,   0,   149}};
        new Matrix(t9).reverseColumns().displayMatrix();

        QRDecomposition QRObj3 = new QRDecomposition(new Matrix(t9));

        System.out.println("Q: ");
        QRObj3.Q().displayMatrix();
        System.out.println("Q'Q = I: ");
        QRObj3.Q().transpose().matrixMultiplication(QRObj3.Q()).displayMatrix(); // Q' Q = I
        System.out.println("R: ");
        QRObj3.R().displayMatrix();
        System.out.println("R = Q'A");
        QRObj3.Q().transpose().matrixMultiplication(new Matrix(t4)).displayMatrix(); // R = Q' A

    }


}
