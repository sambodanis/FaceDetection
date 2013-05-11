/**
 * This is a simple matrix library that supports most
 * useful matrix functions along with some other specialised
 * functions for selecting columns, collapsing matrices and
 * emulating writing.
 */

import java.awt.*;
import java.util.ArrayList;
import java.util.Random;

public class Matrix {

    private static final String hadamardProduct = "hadamardProduct";
    private static final String elementwiseDivision = "elementwiseDivision";
    private static final String elementwiseLog = "elementwiseLog";
    private static final String matrixAddition = "elementwiseAddition";
    private static final String matrixMultiplication = "matrixMultiplication";
    private static final String scalarAddition = "scalarAddition";
    private static final String scalarMultiplication = "scalarMultiplication";
    private static final String nMinusMatrix = "nMinusMatrix";
    private static final int newEdgeLength = 50;
    private static final int maxDrawingEmulationRadius = 5;

    // Sets matrix to null unless a 2d array is passed in.
    public Matrix(double[][] arrayForm) {
        matrixInternal = null;
        rows = 0;
        columns = 0;
        if (arrayForm != null) {
            initMatrix(arrayForm);
        }
    }

    // The user's input drawing is too clean, just a single path of 255s.
    // There is not enough variety in inputs to distinguish between close classes
    // This function adds some more weight to each matrix.
    // This function normally distributes points around a point within a random radius.
    public void emulateWriting() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (matrixInternal[i][j] == 255) {
                    drawAroundPoint(i, j);
                }
            }
        }
    }

    // Reorganises a matrix so that rows are arranged end to end.
    public Matrix matrixToRowVector() {
        Matrix result = new Matrix(new double[1][rows * columns]);
        int k = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result.setObjectAtPoint(0, k++, matrixInternal[i][j]);
            }
        }
        return result;
    }

    public boolean equals1(Matrix m) {
        if (!m.hasEqualSizeTo(this)) return false;
        for (int i = 0; i < this.getRows(); i++) {
            for (int j = 0; j < this.getColumns(); j++) {
                if (this.objectAtPoint(i, j) != m.objectAtPoint(i, j)) return false;
            }
        }
        return true;
    }

    // Turns a single row vector into a matrix.
    public Matrix rowVectorToMatrix(int row) {
        Matrix result = new Matrix(new double[(int)Math.sqrt(columns)][(int)Math.sqrt(columns)]);
        int k = 0;
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                result.setObjectAtPoint(i, j, matrixInternal[row][k++]);
            }
        }
        return result;
    }

    // Scales a matrix down to a specified size.
    // Method: split the old matrix up into n blocks
    // where n is the new matrix's edge size.
    // Average all values within a block and set
    // that as the value for that element of the new matrix.
    public Matrix scaleDown(int newEdgeLength) {
        Matrix result = new Matrix(new double[newEdgeLength][newEdgeLength]);
        int rowStepSize = (int)(rows / newEdgeLength + 0.5);
        int columnStepSize = (int) (columns / newEdgeLength + 0.5);
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                result.setObjectAtPoint(i, j, result.objectAtPoint(i, j) +
                        averageAroundPoint(i * rowStepSize, j * columnStepSize,
                                rowStepSize, columnStepSize));
            }
        }
        return result;
    }

    // Computes the average pixel values of a section of a matrix.
    private double averageAroundPoint(int i, int j, int rowStepSize, int columnStepSize) {
        double result = 0;
        for (int k = i - rowStepSize / 2; k < i + rowStepSize / 2; k++) {
            if (k < rowPadding || k > rows - rowPadding) continue;
            for (int l = j - columnStepSize / 2; l < j + columnStepSize / 2; l++) {
                if (l < colPadding || l > columns - colPadding) continue;
                if (isInBoundsOfMatrix(k, l)) {
                    result += normalDistribution(i, j, k, l) * matrixInternal[k][l] * rowStepSize * columnStepSize;
                }
            }
        }
        return result / (rowStepSize * columnStepSize);
    }

    public boolean isInBoundsOfMatrix(int i, int j) {
        return i >= 0 && i < rows && j >= 0 && j < columns;
    }

    // draws around a point while dropping off with the normal distribution.
    // if the point to be drawn to is not an original part of the shape,
    // its value is averaged with the result of the normal distribution.
    private void drawAroundPoint(int i, int j) {
        rGen = new Random();
        int radius = maxDrawingEmulationRadius;

        for (int k = i - radius; k < i + radius; k++) {
            for (int l = j - radius; l < j + radius; l++) {
                if (isInBoundsOfMatrix(k, l)) {
                    if (matrixInternal[k][l] != 255) {
                        matrixInternal[k][l] += normalDistribution(i, j, k, l)
                                + 1.0 / distanceBetweenPoints(i, j, k, l);
                    }
                }
            }
        }
    }

    // euclidean distance between two points.
    public static double distanceBetweenPoints(int p1x, int p1y, int p2x, int p2y) {
        double xDiff = p2x - p1x;
        double yDiff = p2y - p1y;
        double distance = Math.pow(xDiff, 2.0) + Math.pow(yDiff, 2.0);
        return Math.sqrt(distance);
    }

    // Normal distribution formula.
    private double normalDistribution(int p1x, int p1y, int p2x, int p2y) {
        double distance = distanceBetweenPoints(p2x, p2y, p1x, p1y);
        double sigma = 1.0;
        double leftSide = 1.0 / (sigma * 2.0 * Math.sqrt(2.0 * Math.PI));
        double rightSide = -distance / (2.0 * Math.pow(sigma, 2.0));
        return leftSide * Math.exp(rightSide);
    }

    // Prints out a matrix for debugging on the command line.
    // Maybe cast it to an int for uncluttered viewing
    public void displayMatrix() {
        System.out.println();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
//                if (matrixInternal[i][j] != 255) {
                    System.out.printf("%.2f  ", matrixInternal[i][j]);
//                } else {
//                    System.out.printf("   ");
//                }
            }
            System.out.println();
        }
    }

    // Displays the nth row of a matrix
    public void displayRowN(int n) {
        for (int j = 0; j < columns; j++) {
            System.out.printf("%f\n", matrixInternal[n][j]);
        }
    }

    public int numElements() {
        return this.rows * this.columns;
    }


    // For some reason the eigenvectors come out backwards so in the final step the
    // columns are reversed.
    public Matrix eigenvectors() {
        Matrix B = this;
        Matrix QTemp = null;
        for (int i = 0; i < 30; i++) {
            QRDecomposition qr = new QRDecomposition(B);
            if (i == 0) {
                QTemp = qr.Q();
            } else {
                QTemp = QTemp.matrixMultiplication(qr.Q());
            }
            B = qr.R().matrixMultiplication(qr.Q());
        }
        assert QTemp != null;
        return QTemp.reverseColumns();

//        Matrix eigs = this.eigenvalues();
//        this.eigenvalues().displayMatrix();
////        eigs = new Matrix(test);
//        Matrix result = new Matrix(new double[eigs.getRows()][eigs.getRows()]);
//
//        for (int i = 0; i < eigs.getRows(); i++) {
//            double currentEigenvalue = eigs.objectAtPoint(i, 0);
//            Matrix thisCopy = this.copy();
//
//            Matrix eigenvec = thisCopy.matrixAddition(
//                    (Matrix.identityMatrix(this.getRows(), this.getColumns())
//                            .scalarMultiplication(-1 * currentEigenvalue))); // (A - L .* I)
//
//            eigenvec = eigenvec.cramer(new Matrix(new double[this.getRows()][1])); // solve (A - L .* I)v = 0
//            result = result.insertMatrixIntoThis(0, i, eigenvec);
////            result = result.insertMatrixIntoThis(result, 0, i, eigenvec);
//
//        }
//        return result;
    }

    public Matrix reverseColumns() {
        Matrix result = new Matrix(new double[this.columns][this.rows]);

        for (int i = 0; i < this.getColumns(); i++) {
            int k = result.getColumns() - 1;
            for (int j = 0; j <= k; j++, k--) {
                double left = this.objectAtPoint(i, j);
                result.setObjectAtPoint(i, j, this.objectAtPoint(i, k));
                result.setObjectAtPoint(i, k, left);
            }
        }

        return result;
    }

    // Returns a column vector of the matrix's eigenvalues.
    public Matrix eigenvalues() {
        Matrix B = this;
        for (int i = 0; i < 30; i++) {
            QRDecomposition qr = new QRDecomposition(B);
            B = qr.R().matrixMultiplication(qr.Q());
        }

        Matrix result = new Matrix(new double[this.rows][1]);
        int i = 0;
        int j = 0;
        for (int k = 0; k < result.getRows(); k++) {
            result.setObjectAtPoint(k, 0, B.objectAtPoint(i++, j++));
        }
        return result;
    }

    // Finds the pseudo inverse of a matrix A where:
    // A = USV, pinv(A) = V * S.reciprocalTranspose * U.transpose
    public Matrix pseudoInverse() {
        SVD svd = new SVD(this);
        Matrix newS = new Matrix(svd.S().getMatrix());
        for (int i = 0; i < svd.S().getRows(); i++) {
            for (int j = 0; j < svd.S().getColumns(); j++) {
                if (i == j && svd.S().objectAtPoint(i, j) != 0) {
                    newS.setObjectAtPoint(i, j, 1 / svd.S().objectAtPoint(i, j));
                }
            }
        }
        return svd.V().matrixMultiplication(newS.transpose()).matrixMultiplication(svd.U().transpose());
    }

    public Matrix cramer(Matrix B) {
        Matrix vars = new Matrix(new double[this.getRows()][1]);
        for (int i = 0; i < vars.rows; i++) {
//            Matrix matReplacedAtIndex = this.insertMatrixIntoThis(0, i, B);
            Matrix matReplacedAtIndex = this.copy().insertMatrixIntoThis(0, i, B);
//            Matrix matReplacedAtIndex = this.insertMatrixIntoThis(this, 0, i, B);
//            matReplacedAtIndex.insertMatrixIntoThis(matReplacedAtIndex, 0, i, B);

//            System.out.println("Top: " + matReplacedAtIndex.determinant());
//            System.out.println("Bottom: " + this.determinant());
//            matReplacedAtIndex.displayMatrix();
//            System.out.println("------------");
//            this.displayMatrix();
//            System.out.println("\n");
            double varSolution = matReplacedAtIndex.determinant() / this.determinant();
            vars.setObjectAtPoint(i, 0, varSolution);
        }
        return vars;
    }

    public double determinant() {
        Matrix eig = this.eigenvalues();
//        eig.displayMatrix();
        double eigProduct = 1;
        for (int i = 0; i < eig.getRows(); i++) {
            eigProduct *= eig.objectAtPoint(i, 0);
        }
        return eigProduct;
    }

    public static Matrix projection(Matrix e, Matrix a) {
        Matrix top = e.transpose().matrixMultiplication(a);
        Matrix bottom = e.transpose().matrixMultiplication(e);
        return bottom.objectAtPoint(0,0) != 0 ?
                e.scalarMultiplication(top.objectAtPoint(0,0) /
                                       bottom.objectAtPoint(0,0))
                : new Matrix(e.getMatrix());
    }

    public static Matrix normalise(Matrix e) {
        double sumSquares = 0;
        for (int i = 0; i < e.getRows(); i++) {
            for (int j = 0; j < e.getColumns(); j++) {
                sumSquares += Math.pow(e.objectAtPoint(i, j), 2);
            }
        }

        return e.scalarDivision(Math.sqrt(sumSquares));
    }


    public double trace() {
        if (this.rows != this.columns) throw new Error("Cannot perform trace on non square matrix: "
                + this.rows + "," + this.columns);
        int sumResult = 0;
        for (int i = 0; i < this.rows; i++) {
            sumResult += this.objectAtPoint(i,i);
        }
        return sumResult;
    }

    public double mean() {
        double elemSum = 0;
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                elemSum += this.objectAtPoint(i,j);
            }
        }
        return elemSum / this.numElements();
    }

    public Matrix copy() {
        Matrix result = new Matrix(new double[this.getRows()][this.getColumns()]);
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                result.setObjectAtPoint(i, j, this.objectAtPoint(i, j));
            }
        }
        return result;
    }

    public Matrix insertMatrixIntoThis(int iRow, int iCol, Matrix matrixToInsert) {
        Matrix result = this.copy();
        int k = 0;
        int l = 0;
        for (int i = iRow; k < matrixToInsert.getRows(); i++) {
            for (int j = iCol; l < matrixToInsert.getColumns(); j++) {
                double valToInsert = matrixToInsert.objectAtPoint(k, l++);
                result.setObjectAtPoint(i, j, valToInsert);
            }
            k++;
            l = 0;
        }
        return result;
    }

    public double frobeniusNormal() {
        int resultSum = 0;
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.columns; j++) {
                resultSum += Math.pow(this.objectAtPoint(i,j), 2);
            }
        }
        return Math.sqrt(resultSum);
    }

    // Returns the maximum abs value of all elements in the matrix
    public double maxAbsElem() {
        double max = 0;
        for (int i = 0; i < this.getRows(); i++) {
            for (int j = 0; j < this.getColumns(); j++) {
                if (Math.abs(this.objectAtPoint(i,j)) > max) max = Math.abs(this.objectAtPoint(i,j));
            }
        }
        return max;
    }

    // Extracts a sub-matrix from a full matrix.
    public Matrix extractMatrix(int rowStart, int rowEnd, int columnStart, int columnEnd) {
        Matrix result = new Matrix(new double[rowEnd - rowStart][columnEnd - columnStart]);
        for (int i = rowStart; i < rowEnd; i++) {
            for (int j = columnStart; j < columnEnd; j++) {
                result.setObjectAtPoint(i - rowStart,j - columnStart, this.objectAtPoint(i,j));
            }
        }
        return result;
    }

    // Creates a matrix out of a 2d array.
    private void initMatrix(double[][] arrayForm) {
        rows = arrayForm.length;
        columns = arrayForm[0].length;
        matrixInternal = arrayForm;
    }

    // Gets the element at a certain index of a matrix.
    public double objectAtPoint(int row, int col) {
        if (row < 0 || row > rows || col < 0 || col > columns) {
            System.err.println("Index outside matrix boundaries: " + row + ", " + col);
        }
        return this.matrixInternal[row][col];
    }

    // Sets the element at a certain index of a matrix.
    public void setObjectAtPoint(int row, int col, double value) {
        if (row < 0 || row > rows || col < 0 || col > columns) {
            System.err.println("Index outside matrix boundaries: " + row + ", " + col);
        }
        matrixInternal[row][col] = value;
    }

    // Returns I, the identity matrix with 1s along the diagonal.
    public static Matrix identityMatrix(int nRows, int nColumns) {
        Matrix result = new Matrix(new double[nRows][nColumns]);
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                if (i == j) {
                    result.setObjectAtPoint(i, j, 1);
                }
            }
        }
        return result;
    }

    // Simply appends one matrix to another.
    public Matrix stickToMatrix(Matrix snd) {
        if (snd.getRows() != rows) {
            System.err.println("Second matrix rows (" + snd.getRows() + ") not equal to first matrix rows (" + rows + ")");
            return null;
        }
        Matrix result = new Matrix(new double[rows][columns + snd.getColumns()]);
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                if (j < columns) {
                    result.setObjectAtPoint(i, j, matrixInternal[i][j]);
                } else {
                    result.setObjectAtPoint(i, j, snd.objectAtPoint(i, j - columns));
                }
            }
        }
        return result;
    }

    // Returns a single specified column of a greater matrix.
    public Matrix extractColumn(int index) {
        if (index < 0 || index > columns) {
            System.err.println("Specified index (" + index + ") is greater than num columns (" + columns + ")");
            return null;
        }
        Matrix result = new Matrix(new double[rows][1]);
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < this.getColumns(); j++) {
                if (j == index) {
                    result.setObjectAtPoint(i, 0, matrixInternal[i][j]);
                }
            }
        }
        return result;
    }

    // Removes the Nth column from a matrix.
    public Matrix cutColumnN(int n) {
        if (n < 0 || n > columns) {
            System.err.println("Specified index (" + n + ") is greater than num columns (" + columns + ")");
            return null;
        }
        Matrix result = new Matrix(new double[rows][columns - 1]);
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                if (j < n) {
                    result.setObjectAtPoint(i, j, matrixInternal[i][j]);
                } else {
                    result.setObjectAtPoint(i, j, matrixInternal[i][j + 1]);
                }
            }
        }
        return result;
    }

    // Prepends a column of this.getRows() elements of the specified value to the matrix.
    public Matrix prependColumnOfValue(int value) {
        Matrix result = new Matrix(new double[rows][columns + 1]);
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                if (j == 0) {
                    result.setObjectAtPoint(i, j, value);
                } else {
                    result.setObjectAtPoint(i, j, matrixInternal[i][j - 1]);
                }
            }
        }
        return result;
    }

    // Creates a row vector where each element is the sum
    // of each column in the original matrix.
    public Matrix sigmaSumColumnsToRowVector() {
        Matrix result = new Matrix(new double[1][columns]);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result.setObjectAtPoint(0, j, result.objectAtPoint(0, j) + matrixInternal[i][j]);
            }
        }
        return result;
    }

    // Transposes a matrix, A'(i,j) == A(j,i).
    public Matrix transpose() {
        Matrix result = new Matrix(new double[columns][rows]);
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                result.setObjectAtPoint(i, j, matrixInternal[j][i]);
            }
        }
        return result;
    }

    // Performs elementwise multiplication rather than full matrix multiplication.
    public Matrix hadamardProduct(Matrix sndMatrix) {
        if (!this.hasEqualSizeTo(sndMatrix)) {
            return null;
        }
        return matrixElementWiseFunction(sndMatrix, hadamardProduct, 0);
    }

    // Performs elementwise division
    public Matrix elementwiseDivision(Matrix sndMatrix) {
        if (!this.hasEqualSizeTo(sndMatrix)) {
            return null;
        }
        return matrixElementWiseFunction(sndMatrix, elementwiseDivision, 0);

    }

    // Computes the log of every element of the input matrix.
    public Matrix elementwiseLog() {
        return matrixElementWiseFunction(null, elementwiseLog, 0);
    }

    // Adds the (i,j)th elements of two matrices together
    public Matrix matrixAddition(Matrix sndMatrix) {
        if (!this.hasEqualSizeTo(sndMatrix)) {
            return null;
        }
        return matrixElementWiseFunction(sndMatrix, matrixAddition, 0);
    }

    // Regular matrix multiplication
    public Matrix matrixMultiplication(Matrix sndMatrix) {
        if (this.columns != sndMatrix.getRows()) {
            System.out.println("Incompatable matrix size: " + this.rows + "*" + this.columns +
                    " x " + sndMatrix.getRows() + "*" + sndMatrix.getColumns());
            return null;
        }
        Matrix result = new Matrix(new double[rows][sndMatrix.getColumns()]);
        for (int i = 0; i < result.getRows(); i++) {
//            System.out.println("I = " + i);
            for (int j = 0; j < result.getColumns(); j++) {
//                System.out.println("i,j = " + i + "," + j);
                for (int k = 0; k < columns; k++) {
                    try {
                        double newValue = result.objectAtPoint(i, j) + matrixInternal[i][k] * sndMatrix.objectAtPoint(k, j);
                        result.setObjectAtPoint(i, j, newValue);
                    } catch (ArrayIndexOutOfBoundsException e) {
                        e.printStackTrace();
                        printSimp("Failed with this: " + rows + "," + columns + " snd: " +
                                sndMatrix.getRows() + "," + sndMatrix.getColumns() + " failIndices this: " + i + "," + k +
                                " snd: " + k + "," + j);
                        System.exit(1);
                    }
                }
            }
        }

        return result;
    }

    // Scalar adds every element of the matrix
    public Matrix scalarAddition(double scalar) {
        return matrixElementWiseFunction(null, scalarAddition, scalar);
    }

    // Scalar multiplies every element of the matrix
    public Matrix scalarMultiplication(double scalar) {
        return matrixElementWiseFunction(null, scalarMultiplication, scalar);
    }

    public Matrix scalarDivision(double scalar) {
        return scalar != 0 ? scalarMultiplication(1/scalar) : new Matrix(this.getMatrix());
    }

    // Minuses the matrix value i,j from a scalar at each index.
    public Matrix nMinusMatrix(double scalar) {
        return matrixElementWiseFunction(null, nMinusMatrix, scalar);

    }

    // Checks that 2 matrices are the same size for elementwise operations.
    private boolean hasEqualSizeTo(Matrix sndMatrix) {
        if (this.rows != sndMatrix.getRows() || this.columns != sndMatrix.getColumns()) {
            System.out.println("Matrix dimensions incompatible: " + this.rows + " != " +
                    sndMatrix.getRows() + " || " + this.columns + " != " + sndMatrix.getColumns());
            return false;
        }
        return true;
    }

    // Shorthand printing function
    private void printSimp(String str) {
        System.out.println(str);
    }

    // Performs specified matrix function.
    private Matrix matrixElementWiseFunction(Matrix sndMatrix, String function, double scalar) {
        double[][] resultArr = new double[rows][columns];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (function.equals(matrixAddition)) {
                    resultArr[i][j] = resultArr[i][j] +
                            matrixInternal[i][j] + sndMatrix.objectAtPoint(i, j);

                } else if (function.equals(hadamardProduct)) {

                    resultArr[i][j] = resultArr[i][j] +
                            matrixInternal[i][j] * sndMatrix.objectAtPoint(i, j);

                } else if (function.equals(elementwiseDivision) &&
                        matrixInternal[i][j] != 0 &&
                        sndMatrix.objectAtPoint(i, j) != 0) {

                    resultArr[i][j] = resultArr[i][j] +
                            matrixInternal[i][j] / sndMatrix.objectAtPoint(i, j);

                } else if (function.equals(elementwiseLog) && matrixInternal[i][j] != 0) {

                    resultArr[i][j] = Math.log(matrixInternal[i][j]);

                } else if (function.equals(scalarAddition)) {

                    resultArr[i][j] = matrixInternal[i][j] + scalar;

                } else if (function.equals(scalarMultiplication)) {

                    resultArr[i][j] = matrixInternal[i][j] * scalar;

                } else if (function.equals(nMinusMatrix)) {

                    resultArr[i][j] = scalar - matrixInternal[i][j];

                }
            }
        }
        return new Matrix(resultArr);
    }

    // Returns all internal matrix values.
    public double[][] getMatrix() {
        return matrixInternal;
    }

    public int getRows() {
        return rows;
    }

    public int getColumns() {
        return columns;
    }

    // From an image's bounding box, computes how many rows and columns a matrix should have.
    private void computeRC(Point[] boundingBox) {
        columns = boundingBox[1].x - boundingBox[0].x + 1;
        rows = boundingBox[0].y - boundingBox[1].y + 1;
    }

    // Given a String list of elements of a matrix arranged (i_0, j_0), (i_0, j_1), ... ,,(i_m, j_n)
    // Turns it into an matrix of size m*n.
    public Matrix makeMatrixFromStringArrayList(ArrayList<String> allLines) {
        double[][] result = new double[allLines.size()][allLines.get(0).split(" ").length];

        for (int i = 0; i < result.length; i++) {
            String[] currentRow = allLines.get(i).split(" ");
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] = (Double.parseDouble(currentRow[j]));
            }
        }
        return new Matrix(result);
    }

    // Given an arrayList of cartesian coordinates, builds a matrix from those points
    // By first scaling them so that they are docked to the origin.
    // Row and Col padding are added to give matrices a border of zeros.
    public void makeMatrixFromPointsArrayList(ArrayList<Point> numberArray, Point[] boundingBox) {
        computeRC(boundingBox);
        rowPadding = 6 * rows / newEdgeLength;
        colPadding = 6 * columns / newEdgeLength;
        rows += 2 * rowPadding;
        columns += 2 * colPadding;
        matrixInternal = new double[rows][columns];
        int fillOutRadius = 10;

        for (Point elem : numberArray) {
            int xPos = rowPadding + elem.y - boundingBox[1].y;
            int yPos = colPadding + elem.x - boundingBox[0].x;
            matrixInternal[xPos][yPos] = 255;

            for (int i = xPos - fillOutRadius / 2; i < xPos + fillOutRadius / 2; i++) {
                for (int j = yPos - fillOutRadius / 2; j < yPos + fillOutRadius / 2; j++) {
                    if (isInBoundsOfMatrix(i, j)) matrixInternal[i][j] = 255;
                }
            }
        }
    }


    private double[][] matrixInternal;
    private int rows;
    private int columns;
    private int rowPadding;
    private int colPadding;
    private Random rGen;
}
