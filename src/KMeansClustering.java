import java.util.ArrayList;
import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 * User: sambodanis
 * Date: 05/05/2013
 * Time: 17:38
 * To change this template use File | Settings | File Templates.
 */
public class KMeansClustering {

    public KMeansClustering(int clusters, ArrayList<int[]> faceLocations) {
        centroids = new ArrayList<int[]>();
        Random rgen = new Random();
        for (int i = 0; i < clusters; i++) {
            int[] centroidStart = faceLocations.get(rgen.nextInt(faceLocations.size()));
            while (centroids.contains(centroidStart)) {
                centroidStart = faceLocations.get(rgen.nextInt(faceLocations.size()));
            }
            centroids.add(centroidStart);
        }
        for (int i = 0; i < 20; i++) {
            ArrayList<int[]> facesAndCentroids = assignClustersToCentroids(faceLocations);
            moveCentroids(facesAndCentroids);
        }
    }

    // Moves each centroid to the average position of
    // all points belonging to it.
    private void moveCentroids(ArrayList<int[]> facesAndCentroids) {
        ArrayList<int[]> newCentroids = new ArrayList<int[]>();
        for (int i = 0; i < centroids.size(); i++) {
            int numAssignedToCentroid = 1;
            int sumX = 0;
            int sumY = 0;
            for (int[] point : facesAndCentroids) {
                if (point[2] == i) {
                    numAssignedToCentroid++;
                    sumX += point[0];
                    sumY += point[1];
                }
            }
            int[] newCentroidLocation = {sumX / numAssignedToCentroid, sumY / numAssignedToCentroid};
            newCentroids.add(newCentroidLocation);
        }
        centroids = newCentroids;
    }

    // Assigns each point to the centroid nearest to it.
    private ArrayList<int[]> assignClustersToCentroids(ArrayList<int[]> faceLocations) {
        ArrayList<int[]> facesAndCentroids = new ArrayList<int[]>();
        for (int[] coord : faceLocations) {
            double minDist = Double.MAX_VALUE;
            int minDIndex = 0;
            for (int i = 0; i < centroids.size(); i++) {
                int[] centroid = centroids.get(i);
                if (Matrix.distanceBetweenPoints(centroid[0], centroid[1], coord[0], coord[1]) < minDist) {
                    minDist = Matrix.distanceBetweenPoints(centroid[0], centroid[1], coord[0], coord[1]);
                    minDIndex = i;
                }
            }
            int[] coordAndCentroid = {coord[0], coord[1], minDIndex};
            facesAndCentroids.add(coordAndCentroid);
        }
        return facesAndCentroids;
    }

    public ArrayList<int[]> getCentroids() {
        return centroids;
    }

    private ArrayList<int[]> centroids;
}
