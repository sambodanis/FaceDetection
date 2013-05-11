import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

// Image drawing code with inspiration from
// http://stackoverflow.com/questions/12918780/draw-accurately-on-an-image-in-java?rq=1
// I was having trouble overlaying rectangles on top of a buffered image.

public class faceLocationFrame {

    public faceLocationFrame(ArrayList<int[]> clusters, BufferedImage image) {
        this.clusters = clusters;
        this.image = image;
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(image.getWidth(), image.getHeight());
        drawRectangles(frame);
        frame.pack();
        frame.setVisible(true);
    }

    private void drawRectangles(JFrame frame) {
        faceBox box = new faceBox(image);
        frame.add(box);
        for (int[] point : clusters) {
            box.drawBox(point[1], point[0]);
        }

    }

    private BufferedImage image;
    private ArrayList<int[]> clusters;
}

class faceBox extends JPanel {

    private static final int boxEdgeLength = 60;

    public faceBox(BufferedImage image) {
        this.image = image;
    }

    public boolean drawBox(int x, int y) {
        x = (x + boxEdgeLength) / 2;
        y = (y + boxEdgeLength) / 2;
        for (int i = x; i < x + boxEdgeLength / 3; i++) {
            for (int j = 0; j < y + boxEdgeLength / 3; j++){
                image.setRGB(j, i, 0x000000);
            }
        }
        repaint();
        return true;
    }


    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2d = (Graphics2D) g;
        g2d.drawImage(image, 0, 0, null);
    }

    private BufferedImage image;
}
