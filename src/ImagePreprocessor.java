import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.io.File;
import java.io.IOException;

public class ImagePreprocessor {

    public ImagePreprocessor(String file) {
        buffImage = null;

        try {
            buffImage = ImageIO.read(new File(file));
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
        buffImage = toGrayScale(buffImage);
        int[] pixels = getImagePixels(buffImage);
        image = matrixFromPixelArray(pixels, buffImage.getHeight(), buffImage.getWidth());

    }

    public void displayImage(Matrix im) {
        BufferedImage bufferedImage = new BufferedImage(im.getColumns(), im.getRows(), BufferedImage.TYPE_INT_RGB);
        for (int i = 0; i < im.getRows(); i++) {
            for (int j = 0; j < im.getColumns(); j++) {
                int currPixel = (int)im.objectAtPoint(i,j) << 16 |
                        (int)im.objectAtPoint(i,j) << 8 |
                        (int)im.objectAtPoint(i,j);
                bufferedImage.setRGB(j,i, currPixel);
            }
        }

        JFrame frame = new JFrame();
        frame.getContentPane().setLayout(new FlowLayout());
        frame.getContentPane().add(new JLabel(new ImageIcon(bufferedImage)));
        frame.pack();
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    public BufferedImage toGrayScale(BufferedImage startPic) {
        BufferedImage image = new BufferedImage(startPic.getWidth(), startPic.getHeight(),
                BufferedImage.TYPE_BYTE_GRAY);
        Graphics g = image.getGraphics();
        g.drawImage(startPic, 0, 0, null);
        g.dispose();

        return image;
    }

    private Matrix matrixFromPixelArray(int[] pixels, int height, int width) {
        Matrix result = new Matrix(new double[height][width]);
        int k = 0;
        for (int i = 0; i < result.getRows(); i++) {
            for (int j = 0; j < result.getColumns(); j++) {
                result.setObjectAtPoint(i, j, pixels[k++]);
            }
        }
        return result;
    }

    private int[] getImagePixels(BufferedImage img) {
        int [] temp = null;
        Raster pixelData = img.getData();
        return pixelData.getPixels(0, 0, img.getWidth(), img.getHeight(), temp);
    }

    public Matrix getImageData() {
        return this.image;
    }

    public BufferedImage getBuffImage() {
        return buffImage;
    }

    private Matrix image;
    private BufferedImage buffImage;

}
