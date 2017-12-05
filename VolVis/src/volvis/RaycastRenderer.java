/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunction2DEditor.TriangleWidget;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import util.TFChangeListener;
import util.TrackballInteractor;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    public enum Method {
        SLICER, MIP, COMPOSITING, TRANSFER;
    }
    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    Method method;
    boolean shading;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
        method = Method.SLICER;
    }
    
    public void setMethod(Method method){
        this.method = method;
        changed();
    }
    
    public void applyShading(boolean set) {
        shading = set;
        changed();
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();    
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     

    short getVoxel(double[] coord, boolean interpolation) {
        if (!interpolation) {
            if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                    || coord[2] < 0 || coord[2] > volume.getDimZ()) {
                return 0;
            }

            int x = (int) Math.floor(coord[0]);
            int y = (int) Math.floor(coord[1]);
            int z = (int) Math.floor(coord[2]);

            return volume.getVoxel(x, y, z);
        } else { // tri linear interpolation
            double x = coord[0];
            double y = coord[1];
            double z = coord[2];

            if (x < 0 || x > volume.getDimX() - 1 || y < 0 || y > volume.getDimY() - 1
                    || z < 0 || z > volume.getDimZ() - 1) {
                return 0;
            }

            double xd = (x - Math.floor(x)) / (Math.ceil(x) - Math.floor(x));
            double yd = (y - Math.floor(y)) / (Math.ceil(y) - Math.floor(y));
            double zd = (z - Math.floor(z)) / (Math.ceil(z) - Math.floor(z));

            double c000 = volume.getVoxel((int) Math.floor(x), (int) Math.floor(y), (int) Math.floor(z));
            double c100 = volume.getVoxel((int) Math.ceil(x), (int) Math.floor(y), (int) Math.floor(z));

            double c001 = volume.getVoxel((int) Math.floor(x), (int) Math.floor(y), (int) Math.ceil(z));
            double c101 = volume.getVoxel((int) Math.ceil(x), (int) Math.floor(y), (int) Math.ceil(z));

            double c010 = volume.getVoxel((int) Math.floor(x), (int) Math.ceil(y), (int) Math.floor(z));
            double c110 = volume.getVoxel((int) Math.ceil(x), (int) Math.ceil(y), (int) Math.floor(z));

            double c011 = volume.getVoxel((int) Math.floor(x), (int) Math.ceil(y), (int) Math.ceil(z));
            double c111 = volume.getVoxel((int) Math.ceil(x), (int) Math.ceil(y), (int) Math.ceil(z));

            double c00 = (c000 * (1 - xd)) + (c100 * xd);
            double c01 = (c001 * (1 - xd)) + (c101 * xd);
            double c10 = (c010 * (1 - xd)) + (c110 * xd);
            double c11 = (c011 * (1 - xd)) + (c111 * xd);

            double c0 = (c00 * (1 - yd)) + (c10 * yd);
            double c1 = (c01 * (1 - yd)) + (c11 * yd);

            return (short) Math.round((c0 * (1 - zd)) + (c1 * zd));
        }
    }
    
    void slicer(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord, true);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
               
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    void mip(double[] viewMatrix, boolean moreResponsive) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        int precision;
        
        if (moreResponsive) {
            precision = 20;
        } else {
            precision = 1;
        }
        
        // Diagonal of the image
        double temp2 = (volume.getDimX() * volume.getDimX()) + (volume.getDimY() * volume.getDimY()); 
        double temp = Math.sqrt(temp2);
        double d2 = (temp * temp) + (volume.getDimZ() + volume.getDimZ());
        double d = Math.sqrt(d2);
      
        int diagonal = (int) Math.ceil(d);
        
        // sample on a plane through all slices of the volume data 
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                int val = 0; 
                
                // find the maximum value along the viewing ray
                for (int k = -(diagonal / 2); k < (diagonal / 2); k = k + precision) {   
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + viewVec[0]*k;
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + viewVec[1]*k;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + viewVec[2]*k;

                    int val2 = getVoxel(pixelCoord, true);
                    
                    if (val2 > val) {
                        val = val2;
                    }
                }
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
               
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }    
    }
    
    public void compositing(double[] viewMatrix, boolean moreResponsive) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        TFColor voxelColor;
        
        int precision;
        
        if (moreResponsive) {
            precision = 20;
        } else {
            precision = 1;
        }
        
        // Diagonal of the image
        double temp2 = (volume.getDimX() * volume.getDimX()) + (volume.getDimY() * volume.getDimY()); 
        double temp = Math.sqrt(temp2);
        double d2 = (temp * temp) + (volume.getDimZ() + volume.getDimZ());
        double d = Math.sqrt(d2);
      
        int diagonal = (int) Math.ceil(d);
  
        // sample on a plane through all slices of the volume data 
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                // Initial color along the viewing ray
                voxelColor = new TFColor(0, 0, 0, 1);
                
                // find the composite value along the viewing ray
                for (int k = -(diagonal / 2); k < (diagonal / 2); k = k + precision) {   
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + viewVec[0]*k;
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + viewVec[1]*k;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + viewVec[2]*k;

                    TFColor color = tFunc.getColor(getVoxel(pixelCoord, true));
                    voxelColor.r = (color.r * color.a) + (voxelColor.r * (1 - color.a));
                    voxelColor.g = (color.g * color.a) + (voxelColor.g * (1 - color.a));
                    voxelColor.b = (color.b * color.a) + (voxelColor.b * (1 - color.a));
                }
                
                // Shading -- with L = V
                if (shading) {
                    /*
                    // Define parameters
                    double kAmbient = 0.1;
                    double kDiff = 0.7;
                    double kSpec = 0.2;
                    double alpha = 10;
                    
                    // Compute V
                    double[] v = new double[3];
                    for (int m = 0; m < 3; m++) {
                        v[m] = viewVec[m] / VectorMath.length(viewVec);
                    }
                    
                    // Compute N
                    GradientVolume gvol = new GradientVolume(volume);
                    VoxelGradient vg = gvol.getGradient(i, j, 0);
                    
                    double[] n = new double[3];
                    n[0] = vg.x / vg.mag; 
                    n[1] = vg.y / vg.mag; 
                    n[2] = vg.z / vg.mag; 
                    
                    // Compute L
                    double[] l = v;
                    
                    // Dot product L and N
                    double dotLN = VectorMath.dotproduct(l, n);
                    
                    // Compute R
                    double[] r = new double[3];
                    
                    for (int m = 0; m < 3; m++) {
                        r[m] = (2 * dotLN * n[m]) - viewVec[m];
                    }
                    
                    double dotVR = VectorMath.dotproduct(viewVec, r);
                    
                    // Apply illumination
                    voxelColor.r = (voxelColor.r * kAmbient) + (voxelColor.r * kDiff * dotLN) + (voxelColor.r * kSpec * Math.pow(dotVR, alpha));
                    voxelColor.g = (voxelColor.g * kAmbient) + (voxelColor.g * kDiff * dotLN) + (voxelColor.g * kSpec * Math.pow(dotVR, alpha));
                    voxelColor.b = (voxelColor.b * kAmbient) + (voxelColor.b * kDiff * dotLN) + (voxelColor.b * kSpec * Math.pow(dotVR, alpha));
                    */
                }
               
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }     
    }

    
    void transfer(double[] viewMatrix, boolean moreResponsive) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        int precision;
        
        if (moreResponsive) {
            precision = 20;
        } else {
            precision = 1;
        }
        
        // Diagonal of the image
        double temp2 = (volume.getDimX() * volume.getDimX()) + (volume.getDimY() * volume.getDimY()); 
        double temp = Math.sqrt(temp2);
        double d2 = (temp * temp) + (volume.getDimZ() + volume.getDimZ());
        double d = Math.sqrt(d2);
      
        int diagonal = (int) Math.ceil(d);
        
        // sample on a plane through all slices of the volume data 
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                // Initial color along the viewing ray
                voxelColor = new TFColor(0, 0, 0, 1);
                
                // find the composite value along the viewing ray
                for (int k = -(diagonal / 2); k < (diagonal / 2); k = k + precision) {   
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + viewVec[0]*k;
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + viewVec[1]*k;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + viewVec[2]*k;

                    int val = getVoxel(pixelCoord, true);
                    
                    if (val == 0) {
                        continue;
                    }
                    
                    TFColor color = tFunc.getColor(getVoxel(pixelCoord, true));
                    
                    TransferFunction2DEditor.TriangleWidget tw = getTF2DPanel().triangleWidget;
                    short baseIntensity = tw.baseIntensity;
                    double radius = tw.radius;
                    TFColor setColor = tw.color;
                    double minMagnitude = tw.minMagnitude;
                    double maxMagnitude = tw.maxMagnitude;
                
                    VoxelGradient gradient = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]);        
                     
                    color.r = setColor.r;
                    color.g = setColor.g;
                    color.b = setColor.b;
                    
                    // gradient-based opacity weighting
                    if (minMagnitude <= gradient.mag && maxMagnitude >= gradient.mag) {
                        if (val == baseIntensity && gradient.mag == 0){
                            color.a = setColor.a;
                        } else if (gradient.mag > 0 &&
                                val - (radius * gradient.mag) <= baseIntensity && 
                                baseIntensity <= val + (radius * gradient.mag)) {
                            color.a = setColor.a * (1 - (1 / radius) * Math.abs((baseIntensity - val) / gradient.mag));
                        } else {
                            color.a = 0;
                        }
                    } else {
                        color.a = 0;
                    }
                    
                    voxelColor.r = (color.r * color.a) + (voxelColor.r * (1 - color.a));
                    voxelColor.g = (color.g * color.a) + (voxelColor.g * (1 - color.a));
                    voxelColor.b = (color.b * color.a) + (voxelColor.b * (1 - color.a));
                }          
               
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }    
    }
    
    
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {
        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        
        boolean moreResponsive = true;
        
        switch(method) {
            case SLICER:
                slicer(viewMatrix);
                break;
            case MIP:
                mip(viewMatrix, moreResponsive);
                break;
            case COMPOSITING:
                compositing(viewMatrix, moreResponsive);
                break;
            case TRANSFER:
                transfer(viewMatrix, moreResponsive);
                break;
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
