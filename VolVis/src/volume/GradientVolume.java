/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

    public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

    
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    
    private void compute() {    
        for (int i = 0; i < data.length; i++) {
            data[i] = zero;
        }
        
        // Loop over all voxels
        for (int x = 1; x < volume.getDimX() - 1; x++) {
            for (int y = 1; y < volume.getDimY() - 1; y++ ) {
                for (int z = 1; z < volume.getDimZ() - 1; z++) {
                    // Compute the gradient of the voxel
                    double xCoord = 0.5 * (volume.getVoxel(x + 1, y, z) - volume.getVoxel(x - 1, y, z));
                    double yCoord = 0.5 * (volume.getVoxel(x, y + 1, z) - volume.getVoxel(x, y - 1, z));
                    double zCoord = 0.5 * (volume.getVoxel(x, y, z + 1) - volume.getVoxel(x, y, z - 1));
            
                    setGradient(x, y, z, new VoxelGradient((float) xCoord, (float) yCoord, (float) zCoord));
                }
            }
        }
    }
    

    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                    magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
}
