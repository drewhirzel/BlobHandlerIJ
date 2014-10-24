import ij.IJ;
import ij.plugin.PlugIn;


public class About_BlobHandler implements PlugIn {

  public void run(String arg) { 
   IJ.showMessage("BlobHandler\n" +
    " \n" +
    "BlobHandler provides particle analysis similar to that of ImageJ's \n" +
    "built-in Particle Analyzer.\n" +
    " \n" +
    "The implementation was done by Ben Tupper and Mike Sieracki at\n" +
    "Bigelow Laboratory for Ocean Science\n"+
    "West Boothbay Harbor, Maine 04575\n" + 
    "www.bigelow.org\n" + 
    " \n" +
    "See the documentation included with the distribution for detailed information.\n");  
   }
    
}