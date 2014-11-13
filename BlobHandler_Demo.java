import ij.*;
import ij.plugin.filter.*;
import ij.process.*;
import ij.measure.*;
import ij.gui.*;
import ij.io.Opener;
import java.awt.*;
import ij.plugin.*;
import java.net.*; //for plugin path
import java.io.*;

/**
   * This BlobHandler demo uses an image of a set of tools from the 
   * website for Burger and Burge "The Imaging Book" 
   * http://www.imagingbook.com/index.php?id=100 
   *
   * The demo shows how to programmatically work with blobs - combine, draw, delete, etc.
   * Note that the class BlobHandler_ is *not* used.  That is because BlobHandler_
   * is simply a wrapper around BlobHandler.  Having access and control of the 
   * blobs directly is one of the reasons why we set out to recreate ImageJ's 
   * particle analyzer.  In a sense, particle analyzer 'disposes' of the necessary
   * access once it has measured each blob.  BlobHandler does not do this and
   * thus allows the programmer a chance to get in there are muck around.
   */
public class BlobHandler_Demo implements PlugIn {
   
   private String path = "";
   
   ImagePlus imp; //orig
   ImageProcessor ip;
   ImagePlus limp; //labels
   ImageProcessor lip;
   ImagePlus cimp; //contours
   ImageProcessor cip;
   BlobHandler bh ;
   int w, h;
   double combineIfCloserThan = 40.0;
   
   public void run(String arg) {
         
    long startTime = System.currentTimeMillis();
     imp = loadExample();
     IJ.run(imp, "Invert LUT", "");
     IJ.run(imp, "Fill Holes", "");
     if (imp == null) { 
      IJ.showMessage("Error loading example image");
      return;
    }
     IJ.log("\\Clear");
     IJ.log("**** BlobHandler Demo ****");
     w = imp.getWidth();
     h = imp.getHeight();

     ip = imp.getProcessor();
     
    //make the BlobHandler
      //bh = new BlobHandler(ip, imp, BHBlob.BASIC + BHBlob.SHAPE + BHBlob.LOCATION + BHBlob.ELLIPSE);
      bh = new BlobHandler(ip, imp, BHBlob.BASIC + BHBlob.SHAPE + BHBlob.LOCATION + BHBlob.ELLIPSE);
      //let it run
      bh.process();
      // get the labels of the blobs it found (1,2,3,...)
      int[]  labels = bh.blobs.getLabels(); 
      IJ.log("The 'BlobHandler Test Image' is a collage made from 8 different\n"+ 
         "binary masks of plankton imaged using FlowCam Imaging Flow Cytometer\n" +
         "See www.fluidimaging.com for details on that instrument\n\n" + 
         "The plankton show a variety of shapes and some are multi-parted.\n" + 
         "BlobHandler permits arbitrary association among the individual blobs\n\n");
         
      IJ.log("Originally there are 46 obects, but we can easily combine objects\n"+
         "within a specified distance of each other.\n" + 
         "These will be treated as one multi-parted object - aka 'CompoundBlobs'\n"+
         "The BlobHandler class will compute a simple table of the associations\n"+
         "to combine like this...\n\n");
         
      IJ.log("      double combineIfCloserThan = 40.0;\n"+
            "      int[][] table = bh.computeMergeTable(combineIfCloserThan);\n"+
            "      for (int i = 0; i < table.length; i++){\n"+
            "         if (table[i].length > 1) { bh.blobs.combine(table[i]); } \n"+
            "      }\n");
      
      //the table is a 'ragged' 2d array - each row represents a potential cluster
      // of blobs identified by their unique blob labels.
      //if the potential cluster has more than one item in it then it can be 
      //combined - the combine method utilizes the labels found to sort through 
      //the combinations
      int[][] table = bh.blobs.computeMergeTable(combineIfCloserThan);
      for (int i = 0; i < table.length; i++){
         if (table[i].length > 1) {
            bh.blobs.combine(table[i]);
         } //table shows two or more to combine
      }//setp through the table rows
    
      //remove a compund blob for any old reason
      IJ.log("\nWe can also remove blobs, in this case the laboea cell (bullet-shaped) \n"+
      " at lower right in the original image like this using the unique blob label...\n\n");
      IJ.log("    bh.blobs.removeBlob(35);\n");
      bh.blobs.removeBlob(35);
         
    //make a results table to show results
      ResultsTable rt = new ResultsTable();
      rt.reset();
      int x, y;
      BHPolygon cont;
      BHCompoundBlob blob;
      int row = 0;
      labels = bh.blobs.getLabels();
      for (int i = 0; i < labels.length; i++){
         blob = bh.blobs.get(labels[i]);
        if (blob != null){
         rt.incrementCounter(); blob.showInfo(rt, row); row++;
        } else {
          IJ.log("Blob[" + i + "]= null");
        }
      }
      rt.show("Results");

    long endTime = System.currentTimeMillis();
      
      //make masks and label images
      lip = bh.createCountMask();
      lip.resetMinAndMax();
      limp = new ImagePlus("Labeled", lip);
      limp.setCalibration(imp.getCalibration());
      limp.show();
      WindowManager.setTempCurrentImage(limp);
      IJ.run("3-3-2 RGB");
      
      //make masks and label images
      ImageProcessor cip = bh.createContourMask();
      cip.resetMinAndMax();
      ImagePlus cimp = new ImagePlus("Contours", cip);
      cimp.setCalibration(imp.getCalibration());
      cimp.show();
      WindowManager.setTempCurrentImage(cimp);
      IJ.run("3-3-2 RGB");
      
      

    IJ.log("\nShowing the convex hulls of compound objects 'Labeled' image (limp)\n\n");
 
    IJ.log("  Overlay limpOverlay = limp.getOverlay();\n"+
      "  if (limpOverlay == null) {limpOverlay = new Overlay();}\n"+
      "  limpOverlay = bh.blobs.showAllContours(BHPolygon.CHULL, limpOverlay);\n"+
      "  limp.setOverlay(limpOverlay);\n"+
      "  limp.updateAndDraw();\n");
           
      Overlay limpOverlay = limp.getOverlay();
      if (limpOverlay == null) {limpOverlay = new Overlay();}
      limpOverlay = bh.blobs.showAllContours(BHPolygon.CHULL, limpOverlay);
      limp.setOverlay(limpOverlay);
      limp.updateAndDraw();
        
    IJ.log("\nShowing the equivalent ellipse of Thallasionema cells (star-like)\nlabeled 4 on the 'Contours' image (cimp)\n\n");

    IJ.log("        blob = (BHCompoundBlob) bh.blobs.get(4);\n"+
         "        PolygonRoi ellipse = blob.getEllipseRoi();\n"+
         "        cimp.setRoi(ellipse);\n"+
         "        cimp.updateAndDraw();\n");
            
//        blob = (BHCompoundBlob) bh.blobs.get(4);
//        PolygonRoi ellipse = blob.getEllipseRoi();
//        cimp.setRoi(ellipse);
         Overlay cimpOverlay = cimp.getOverlay();
         if (cimpOverlay == null) { cimpOverlay = new Overlay();}
         cimpOverlay = bh.blobs.showAllEllipses(cimpOverlay);
         cimp.setOverlay(cimpOverlay);
         cimp.updateAndDraw();
         
      
      IJ.log("\n----  Demo complete   -----");
      String rate = IJ.d2s((w*h) * 1000.0 /(endTime-startTime));
      IJ.log("  Elapsed time = " +IJ.d2s( ((endTime-startTime) /1000.0)) + "s");
      IJ.log("  " + rate + " pixels per second");
      
   }


   // see http://rsb.info.nih.gov/ij/plugins/download/JAR_Resources_Demo.java
   // for details on loading images from the jar
   private ImagePlus loadExample() {
      ImagePlus imp = null;
      InputStream is = getClass().getResourceAsStream(path+"extras/test.tif");
      if (is!=null) {
          Opener opener = new Opener();
          imp = opener.openTiff(is, "BlobHandler Test Image");
          if (imp!=null) imp.show();
      }
      return imp;
   }

    
}
