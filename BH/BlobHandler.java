import ij.*;
import ij.process.*;
import ij.measure.*;
import ij.gui.*;
import java.util.Vector;
/**
  * This class implements the labeling and contouring of a bilevel image
  * as described in REF5.  This closely follows the implementation 
  * found in REF4.
  * <p>
  * After the first pass, one BHCompoundBlob class is instantiated for each unique component in 
  * the mask.  Each BHCompoundBlob provides enough information for the reconstruction of 
  * the original image plus its out contour and zero or more inner contours.
  * Here, contours are defined as the boundary pixels within the blob.  This 
  * definition is slight different than the perimeter camputed by ImageJ which
  * is akin to shrink wrapping or following the crack between pixels.
  * <p>
  * Blobs can be subsequently merged or deleleted based upon criteria such as 
  * area, proximity, etc.
  * <p>
  * Contact: Mike Sieracki, msieracki@bigelow.org <br>
  * Author:  Ben Tupper, btupper@bigelow.org <br>
  * History: Developed over the course of 2009 and 2010 using ImageJ 1.43 and Mac OSX 10.5.8 <br>
  *   at Bigelow Laboratory for Ocean Science, www.bigelow.org  */
  
 
public class BlobHandler {
  
  /** If true output meaningful info */
  public boolean verbose = true;
  /** If true, includes measurment info about the background */
  public boolean includeBackground = false;
  /** Set to the combination of mearements to make using BHBlob.BASIC + BHBlob.LOCATION + whatever...
    */
  public int measureFlag = BHBlob.BASIC;
  //public int volMethod = BHBlob.VOL_AUTO;
  private ImageProcessor bip; //binary image processor
  private ImageProcessor gip; //grayscale image processor
  private boolean hasGray = false;
  private Calibration gcal; //grayscale calibration factors
  private ImageProcessor lip; //label mask as short processor
  public int w, h; //input height and width
  private int pw, ph; //padded height and width
  private byte[][] origMask; //original as 0/1 mask with pad
  private int[][] labelMask; //label mask with pad
 
  /** Used during component labeling*/
  public static final int BACKGROUND = 0;
  /** Used during component labeling*/
  public static final int FOREGROUND = 1;  
  /** Used during component labeling*/
  public static final int VISITED = -1;
  /** The container for the generated blobs */
  public BHBlobGroup blobs = new BHBlobGroup(); //container for the blobs
  /** The region label used during component labeling */
  private int regionID = 0; //used for incrementing blob labels
  /** The slice  number (1,2,3,...) */
  public int slice = 1;
  /** for search clockwise from right */
  private int[][] delta = {{1,0},{1,1}, {0,1}, 
      {-1,1}, {-1,0}, {-1,-1}, {0,-1}, {1,-1}};
    
  private Wand wand;
  
/**
  * A dummy constructor
  */
  public BlobHandler(){
    
  }
 
/**
  * Constructs a set of grayscale blobs
  * @param thebip the bilevel mask
  * @param thegimp the associated grayscale image (8, 16 or 32 bit)
  * @param theFlags the measurement flags (see BHBlob)
  */
  public BlobHandler(ImageProcessor thebip, ImagePlus thegimp, int theFlags){
    ImageProcessor grayip = null;
    Calibration thecal = new Calibration();
    if (thegimp == null){ 
      grayip = thegimp.getProcessor();
      thecal = thegimp.getCalibration();
    }
    
    if (setup(thebip, grayip, thecal, theFlags, slice) == false) {
       return;
    }
    
  }
 
/**
  * Constructs a set of grayscale blobs using the default measurements
  * @param thebip the bilevel mask
  * @param thegimp the associated grayscale image (8, 16 or 32 bit)
  */
  public BlobHandler(ImageProcessor thebip, ImagePlus thegimp){

    ImageProcessor grayip = null;
    Calibration thecal = new Calibration();
    if (thegimp != null){ 
      grayip = thegimp.getProcessor();
      thecal = thegimp.getCalibration();
    }
    
    if (setup(thebip, grayip, thecal, measureFlag, slice) == false) {
       return;
    }
  }
 
/**
  * Constructs a set of grayscale blobs
  * @param thebip the bilevel mask
  * @param thegip the associated grayscale image (8, 16 or 32 bit)
  * @param cal the Calibration to use
  * @param theFlags the measurements flag (see BHBlob)
  */
  public BlobHandler(ImageProcessor thebip, ImageProcessor thegip, Calibration cal, int theFlags){
    if (setup(thebip, thegip, cal, theFlags, slice) == false) { 
      return;
    }
  }

/**
  * Constructs a set of grayscale blobs using the default measurements
  * @param thebip the bilevel mask
  * @param thegip the associated grayscale image (8, 16 or 32 bit)
  * @param cal the Calibration to use
  */
  public BlobHandler(ImageProcessor thebip, ImageProcessor thegip, Calibration cal){
    if (setup(thebip, thegip, cal, measureFlag, slice) == false) { 
      return;
    }
  }
    
/**
  * Constructs a set of grayscale blobs with the default measurements
  * @param thebip the bilevel mask
  * @param thegip the associated grayscale image (8, 16 or 32 bit)
  */
  public BlobHandler(ImageProcessor thebip, ImageProcessor thegip){
    if (setup(thebip, thegip, new Calibration(), measureFlag, slice) == false) { 
      return;
    }
  }
   
/**
  * Constructs a set of grayscale blobs
  * @param thebip the bilevel mask
  * @param thegip the associated grayscale image (8, 16 or 32 bit)
  * @param theFlags the measurements flag (see BHBlob)
  */
  public BlobHandler(ImageProcessor thebip, ImageProcessor thegip, int theFlags){
    if (setup(thebip, thegip, new Calibration(), theFlags, slice) == false) { 
      return;
    }
  }
  

  public boolean setup(ImageProcessor theBip, ImageProcessor theGip, Calibration theCal, int theMeasurements){
    return setup( theBip,  theGip,  theCal,  theMeasurements, slice);
  }
  
  public boolean setup(ImageProcessor theBip, ImageProcessor theGip, Calibration theCal, int theMeasurements, int theSlice){
    
    if (isBinary(theBip) == false){
       IJ.showMessage("Input processor must be binary");
       return false;
    }
    
    wand = new Wand(theBip);
    wand.setAllPoints(true);
    
    slice = theSlice;
    regionID = 0;
    blobs.clear();
    bip = theBip;
    gip = theGip;
    hasGray = (gip != null);
    
    
    if (theCal != null) {
      gcal = theCal;
    } else {
      gcal = new Calibration();
    }
    
    measureFlag = theMeasurements;
    
    h = bip.getHeight();
    w = bip.getWidth();
    pw = w + 2;
    ph = h + 2;
    
    return true;    
  }  
    
  public boolean process(){
    
//    long startTime = System.currentTimeMillis();

    if (isBinary(bip) == false){
      IJ.log("First argument must be binary");
      return false;
    } 
    //if (verbose) {IJ.log("First argument is binary");}
       
    if (constructBlankLabel() == false){
      IJ.log("Error constructing the mask and label arrays");
      return false;
    }
    
    //if (verbose) {IJ.log("Masks created");}
    
    if (labelAndContour() == false){
      IJ.log("Error crafting the label and/or contours");
      return false;
    }
    
    //blobs.setVolumeMethod(volMethod);
//    long endTime = System.currentTimeMillis();

//    if (verbose) {
//      IJ.log("Labels, contours and blob measurements complete");
//      String rate = IJ.d2s((w*h) * 1000.0 /(endTime-startTime));
//      IJ.log("Elapsed time = " +IJ.d2s( ((endTime-startTime) /1000.0)) + "s");
//      IJ.log(rate + " pixels per second");
//      }
    return true;
  } //process 
 
/**
  * count the number of objects
  * @return the number of blobs - conveniently filters out placeholder empty blobs
  */
  public int countBlobs(){
    return blobs.countBlobs();
  }
  
/**
  * Constructs the labelMask and copies pixels from the 
  * input binary mask to origMask
  * @return returns true upon succesful completion
  */  
  private boolean constructBlankLabel(){
    
    labelMask = new int[ph][pw];
    origMask = new byte[ph][pw];    
    
    for (int y = 0; y < h ; y++){
      for (int x = 0; x < w; x++){
        if (bip.get(x,y) == 0){
          origMask[y+1][x+1] = BACKGROUND;
        } else {
          origMask[y+1][x+1] = FOREGROUND;
        }
      }//x loop
    }// y loop
    return true;
  }  
   
/**
  * This method handles the scanning, labeling and contouring
  * @return Returns true upon successful completion
  */
  
  private boolean labelAndContour(){
    BHSingleBlob singleBlob;
    int label = 0;
    if (includeBackground){
      singleBlob = new BHSingleBlob(label, gcal, measureFlag);
      singleBlob.setSlice(slice);
      blobs.add(singleBlob);
    }

    IJ.showStatus("Scanning for blobs");
    for (int y = 1; y < (ph-1); y++){
      IJ.showProgress(y, ph-1);
      label = 0; //at each row reset the label
     
      for (int x = 1; x < (pw-1); x++){

        if (origMask[y][x] == FOREGROUND) {
          //encounter a foreground pixel?
          if (label != 0) {
            // assign the label to the labelMask
            // this is a continuation of a scan
            // but only mark if labelMask is zero
            assignLabel(x,y, label); 
          } else {  //label != 0
            label = labelMask[y][x];
            if (label == 0){ //new outside contour
              regionID++;
              label = regionID;
              singleBlob = new BHSingleBlob(label, gcal, measureFlag);
              singleBlob.setSlice(slice);
              blobs.add(singleBlob);
              assignLabel(x,y,label);
              blobs.addContour(traceOuterContour(x, y, label));
              labelMask[y][x] = label; //assignLabel(x,y,label);
              blobs.addContour(doOuterWand(x,y,label));
            } // label == 0
          }//
          
        } else { // origMask == background
          if (label != 0){
            if (labelMask[y][x] == 0){
              blobs.addContour(traceInnerContour(x-1,y, label));
              //inner crack boundary doesn't work if there is a blob within
              //a blob so we'llhave to figure out how to get the inner wanding
              //to work without scanning across the object
              blobs.addContour(doInnerWand(x-1, y, label));
            }
            label = 0;
          } //label != 0
          if (includeBackground){assignLabel(x,y, BACKGROUND);}
      
        }// pixel is background or foreground?
        
      }// x-loop
      
    } // y-loop
    IJ.showStatus(" ");
    IJ.showProgress(2.0);
    
    blobs.translate(-1,-1);
    
    blobs.setSideFlag(0,0,w, h);
    return true;
  }
   
   
/**
  * This method follows the boundary, labeling the boundary pixels (required for CCL), but
  * the boundary contour extracted is the crack between pixels
  * @param xS The starting x coordinate
  * @param yS The starting y coordinate
  * @param label The label for the contour
  * @param c  The BHPolygon object
  * @return Returns the same contour object
  */
   private BHPolygon doWand(int startX, int startY, int label, BHPolygon c){
      int dx = 0, dy = 0;
      boolean ok = false;
      final int startDirection;
      
      
      if (inside(startX,startY))  {    // inside at left, outside right
           startDirection = 1;         // starting in direction 1 = up
      } else {
         startDirection = 3;         // starting in direction 3 = down
         startY++;             // continue after the boundary that has direction 3
      }
     

      int x = startX;
      int y = startY;
      int direction = startDirection;
      do {
          int newDirection;
          newDirection = direction + 1;
          do {
          if (inside(x, y, newDirection)) {break;}
              newDirection--;
          } while (newDirection >= direction);
          
          c.add(x+dx,y+dy);  //we have to add one because we are working on the 
                    
          switch (newDirection & 3) { // '& 3' is remainder modulo 4
              case 0: x++; break;
              case 1: y--; break;
              case 2: x--; break;
              case 3: y++; break;
          }
          direction = newDirection;
      } while (x!=startX || y!=startY || (direction&3)!=startDirection);
      if (startX!=x) {           // if the start point = end point is a corner: add to list
                c.add(x+dx, y+dy);
      }
      return c;
   }
/**
  * Trace a contour from the provided starting point extracting the boundary
  * locations as pixels (upper left)
  * @param xS The starting x coordinate
  * @param yS The starting y coordinate
  * @param label The label for the contour
  * @param c  The BHPolygon object
  * @return Returns the same contour object
  */
  private BHPolygon traceContour(int xS, int yS, int label, BHPolygon c){
    int dStart = (c.isOuter() ? 0 : 1);// outer = 0, inner = 1;
    int dNext, dSearch;
    BHPoint point = new BHPoint(xS, yS,label );
    
    int xT, yT; //successor
    int xP, yP; //previous
    int xC, yC; //current
    boolean done;
    
    dNext = findNextPoint(point, dStart);
    c.add(point);
    xP = xS; yP = yS;   //  Previous <- Start
    xC = xT = point.x;  //  Current <- Succesor <- point
    yC = yT = point.y;
    done = ((xS==xT) && (yS == yT));
    //if (done) { //isolated point?
    //  assignLabel(xS, yS, label);
    //} else {
      while (!done) {
        assignLabel(xC, yC, label);
        dSearch = (dNext + 6) % 8;
        point = new BHPoint(xC, yC, label);
        dNext = findNextPoint(point, dSearch);
        xP = xC; yP = yC; // Previous <- Current
        xC = point.x; yC = point.y; // Current <- Next
        
        done = ((xP==xS) && (yP== yS) && (xC== xT) && (yC==yT));
        if (!done) {
          assignLabel(xC, yC, label);
          c.add(point);
        }
        
      } //endwhile
    //} //isolated point?
        
    return c;
  }//traceContour
    
  
    // check pixel at (x,y), whether it is inside traced area
    private boolean inside(int x, int y) {
        return origMask[y][x] != BACKGROUND;
    }

    // check pixel in a given direction from vertex (x,y)
    private boolean inside(int x, int y, int direction) {
        switch(direction & 3) {         // '& 3' is remainder modulo 4
            case 0: return inside(x, y);
            case 1: return inside(x, y-1);
            case 2: return inside(x-1, y-1);
            case 3: return inside(x-1, y);
        }
        return false; //will never occur, needed for the compiler
    }

  
/**
  * Assign the label to the labelMask and adds the pixel to the 
  * blob indicated by label. If the pixel is already assigned then it 
  * does nothing.
  * @param x the pixel's x coordiante
  * @param y the pixels'y y coordinate
  * @param label the label to assign to this pixel
  * @return Returns true if the pixel is assigned and false otherwise
  */
  private boolean assignLabel(int x, int y, int label){
   boolean ok = false;
   if (labelMask[y][x] == BACKGROUND){
     labelMask[y][x] = label;
     double val = (hasGray) ? gip.getPixelValue(x-1,y-1) : 0.0;
     blobs.addPoint(label, x, y, val);
     ok = true;
   } else {
     ok = false;
   }
   return ok;
  }
  
 
/**
  * Trace an inner contour
  * @return the BHPolygon traced
  */
  private BHPolygon traceInnerContour(int x, int y, int label){
    BHPolygon c = new BHPolygon(BHPolygon.INNER, label, gcal);
    traceContour(x,y, label,c);
    return c;    
  }
  
/**
  * Trace an inner crack perimeter - note that the specified start location
  * is not going to be the same as the first location in the output polygon<br>
  * Inner crack contours must start with the background to the left
  * so the doWand method will scan across the inner hole from the start 
  * location specified until it reaches the other side of the hole.
  * 
  * @param x the starting pixel location
  * @param y the starting pixel location
  * @param label the label to assign to this contour
  * @return the BHPolygon traced
  */ 
  private BHPolygon doInnerWand(int x, int y, int label){
    BHPolygon c = new BHPolygon(BHPolygon.INNER + BHPolygon.CRACK, label, gcal);
    doWand(x+1,y, label, c);
    return c; 
  }
  
/**
  * Trace an outer contour
  * @param x the starting pixel location
  * @param y the starting pixel location
  * @param label the label to assign to this contour
  * @return the BHPolygon traced
  */
  private BHPolygon traceOuterContour(int x, int y, int label){
    BHPolygon c = new BHPolygon(BHPolygon.OUTER, label, gcal);
    traceContour(x,y, label, c);
    return c;
  }
  
/**
  * Trace an outer crack perimeter 
  * @param x the starting pixel location
  * @param y the starting pixel location
  * @param label the label to assign to this contour
  * @return the BHPolygon traced
  */ 
  private BHPolygon doOuterWand(int x, int y, int label){
    BHPolygon c = new BHPolygon(BHPolygon.OUTER + BHPolygon.CRACK, label, gcal);
    doWand(x,y, label, c);
    return c; 
  }
/**
  * Find the next contour point - alters point to next in line on contour
  * and returns next direction.
  * @param point The point under consideration
  * @param dir The initial search direction
  * @param the direction to turn next
  */
  private int findNextPoint(BHPoint point, int dir){
    int x;
    int y;
    for (int i = 0; i< 7; i++){
      x = point.x + delta[dir][0];
      y = point.y + delta[dir][1];
      if (origMask[y][x] == BACKGROUND) {
        labelMask[y][x] = VISITED;
        dir = (dir + 1) % 8; 
      } else {
         point.x = x;
         point.y = y;
         break;
      } //background or foreground?
    }//i loop
    return dir;
  }
///**
//  * Converts xy to index
//  * 
//  */
//  private int xyToIndex(int x, int y){
//    return y*w + x ;
//  }
  
///**
//  * Converts xy to pindex (padded index)
//  */
//  private int xyToPindex (int x, int y){
//    return (y+1)*pw + (x+1);
//  }
//  
/**
  * Converts index to xy
  */
  private int[] indexToxy(int index){
    int[] xy = new int[2]; 
    xy[1] = index/w;
    xy[0] = index % w;
    return xy;
  }
 
/**
  * Converts index to pindex
  */
  public int indexToPindex(int index){
    
    int[] xy = indexToxy(index);
    return (xy[1]+1)*pw + (xy[0]+1);
  }

/**
  * Converts pindex to pxy
  */
  public int[] pindexToPxy(int pindex){
    int[] pxy = new int[2];
    pxy[1] = pindex/pw;
    pxy[0] = pindex % pw;
    return pxy; 
  }
  
/**
  * Converts pindex to index
  */
  public int pindexToIndex(int pindex){
    int[] pxy = pindexToPxy(pindex);
    return (pxy[1]-1) * w + (pxy[0]-1);
  } 
       
  /*
   * Tests the binary processor to make sure it really is bilevel (or binary)
   * @param the ImageProcessor to test
   * @return returns true if the processor is bilevel
   */
  private boolean isBinary(ImageProcessor theIp){
    if ( theIp instanceof BinaryProcessor) { return true;}
    int[] hist;
    
    try { 
      hist = theIp.getHistogram(); 
    } catch (NullPointerException e) {
      return false;
    }
    return (hist[0] + hist[255]) == (theIp.getWidth()* theIp.getHeight());
  }
    
/**
  * Test that the image is Gray8, 16 or 32 bit
  */
  private boolean isGray(ImageProcessor theIp){
    boolean isG = false;
    if (theIp instanceof ByteProcessor) {
      isG = true;
    } else if (theIp instanceof ShortProcessor) {
      isG = true;
    } else if (theIp instanceof FloatProcessor) {
      isG = true;
    } else {
      
    }
         
    return isG;
  }     
  
 
///**
//  * This returns a int[nblob][nblob] table of 
//  * connectivity based upon the distance between outer
//  * contour(s) of a blob with its neighbors in the group.
//  * Associations are made to the lowest numbered label in each
//  * connectivity neighborhood.  Use the compressMergeTable()
//  * method to reduce the complexity of the table.
//  * Here's an example...  Label 1 is belongs to a group
//  * of 8 blobs (see first element of first row).  This
//  * group includes labels 1, 3, 6, 9, 10, 11, 13, and 15.
//  * Label 7 is all by itself.  And label 12 belongs to 
//  * a group of three (12, 14 and 16).
//  * Label
//  *  1	8	0	3	0	0	6	0	0	9	10	11	0	13	0	15	0
//  *  2	0	4	0	4	5	0	0	8	0	0	0	0	0	0	0	0
//  *  3	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  6	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  7	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0
//  *  8	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  9	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  10	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  11	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  12	0	0	0	0	0	0	0	0	0	0	0	3	0	14	0	16
//  *  13	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  14	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  15	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  *  16	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
//  * @param mergeIfCloserThan a distance threshold - blobs
//  *   closer (exclusive) than this value are merged on the table
//  *   but not in reality - use the BHBlob.add(otherBlob) functionality
//  *   for that.
//  **/
//  public int[][] computeMergeTable(double mergeIfCloserThan){
//    
//    if (blobs.count() == 0){
//      return null;
//    }
//    double[][] dist = blobs.computeMinimumDistanceTable();
//    int[] n = new int[dist.length];
//    int[][] idx = new int[dist.length][dist[0].length];
//    int[] labels = blobs.getLabels();//this is dist.length long
//    int[] lut = blobs.getLabelsLUT();
//    // this is a pointer of label associations
//    int[] ptr = new int[idx.length]; 
//       
//    //first make a look-up table based upon labels
//    //rather than distance
//    for (int y = 0; y < idx.length; y++){
//      ptr[y] = labels[y];
//      idx[y][0] = labels[y];
//      for (int x = y+1; x < idx[0].length; x++){ //was from y+1
//        if ((dist[y][x] > 0.0) && (dist[y][x] < mergeIfCloserThan)){
//          idx[y][x] = labels[x];
//          idx[x][y] = labels[x];//symmetric
//        }
//      }
//    }
//    
//    //make the follwoing true to spit out the raw connectivity table
//    if (true) {
//       IJ.log("idx table in bh.computeMergeTable() BEFORE");
//       for (int y = 0; y < idx.length; y++){
//         String s = "";
//         for (int x = 0; x < idx[0].length; x++){
//         s = s + idx[y][x] + " ";
//         }
//         IJ.log(s);
//      }
//    }
//  
//    // idx holds a list of indices?
//    for (int y = 0; y < idx.length; y++){
//      for (int x = (y+1); x < idx[y].length; x++){
//        if (idx[y][x] > 0){ 
//          int lutidx = lut[idx[y][x]];
//          if (idx[lutidx][0] == idx[y][x]){
//            idx[lutidx][0] = idx[y][0];
//            ptr[lutidx] = idx[y][0];
//          } else {
//            idx[y][0] = idx[lutidx][0];
//            ptr[y] = idx[lutidx][0];
//          }
//        } //if
//      } //x-loop
//     n[lut[idx[y][0]]]++; //increment the bin count
//    }//y-loop
//    
//    
//    //populate the return framework
//    int[][] ret = new int[idx.length][idx.length];
//    for (int y = 0; y< idx.length; y++){
//     ret[lut[idx[y][0]]][lut[labels[y]]] = labels[y];
//    }
//    
//    //count the number of connected blobs
//    int connCount = 0;
//    for (int y = 0; y<idx.length; y++){
//      if (n[y] > 0) { connCount++;}
//    }
//
//
//   if (true){
//      IJ.log("idx table after initial pass");
//      for (int i=0; i < idx.length; i++){
//         String s = "";
//         for (int j=0; j < idx[0].length; j++){
//            s = s + idx[i][j] + " ";
//         }
//         IJ.log(s);
//      }
//   }
//    //now reduce the return framework to a smaller ragged array
//    //with just one row per connected components
//    //each row will have just the labels of the blobs belonging to 
//    //this component
//    int[][] ret2 = new int[connCount][1];
//    int nconn = 0;
//    for (int y = 0 ; y < ret.length; y++){
//      if (n[y] > 0){
//        ret2[nconn] = new int[n[y]];
//        int nx = 0;
//        ret2[nconn][nx] = ret[y][0];
//        for (int x = 0; x< ret[y].length; x++){
//           if (ret[y][x] != 0){
//             ret2[nconn][nx] = ret[y][x];
//             nx++;
//           }
//        }//x-loop 
//        nconn++;  
//      } //if n[y] > 0
//    }//y-loop    
//    
////    IJ.log("computeMergeTable() complete");
//    //assign the count to the first element of each row
//    //for (int i = 0; i < n.length; i++){idx[i][i]=n[i];}
//    return ret2;
//  }
  


  
//  /** Compresses the merge table to a ragged array
//    * Using the example from computeMergeTable the result
//    * would look like...
//    * 1  3  6  9  10  11  13  15
//    * 2  4  5  8
//    * 7
//    * 3  12  14  16 
//    * @param table merge array from computeMergeTable()
//    */
//  public int[][] compressMergeTable(int[][] table){
//    
//    int[] labels = blobs.getLabels();
//    int[] lut = blobs.getLabelsLUT();
//    //the counts per row
//    int[] n = new int[table.length]; 
//    for (int y = 0; y< table.length; y++){
//      n[lut[table[y][0]]]++;
//    }
//    
//    
//  
//    int[][] ret = new int[table.length][table.length];
//    for (int y = 0; y< table.length; y++){
//     ret[lut[table[y][0]]][lut[labels[y]]] = labels[y];
//    }
//    
//    int connCount = 0;
//    for (int y = 0; y<table.length; y++){
//      if (n[y] > 0) { connCount++;}
//      String s = "n[" + labels[y] + "] = " + n[y] + " ->";   
//      for (int x = 0; x<table[y].length;x++){
//        s = s + ret[y][x] + " ";
//      }
//      IJ.log(s);
//    }
//    IJ.log("connCount=" + connCount);
//    
//    int[][] ret2 = new int[connCount][1];
//    int nconn = 0;
//    for (int y = 0 ; y < ret.length; y++){
//      if (n[y] > 0){
//        ret2[nconn] = new int[n[y]];
//        int nx = 0;
//        ret2[nconn][nx] = ret[y][0];
//        for (int x = 0; x< ret[y].length; x++){
//           if (ret[y][x] != 0){
//             ret2[nconn][nx] = ret[y][x];
//             nx++;
//           }
//        }//x-loop 
//        nconn++;  
//      } //if n[y] > 0
//    }//y-loop
//      
//    return ret2;  
//  }
// 
// 
// 
///**
//  * Returns the compressed merge table given the threshold
//  * distance.
//  * @param mergeIfCloserThan  The exclsuive distance that blobs must
//  * be within to be merged
//  */
//  public int[][] compressMergeTable(double mergeIfCloserThan){
//      int[][] table = computeMergeTable(mergeIfCloserThan);
//      return compressMergeTable(table);
//  }
    
    

  public ImageProcessor createContourMask(){
   return createContourMask(BHPolygon.OUTER);
  }
  
  public ImageProcessor createContourMask(int typeOfContour){
      ImageProcessor contourmask = new ShortProcessor(w,h);
//      blobs.labelContoursProcessor(contourmask, typeOf);
//      return contourmask;
      return createContourMask(contourmask, BHPolygon.OUTER);
     
  } 
  
  public ImageProcessor createContourMask(ImageProcessor cmask, int typeOfContour){
      blobs.labelContoursProcessor(cmask, typeOfContour);
      return cmask;     
  }
  
  public ImageProcessor createCountMask(){ 
   ImageProcessor countmask = new ShortProcessor(w,h);
   blobs.labelPixelsProcessor(countmask);
   return countmask; 
  }
} //BlobHandler class