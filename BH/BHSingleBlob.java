import java.util.Vector;
import ij.process.*;
import ij.*;
import ij.measure.*;
import ij.gui.*;
import ij.plugin.filter.*;
/**
  * Manages the blob statistics - points, contours, and binary and gray stats -
  * for a single contiguous group of pixels.  If I had to name these again I 
  * would have BlobBase, Bloblet and Blob where this is a Bloblet. 
  * <p>
  * This is where the perimeter contours (outer and inner) are stored as well
  * as the convex hull.  The convex hull is not computed until requested, and then 
  * it is stored so it doesn't have to be recomputed. 
  * <p> 
  * Provisions have been made for the storage of the perimeter that follows the
  * crack between pixels, like ImageJ's Wand produces.  At this time the 
  * crack contours are not computed.
  * <p>
  * Contact: Mike Sieracki, msieracki@bigelow.org <br>
  * Author:  Ben Tupper, btupper@bigelow.org <br>
  * History: Developed over the course of 2009 and 2010 using ImageJ 1.43 and Mac OSX 10.5.8 <br>
  *   at Bigelow Laboratory for Ocean Science, www.bigelow.org<p><p>
  * 2010-02-27 BTT now uses Calibration for coordinate conversions
  * 2010-04-07 BTT added "crack" boundary stuff
  */


public class BHSingleBlob extends BHBlob {
    /** The pixles as point */
    public BHPointGroup points;//the pixels as points
    /** The convex hull */
    public BHPolygon chull;
    /** The polygon from CCL describing the outer contour as pixel locations (upper left) */
    public Vector<BHPolygon> opixel; 
    /** The polygon(s) from CCL describing the inside contours(s) as pixel locations (upper left) */
    public Vector<BHPolygon> ipixel;
    /** For the outer crack contour*/
    public Vector<BHPolygon> ocrack; 
    /** For the inner crack contour(s) - not used */
    public Vector<BHPolygon> icrack;
        
/**
  * Initialize
  */ 
  public BHSingleBlob(){
    initialize(new Calibration());
  }
  public BHSingleBlob(Calibration cal, int measurements){
   measureFlag = measurements;
   initialize(cal); 
  } 
  public BHSingleBlob(int theLabel){
    label = theLabel;
    initialize(new Calibration());
  }
  public BHSingleBlob(int theLabel, int measurements){
    label = theLabel;
    measureFlag = measurements;    
    initialize(new Calibration());
  }
  public BHSingleBlob(int theLabel, Calibration cal){
    label = theLabel;
    initialize(cal);
  }
  public BHSingleBlob(int theLabel, Calibration cal, int measurements){
    label = theLabel;
    measureFlag = measurements;    
    initialize(cal);
  }
/**
  * Initialize - sets everything to 0 and empties the containers
  * It does NOT reset the label.
  */
  protected void initialize(Calibration cal){
    
    super.initialize(cal);

    if (points != null) {
      points.clear();
    } else {
      points = new BHPointGroup();
    }
    
    if (opixel != null){
      opixel.clear();
    } else {
      opixel = new Vector<BHPolygon>();
    }
    if (ocrack != null){
      ocrack.clear();
    } else {
      ocrack = new Vector<BHPolygon>();
    }
        
    if (ipixel != null) {
      ipixel.clear();
    } else {
      ipixel = new Vector<BHPolygon>();
    }
    if (icrack != null) {
      icrack.clear();
    } else {
      icrack = new Vector<BHPolygon>();
    }
  }  


/**
  * Returns the number of points (pixels)
  */
  public int npoints(){
    return points.npoints;
  }
/**
  * Computes the binary and grayscale
  */
  public void computeStats(){

    
    //the features get reset
    resetFeaturesToDefaults();
    BHPoint point;
    
    //this is the area covered by the pixels
    area = (double) points.npoints * cal.pixelWidth * cal.pixelHeight;
    //area = (double) points.npoints * (xscale[0] * yscale[0]);
    
    //outerArea, enclosed by outer contour
    outerArea = computeContourArea(BHPolygon.OUTER + BHPolygon.CRACK);
    //innerArea, this is the void area (may be 0)
    innerArea = computeContourArea(BHPolygon.INNER + BHPolygon.CRACK);
    //voidFraction, this is the ratio of inner/outer
    voidFraction = innerArea/outerArea;
    //nEuler
    nEuler = opixel.size() - ipixel.size();
    
    super.computeMoments(points.getArray());
   
    xc = cal.getX(xsum/points.size()+ 0.5) ;
    yc = cal.getY(ysum/points.size()+ 0.5) ;
      
    gmass = gsum1; 
    gmean = gmass/points.npoints;
    double mean2 = gmean*gmean;
    gvar = gsum2/points.npoints - mean2;
    double sDeviation = Math.sqrt(gvar);      
    gskew = ((gsum3 - 3.0*gmean*gsum2)/points.npoints  + 2.0*gmean*mean2)/(gvar*sDeviation);
    gkurt = (((gsum4 - 4.0*gmean*gsum3 + 6.0*mean2*gsum2)/points.npoints - 
      3.0*mean2*mean2)/(gvar*gvar)-3.0);
    xcm = cal.getX(gxsum/gsum1 + 0.5);
    ycm = cal.getY(gysum/gsum1 + 0.5);
    
    //get the convex hull of the first item
    if ((opixel != null) && (opixel.size() > 0)) {
      chull = getConvexHull();
    }
    
    //perim = getOuterPerimeterLength();
    perim = getBoundaryLength(BHPolygon.OUTER + BHPolygon.CRACK);

    //}
    
    // these match IPLab with very good precision
    m11 = xysum - (xsum*ysum)/ points.npoints;
    m20 = x2sum - (xsum*xsum) / points.npoints;
    m02 = y2sum - (ysum*ysum) / points.npoints;
    mtheta = getAngle(m20, m11, m02);
    
    if (m02 != m20) {
      morient = 0.5 * Math.atan(2*m11/(m20-m02));
    }
      
    computeEllipse(true);
    //if ((measureFlag & ELLIPSE) != 0){
      //see Imaging Book 11.4 Properties of Binary Regions 
      //double a1 = u20 + u02 + Math.sqrt((u20-u02)*(u20-u02) + 4.0 * u11*u11);
      //double a2 = u20 + u02 - Math.sqrt((u20-u02)*(u20-u02) + 4.0 * u11*u11);    
      //mra = cal.getX(Math.sqrt(2.0*a1/points.size()));
      //mrb = cal.getY(Math.sqrt(2.0*a2/points.size()));
      //mecc = a1/a2;
    //}
    
    //if ((measureFlag & SHAPE) != 0){
      feret = computeFeret(ocrack.get(0).getXY(false));
      
          //many of these can be found in any image processing text.
      // e.g. see tbale 10.2 in Russ 3 ed.
      if (npoints() > 1) {
        circ = (perim == 0.0) ? 
          0.0 : 4.0 * Math.PI * area/Math.pow(perim,2.0);
      } else {
        circ = 1.0;
      }
      if (circ > 1.0) {circ = 1.0;}
      if (chull != null ){
        solid = computeArea(BHPolygon.OUTER + BHPolygon.CRACK)/chull.computeArea();
      }
    
      double[] box = getBoundingBox(false);
      compact = Math.sqrt((4.0 / Math.PI) * area) / feret[FERETMAX];
      extent = area/(box[2]*box[3]);  
      //aspect = feret[FERETMAX]/feret[FERETMIN];
      aspect = mra/mrb;
      round = (4.0 * area) / (Math.PI*Math.pow(feret[FERETMAX], 2.0));
     }//computeStats
 
 

  

/**
  * Computes the volume estimate based upon the switch volMethod control
  * Currently then CCL style perimeters are required for the SVW_Volume method.
  * the volume must be calculated for each particle added to this blob
  */ 
  public double computeVolume(){
// [[public static int VOL_SKEDM = 3, VOL_SVW = 2, VOL_ABD = 1, VOL_NONE = 0, VOL_AUTO = -1;
// public int volCutoff = 0.4; //for auto  <cutoff use VOL_SKEDM >= cutoff  use VOL_SVW
    
      if (volMethod == VOL_SVW) {
        if ((usePerimeter & USECCLPERIMETERS) != 0) {
          mvol = computeVol_SVW();
        } else {
          mvol = 0.0;
        }
      } else if (volMethod == VOL_ABD){
        mvol = computeVol_ABD(); 
      } else if (volMethod == VOL_SKEDM) {
        mvol = computeVol_SKEDM();
      } else if (volMethod == VOL_AUTO){
        if (solid < volCutoff) {
          mvol = computeVol_SKEDM(); //use SKEDM
        } else {
          if ((usePerimeter & USECCLPERIMETERS) != 0) {
            mvol = computeVol_SVW();
          } else {
            mvol = 0.0;
          }//perimeters are CCL?
        }//use SVW
      } else {
        mvol = 0.0; //VOL_NONE or unrecognized
      } 
      
     return(mvol);
  }//computeVol
 
 
/**
  * Computes the volume using the product of the skeletonization (SK) and the 
  * Euclidean  Distance Map (EDM). This is proabaly the slowest of all the methods
  * but timing is untested.
  */
  public double computeVol_SKEDM(){
   return super.computeVol_SKEDM(this);
  } //compueVol_SKEDM
  
/**
  * Computes the volume by transforming the area to an 
  * area based radius.  Use that radius is used in
  * vol = 4/3*pi*r^3
  */
  public double computeVol_ABD(){
    mvol = 0.0;
    double r = Math.sqrt(area/Math.PI);
    mvol = 4.0/3.0*Math.PI*r*r*r; 
    return mvol;
  }
  
/**
  * Computes the volume via Sieracki, Viles and Webb (1989). 
  */
  public double computeVol_SVW(){
     boolean debug = false;
     mvol = 0.0;
     //get biovolume and ferets
      if ((opixel == null) || (opixel.size() == 0)) {return mvol;}
       if (debug) {IJ.log("ComputeVol for " +label);}
        //get the *first* outer contour
        BHPolygon c = opixel.get(0);
        //get the rotated points - rotation can cause holes
        double xcPixel = xc/cal.pixelWidth + cal.xOrigin - 0.5;
        double ycPixel = yc/cal.pixelHeight + cal.yOrigin - 0.5;
        double[][] xy = c.getRotatedRaw(xcPixel, ycPixel, getPhi(mtheta));
        if (debug) {for (int i = 0; i < xy.length; i++){IJ.log(xy[i][0]+ ", "+xy[i][1]);}}
        int counter = 0;
        double[] x;
        double[] y;
        if (xy.length > 1) {
          x = new double[xy.length*2]; //make new array big enough
          y = new double[xy.length*2]; //to hold the "filled" contour
          double curx = 0.0;
          double nextx = 0.0;
          int curi = 0;
          int nexti = 0;
          for (int i = 0; i < xy.length; i++){
            //be careful at the end - this wraps the last point
            // back to the next to first using modulo 
            curi = i;   
            nexti = (i+1) % (xy.length); //wrap if needed
            curx = rint(xy[curi][0]);   //round the "cur/next" values
            nextx = rint(xy[nexti][0]);
            dx = nextx - curx;
            dy = xy[nexti][1] - xy[curi][1];
            if (Math.abs(dx) > 1.0) {
              x[counter] = xy[curi][0]; //average of the pair fills in
              x[counter+1] = (dx * 0.5) + xy[curi][0];
              //x[count+1] = (dx * 0.5) + curx;
              y[counter] = xy[curi][1];
              y[counter+1] = (dy * 0.5) + xy[curi][1];
              counter = counter + 2;
              if (debug) {
                IJ.log("gap=" + dx + " at [x1,x2]=" + curx + ", "+ nextx);
                IJ.log("        was = "+ xy[curi][0] + ", " + xy[nexti][0]);
                IJ.log("        curi=" + curi + ", nexti=" + nexti);
                IJ.log("                   [y1,y2]=" +xy[curi][1] + ", "+ xy[nexti][1]); 
              }
            } else {
              if (debug) {IJ.log("no gap = "+ dx + "at [x1,x2]=" + curx + ", "+ nextx);}
              x[counter] = xy[curi][0]; //no gap found
              y[counter] = xy[curi][1];
              counter = counter + 1;            
            } //is there a gap ?
          }//i-loop to fill in the gaps created by the transform
        } else {
          //only one point found
          counter = 1;
          x = new double[1];
          y = new double[1];
          x[0] = xy[0][0];
          y[0] = xy[0][1];
        } //xy.length > 1
     
  
     // find the ranges of x and y in the rotated and filled
     // these are used to determine the ferets (and bins)
     // note that x and y have more elements than valid points
     // thus we find the range only within the 0 to (count-1) range
       double[] mmx = range(x, counter);
       double[] mmy = range(y, counter);
           
       //rounding each matches the result of pmbiovol.c
       vferet = rint(mmy[1]) - rint(mmy[0]) + 1.0;
       hferet = rint(mmx[1]) - rint(mmx[0]) + 1.0;
  
       // compute bins and the min/max y for each bin
       int bins = (int) hferet;
       double[] ymax = new double[bins]; // the upper and lower y values 
       double[] ymin = new double[bins]; // in each bin
       double x0 = rint(mmx[0]); //or should it be an int?
       //intialize the y ranges
       double bignum = Double.MAX_VALUE;
       for (int i = 0; i<bins; i++){
        ymax[i] = 0.0-bignum;
        ymin[i] = bignum;
       }
       
      int xtmp;
      double ytmp;
      //find the min/max in y of each bin
      for (int i = 0; i < counter; i++){
        xtmp = (int) rint(x[i]-x0); // this is the "index" for the y values
        ytmp = y[i];
        //the following is a kluge for a rare issue with
        //out-of-bound indexing on small objects
        if (xtmp < 0) {xtmp = 0;}
        if (xtmp >= ymin.length) {xtmp = ymin.length-1;}
        
        if (debug) {
          IJ.log("[i,xtmp,ytmp]="+ i + ", " + xtmp + ", " + ytmp + ", y[i]=" + y[i]);
          IJ.log("  before [ymin, ymax]=" + ymin[xtmp] + ", " + ymax[xtmp]);
        }
        
        if (ytmp < ymin[xtmp]) {ymin[xtmp] = ytmp;}
        if (ytmp > ymax[xtmp]) {ymax[xtmp] = ytmp;}
        if (debug) {
          IJ.log("  after [ymin, ymax]=" + ymin[xtmp] + ", " + ymax[xtmp]);
        }
      } 
      
     
      double r; 
      for (int i = 0; i < bins; i++){
        if (debug) {
          IJ.log("i, ymax, ymin " + i + ", " + ymax[i] + ", " + ymin[i] );
        }
        //radius is half the y difference for each increment plus one 
        //(because pixels have size 1) 
        r = (rint(ymax[i] - ymin[i]) + 1.0) * 0.5 * cal.pixelHeight; 
        // sum  of pi * r^2 * width
        mvol += r*r*cal.pixelWidth; 
        //mvol += r*r*xscale[0]; 
       } // i loop through the bins
      
       mvol = mvol  * Math.PI;   
       hferet = hferet*cal.pixelWidth;
       vferet = vferet*cal.pixelWidth;
       //hferet = hferet*xscale[0];
       //vferet = vferet*yscale[0];
     return mvol;      
  }//computeVol_SVW
  
///**
//  * Returns the aggregate perimeter for all outer contours associated with
//  * this blob.  This makes no sense if there is only one outer contour
//  */
//  public double getOuterPerimeterLength(boolean allContours){
//   double perim = ocrack.get(0).getLength();
//   if (allContours == true){
//    for (int i = 1; i < ocrack.size(); i++){
//      perim += (ocrack.get(i)).getLength();
//    }
//   }
//   return perim;
//  } 
//  
///**
//  * Returns the aggregate perimeter for all outer contours associated with
//  * this blob
//  */
//  public double getOuterPerimeterLength(){
//   return getOuterPerimeterLength(true);
//  }  
  
/**
  * Returns the boundary length for the specified boundary
  * @param which the integer specifying which boundary to measure (OUTER + CRACK, etc.)
  */
  public double getBoundaryLength(int which){
      BHPolygon[] poly = getContours(which);
      double p = 0.0;
      for (int i = 0; i < poly.length; i++){
         p += poly[i].getLength();
      }
      return p;
  }
  
/**
  * Add a point
  */
  public int addPoint(BHPoint point){
     //add the point
     points.add(point);
     //return the number of points
     return points.size();
  }
  
/**
  * Adds a point
  */  
  public int addPoint(int x, int y){
    return addPoint(new BHPoint(x,y));   
  }
  
  
/**
  * Add a contour - sorted internally by the contour.type
  * It does NOT check that the label matches the blob label
  */
  public void addContour(BHPolygon c){
    if (c.isOuter()) {
        //IJ.log("Adding outer contour to " + c.label);
       if (c.isCrack()) {
         ocrack.add(c);
       } else {
         opixel.add(c);
       }
    } else {
      if (c.isCrack()) {
        icrack.add(c);
      } else {
         ipixel.add(c);
      }
    }
  }//addContour
    
  
/**
  * Label the IP's pixels with the specified color
  */
  public void labelPixels(ImageProcessor ip, int col){
    //BHPoint[] p = getPointsAsArray();
    for (int i = 0; i < points.npoints; i++){
      ip.putPixel(points.points[i].x, points.points[i].y, col);
    }

  }//labelPixels
  
/**
  * Label the IP's pixels with the label value
  */  
  public void labelPixels(ImageProcessor ip){
    labelPixels(ip, label); 
  }//labelPixels
  

/**
  * Returns the upper left most pixel
  */ 
  public BHPoint getStart(){
     return points.getStart();
  }
/**
  * Labels the contour with the specified color
  */
  public void labelContours(ImageProcessor ip, int col, int which){
  
  
    BHPolygon[] poly = getContours(which);
    if ((poly != null) && (poly.length > 0)){
      for (int i = 0; i < poly.length; i++){
       //BHPolygon c = opixel.get(i);
       poly[i].labelContour(ip, col);
      }// i loop
    }    
//    
//    if (which == BHPolygon.OUTER) {
//      if ((opixel != null) && (opixel.size() > 0)){
//        for (int i = 0; i < opixel.size(); i++){
//         BHPolygon c = opixel.get(i);
//         c.labelContour(ip, col);
//        }// i loop
//      }
//    } else {
//      if ((ipixel != null) && (ipixel.size() > 0)){
//        for (int i = 0; i < ipixel.size(); i++){
//         BHPolygon c = ipixel.get(i);
//         c.labelContour(ip, col);
//        }// i loop
//      } 
//      
//    }
    
    
  }//labelOuterContours

/**
  * labels the contex hull with the supplied color
  */
  public void labelConvexHull(ImageProcessor ip, int col){
    if ((chull != null) || (chull.npoints > 0)){
      for (int j = 0; j < chull.npoints; j++){
          ip.putPixel(chull.points[j].x, chull.points[j].y, col);
      } // j loop
    }
  }
    
/**
  * Returns the pixel values for the blob within the imageprocessor
  */
  public float[] getPixelValues(ImageProcessor ip){
    BHPoint[] p = getPointsAsArray();
    float[] values = new float[points.size()];
    for (int i = 0; i < points.size() ; i++){
      values[i] = ip.getPixelValue(p[i].x, p[i].y);
    }
    return values;
  } //getPixelValues

/**
  * Returns the points as an array
  */
   public BHPoint[] getPointsAsArray(){
   
      return points.getArray();
   }
   
/**
  * Moves the blob by dx, dy
  */
  public void translate(int dx, int dy){
    
    if (points.size() != 0) {
      BHPoint[] p = points.get();
      for (int i = 0; i<p.length; i++) {
        p[i].translate(dx,dy);
      }
    }
    
    BHPolygon[] poly;
      
    poly = getContours(BHPolygon.OUTER);
    if (poly != null){
      for (int i = 0; i<poly.length; i++) {
        poly[i].translate(dx,dy);
      }      
    }
    poly = getContours(BHPolygon.OUTER + BHPolygon.CRACK);
    if (poly != null){
      for (int i = 0; i<poly.length; i++) {
        poly[i].translate(dx,dy);
      }      
    }
    poly = getContours(BHPolygon.INNER);
    if (poly != null){
      for (int i = 0; i<poly.length; i++) {
        poly[i].translate(dx,dy);
      }      
    }    
    poly = getContours(BHPolygon.INNER + BHPolygon.CRACK);
    if (poly != null){
      for (int i = 0; i<poly.length; i++) {
        poly[i].translate(dx,dy);
      }      
    }    
       
//    if (ipixel.size() != 0){
//      BHPolygon[] ic = (BHPolygon[]) ipixel.toArray(new BHPolygon[ipixel.size()]);
//      for (int i = 0; i<ic.length; i++) {
//        ic[i].translate(dx,dy);
//      }      
//    }
//
//    if (icrack.size() != 0){
//      BHPolygon[] ic = (BHPolygon[]) ipixel.toArray(new BHPolygon[ipixel.size()]);
//      for (int i = 0; i<ic.length; i++) {
//        ic[i].translate(dx,dy);
//      }      
//    }
   
    //if the location changes it is worth recalculating
    computeStats();
  }
 



/**
  * Shows this blob's info on the specified row of ResultsTable
  */   
  public void showInfo(ResultsTable rt, int row){
    
       BHPoint pointOne = points.getStart();
        if ((measureFlag & BASIC) != 0){
          rt.setValue("Blob", row, (double) label);
          rt.setValue("Count", row, 1.0);  // a bit silly, but consistent
          rt.setValue("Area", row, area);
          rt.setValue("Perim", row, perim);
        }
       if ((measureFlag & SLICE) != 0){
  		   rt.setValue("Slice", row, (double) slice);
  		 }
  		  if ((measureFlag & LOCATION) != 0) {
    		  rt.setValue("BX", row, xrange[0]);
    		  rt.setValue("BY", row, yrange[0]);
    	 	  rt.setValue("BW", row, xrange[1]-xrange[0] + 1);
    		  rt.setValue("BH", row, yrange[1]-yrange[0] + 1);
  		    rt.setValue("XC", row, xc);
  		    rt.setValue("YC", row, yc);
    		  rt.setValue("XStart", row, pointOne.x);
    		  rt.setValue("YStart", row, pointOne.y);
    		  rt.setValue("Side", row, (double) side);
  		  }
  		  	  
  		 // rt.setValue("Hferet", row, hferet);
  		 // rt.setValue("Vferet", row, vferet);
  		 if ((measureFlag & SHAPE) != 0){
        rt.setValue("FeretMax", row, feret[0]);
        rt.setValue("FeretMin", row, feret[2]);
        rt.setValue("FeretMean", row, feret[5]);
        rt.setValue("FeretMed", row, feret[6]);
        rt.setValue("FeretAngle", row, Math.toDegrees(feret[1]));
        rt.setValue("FeretX", row, feret[3]);
        rt.setValue("FeretY", row, feret[4]);
  		  rt.setValue("Theta", row, Math.toDegrees(mtheta));
  		  rt.setValue("Orient", row, Math.toDegrees(morient));
  		  rt.setValue("Circ", row, circ);
  		  rt.setValue("Round", row, round);
  		  rt.setValue("Solid", row, solid);
  		  rt.setValue("Compact", row, compact);
  		  rt.setValue("Extent", row, extent);
  		  rt.setValue("Aspect", row, aspect);
  		  rt.setValue("VoidFraction", row, voidFraction);
  		  rt.setValue("NEuler", row, nEuler);
  		 }

  		 if ((measureFlag & ELLIPSE) != 0){
    		  rt.setValue("Semimajor", row, mra);
    		  rt.setValue("Semiminor", row, mrb);
  	   	  rt.setValue("Ecc", row, mecc);
  		 }
  		   		 
  		 if ((measureFlag & VOLUME) != 0) {
  		  rt.setValue("Vol", row, mvol);
  		  //the following is a temporary measure for the KN193 cruise
  		  //what this computes is the number of SingleBlobs with Solid < 0.4
  		  //and stores it in the results table with the name NTricky
  		  int ntricky = 0;
  		  if (solid < 0.4) {ntricky = 1;}
  		  rt.setValue("NTricky", row, ntricky);
  		  
  		 }
  		 
  		 if ((measureFlag & GRAY) != 0) {
  		    rt.setValue("GMin", row, grange[0]);
  		    rt.setValue("GMax", row, grange[1]);
  		    //rt.setValue("GMass", row, gmass);
  		    rt.setValue("GMean", row, gmean);
  		    rt.setValue("GVar", row, gvar);
  		    rt.setValue("GSkew", row, gskew);
  		    rt.setValue("GKurt", row, gkurt);
  		  
  		 }
  		 
  		 if ((measureFlag & MOMENTS) != 0) {
  		  for (int k = 0; k < nHu; k++){rt.setValue("Hu" + (k+1), row, mHu[k]);}
  		  for (int k = 0; k < nAMI; k++){ rt.setValue("AMI" + (k+1), row, mAMI[k]);}
  		 }

  		 
  }//showInfo (at a specified row)

/**
  * Shows this blob's info on the last row of ResultsTable
  */  
  public void showInfo(ResultsTable rt){
    showInfo(rt, rt.getCounter()); 
  } //showInfo at last row
  
  public void showInfo(){
    IJ.log("Blob = " + label);
    //IJ.log("  xscale=" + xscale[0] + "," + yscale[0]);
    IJ.log("  NPoints = " + points.npoints);
    IJ.log("  [xc, yc] = " + IJ.d2s(xc) + ", " + IJ.d2s(yc));
    //IJ.log("  [xs, ys] = " + points.points[0].x + ", " + points.points[0].y);
    IJ.log("  [xS, yS] = " + IJ.d2s(xrange[0]) + ", " + IJ.d2s(yrange[0]));
    IJ.log("  [wB, hB] = " + (xrange[1]-xrange[0]+1) + ", " + 
      (yrange[1]-yrange[0] + 1));
    IJ.log("  [hferet, vferet] = " + hferet + ", " + vferet);
    IJ.log("  orientation = " + IJ.d2s(mtheta));
    IJ.log("  mvol = " + IJ.d2s(mvol));
  }
  public void showDetailedInfo(){
    showInfo();
    points.printPoints(); 
  }
  
  
  /*
  * Returns the angle phi given the orinetation theta
  */
  public double getPhi(double theta){
      double phi = 0.0;
      if (theta > (Math.PI*0.5)) {
        phi = Math.PI - theta;
      } else {
        phi = 0.0 - theta;
      }
    return phi;
  }//getPhi
  
/*
 * Returns the range from 0 to count in the array v
 */
  public double[] range(double[] v, int count){
    double[] mm = {v[0], v[0]};
    for (int i = 1; i < count ; i++){
      if (v[i] < mm[0]) { mm[0] = v[i];}
      if (v[i] > mm[1]) { mm[1] = v[i];}
    }
    return mm;
  }//range
  
/*
 * Returns the range from 0 to v.length in the array v
 */
  public double[] range(double[] v){
    return range(v, v.length);
  } //range

/*
 * Returns the angle in radians given the m20, m11 and m02 binary 
 * moments.
 */  
  public double getAngle(double m20, double m11, double m02){
   double mtheta;
   if (m02 == m20) {
      if ((m20 * m11) > 0.0) {
         mtheta = Math.PI/4.0;
      } else if (m11 == 0.0) {
         mtheta = 0.0d;
      } else {
         mtheta = 0.0 - Math.PI/4.0;
      }
   } else {
     mtheta = Math.atan(2.0*m11/(m20-m02))/2.0;
   }
   
   //orient relative to x axis
   if (m02 > m20) {
      mtheta = Math.PI/2.0 - mtheta;
   } else if ((m02 < m20) && (m11 > 0.0)) {
      mtheta = Math.PI - mtheta;
   } else {
      mtheta = 0.0-mtheta;
   }
   
   return mtheta;
  }// getAngle


/**
   * Get the convex hull stored in a vector
   */
   public Vector<BHPolygon> getConvexHullVector(){
   
      Vector <BHPolygon> v = new Vector<BHPolygon>(); 
      v.add(getConvexHull());
      return v;
   }
   
/**
   * Get the convex hull stored in a vector
   */
   public Vector<BHPolygon> getConvexHullVector(int which){
   
      Vector <BHPolygon> v = new Vector<BHPolygon>(); 
      v.add(getConvexHull(which));
      return v;
   }
/**
  * Computes the convex hull of the OUTER + CRACK boundary
  */
  public BHPolygon getConvexHull(){
   return getConvexHull(BHPolygon.OUTER + BHPolygon.CRACK);
  }
/**
  * Computes the convex hull enveloping the boundary specified
  */
  public BHPolygon getConvexHull(int which){
   //does it already exist and match the specified one?
   if (chull != null) {
      if (chull.type == which){    
         return chull;
      }
   }
   //can it be made?
   BHPolygon[] poly = getContours(which);
   if ((poly == null) || (poly.length == 0)) { return null;} 
   //make it 
   chull = poly[0].getConvexHull();
   if (chull == null) {return null;}
   //transfer the calibration and label
   chull.label = label;
   chull.cal = cal;
   return chull;
  }// gets the 


/**
  * Computes the aggregate area contained within the inner/outer perimeter(s)
  * @param which The 0/1 flag to indicate which contours to use. (See BHPolygon.OUTER and BHPolygon.INNER)
  * @return The aggregate area of the contour(s)  
  */
  public double computeArea(int which){
   
   double polyarea = 0.0;
   BHPolygon[] poly = getContours(which);
   
   if (poly != null){
     for (int i = 0; i < poly.length; i++){
      polyarea += poly[i].computeArea(); 
     }
   } 
   return polyarea;
 }//computeArea
   
/**
  * Returns a contour of the equivalent ellipse
  */   
 public PolygonRoi getEllipseRoi(){
  return getEllipseRoi(this);
 }


   
/**
  * Computes the minimum distance between the outer contour(s) of this 
  * blobs with that provided. Assumes no overlap/intersection/enclosing
  * etc.  
  */
  public double computeMinimumDistanceTo(BHSingleBlob b){
    //does the other have contours?
    //if (b.countContours(BHPolygon.OUTER) == 0) { return Double.NaN;}
    //does this have contours?
    //if (countContours(BHPolygon.OUTER) == 0) { return Double.NaN;}
    
    double temp;
    double d = Double.MAX_VALUE;
    BHPolygon[] poly1 = getContours(BHPolygon.OUTER + BHPolygon.CRACK);
    BHPolygon[] poly2 = b.getContours(BHPolygon.OUTER + BHPolygon.CRACK);
    
    if ((poly1 == null) || (poly2 == null)){ return Double.NaN;}
    
    for (int i = 0; i < poly1.length; i++){
      for (int j = 0; j < poly2.length; j++){
         temp = poly1[i].computeMinimumDistanceTo(poly2[j]);
         if (temp < d) { d = temp;}
      } // j-loop
    } //i=loop

    return d;
  } //computeMinimumDistanceTo
  

/**
  * Gets the contour vector set specified by OUTER, INNER, CRACK, etc.
  */
  public Vector<BHPolygon> getContourVector(int which){
  
      Vector<BHPolygon> thisV = null;
      switch(which) {
         case BHPolygon.OUTER:  thisV = opixel; break;
         case BHPolygon.INNER:  thisV = ipixel; break;
         case (BHPolygon.OUTER + BHPolygon.CRACK): thisV = ocrack; break;
         case (BHPolygon.INNER + BHPolygon.CRACK): thisV = icrack; break;
         case (BHPolygon.CHULL): thisV = getConvexHullVector(); break;
         case (BHPolygon.CHULL + BHPolygon.CRACK): thisV = getConvexHullVector(which); break;
         default:  break;
      }
      
      return thisV;
  
  }
/**
  * Counts the number of contours including siblings
  */
  public int countContours(int which){
    
    int n = 0;
    Vector<BHPolygon> v = getContourVector(which);
    if (v != null) {n = v.size();}
    return n;
  }
  

/**
  * Returns all of the contours from this and its siblings
  */
  public BHPolygon[] getContours(int which){
     
     Vector <BHPolygon> v = getContourVector(which);
     int n = v.size();
     if (n == 0){
      return null;
     } else {
      return (BHPolygon[]) v.toArray(new BHPolygon[n]);
     }
     
   //     int n = countContours(which);
   //     BHPolygon[] p = new BHPolygon[n];
   //     if (n == 0) {return p;}
   //     int ix = 0;
   //     if (which == BHPolygon.OUTER){
   //       for (int i = 0; i<opixel.size(); i++){
   //          p[ix] = opixel.get(i);
   //          ix++;
   //       }
   //     } else if {
   //       for (int i = 0; i<ipixel.size(); i++){
   //          p[ix] = ipixel.get(i);
   //          ix++;
   //       }      
   //     }
   //    return p; 
  }
  
  
/**
  * Computes the area enclosed by the specified polygons
  */
  public double computeContourArea(int which){
      double parea = 0.0;
      BHPolygon[] p = getContours(which);
      if (p != null){
         for (int i = 0; i < p.length; i++){
            parea += p[i].computeArea();
         }
      }
      return parea;
  } //computeContourArea
  


/**
  * Returns the aggregated contours points as [npoints][2]
  * @param which The type to return such as BHPolygon.OUTER
  */   
  public int[][] getContourAsXY(int which){
      
   BHPolygon[] p = getContours(which);
   int np = 0;
   for (int i = 0; i< p.length; i++){
     np += p[i].npoints;
   }
   
   int[][] xy = new int[np][2];
   int ix = 0;
   for (int i = 0; i < p.length; i++){
      System.arraycopy(p[i].getXY(), 0, xy, ix, p[i].npoints);
      ix += p[i].npoints;
   }
   return xy;
  } //getContourAsXY
/**
  * Creates a FLoatProcessor of the correct width and height and populated
  * with pixel values - background is 0.
  */ 
  public ImageProcessor createProcessor(){
    ByteProcessor ip = (ByteProcessor) super.createProcessor();
    int x0 =  xrange[0];
    int y0 =  yrange[0];
    //x0 = 0; y0 = 0;
    BHPoint[] p = getPointsAsArray();
    //if ((p == null) || (p.length  != w*h)){ return null;}
    for (int i = 0; i < p.length; i++){
        //ip.putPixel(points.points[i].x, points.points[i].y, col);
        ip.putPixelValue(p[i].x  - x0, p[i].y - y0, 255.0);
    }//i-loop
    //IJ.log("ip=" + ip);
    return ip;
  }//createProcessor  
      
}//BHSingleBlob