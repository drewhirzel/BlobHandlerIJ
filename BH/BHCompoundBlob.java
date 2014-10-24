import java.util.Vector;
import ij.process.*;
import ij.plugin.filter.*;
import ij.*;
import ij.measure.*;
import ij.gui.*;


/**
  * This class is a container for related BHSingleBlobs. As things have evolved it became
  * apparent that this is the primary class for user interface.  It allows users
  * to gather groups related single blobs, like the parts of a planktonic chain.
  * <p>
  * The ID label of the CompoundBlob is the ID label of the _first_ single blob contained.
  * <p> 
  * Outer and Inner contours (CCL or crack when implemented) are those of the 
  * constituent SingleBlobs.
  * <p> 
  * The convex hull of the Compound Blob is the convex hull of the constituent 
  * outer boundaries.
  * <p>
  * Contact: Mike Sieracki, msieracki@bigelow.org <br>
  * Author:  Ben Tupper, btupper@bigelow.org <br>
  * History: Developed over the course of 2009 and 2010 using ImageJ 1.43 and Mac OSX 10.5.8 <br>
  *   at Bigelow Laboratory for Ocean Science, www.bigelow.org<p><p>
  * 2010-02-27 BTT convert to using Calibration for coordinate conversion
  * 2010-04-07 BTT added "crack" boundary stuff
  **/
public class BHCompoundBlob extends BHBlob {

    /**  This is where all of the SingleBlobs are stored */
    private Vector<BHSingleBlob> vData = new Vector<BHSingleBlob>();
    
    /** The convex hull of the compound blob */
    public BHPolygon chull;
     
  /**
    * Contructs an empty blob
    */
  public BHCompoundBlob(){
    initialize();
  }
    
  public BHCompoundBlob(int theLabel){
    label = theLabel;
    initialize();
  }

/**
  * Constructs the CompoundBlob and adds its first SingleBlob 
  * @param b the first blob to add - the CompundBlob gets its label and 
  * spatial scaling info from this SingleBlob.
  */
  public BHCompoundBlob(BHSingleBlob b){
    label = b.label;
    initialize();
    add(b);
  }    
    
  

/**
  * Initialize - sets everything to 0 and empties the containers
  * It does _NOT_ reset the label.
  */
  protected void initialize(){
    super.initialize(new Calibration());        
    //clear out the siblings
    if (vData != null){
      vData.clear();
    } else {
      vData = new Vector<BHSingleBlob>();
    } 
  }//initialize  
  
 /**
  * Sets the volume calculation method
  * See BHSingleBlob for details
  */
  public void setVolumeMethod(int theMethod){
    volMethod = theMethod;//set your own method
    BHSingleBlob sb = null;
    for (int i = 0; i < vData.size(); i++){
        sb = vData.get(i);
        if (sb != null){
          sb.volMethod = theMethod;//set your children's
          //sb.computeStats();
        }
    } //i-loop
  }//setVolumeMethod
    

/**
  * This method computes the stats for this CompundBlob.
  * 1. reset the values
  * 2. if 0 SingleBlob then return with values reset
  * 3. if one or more then compute the aggregate values
  *
  * grayscale stats are computes just like in ImageJ 
  * (see ByteStatistics.calculateMoments for example)
  */
  public void computeStats(){

    resetFeaturesToDefaults();
    int nBlobs = count();
    if (nBlobs == 0){return;} //work is done here
    
    BHSingleBlob aBlob = vData.get(0);
    int[] labels = getBlobLabels();
    for (int i = 0; i < labels.length; i++){
      aBlob = getByLabel(labels[i]); 
      area += aBlob.area;
      //outerArea, enclosed by outer contour
      outerArea += aBlob.outerArea;
      //innerArea, this is the void area (may be 0)
      innerArea += aBlob.innerArea;
      //nEuler
      nEuler += aBlob.nEuler;
    }
      //voidFraction, this is the ratio of inner/outer
    voidFraction = innerArea/outerArea;  
 
    aBlob = vData.get(0);
   
    BHPoint[] p = getPoints();
    int np = npoints();

    super.computeMoments(getPoints());
    xc = cal.getX(xsum/np + 0.5);
    yc = cal.getY(ysum/np + 0.5);
      
    gmass = gsum1; 
    gmean = gmass/np;
    double mean2 = gmean*gmean;
    gvar = gsum2/np - mean2;
    double sDeviation = Math.sqrt(gvar);      
    gskew = ((gsum3 - 3.0*gmean*gsum2)/np 
      + 2.0*gmean*mean2)/(gvar*sDeviation);
    gkurt = (((gsum4 - 4.0*gmean*gsum3 + 6.0*mean2*gsum2)/np - 
      3.0*mean2*mean2)/(gvar*gvar)-3.0);
    xcm = cal.getX(gxsum/gsum1 + 0.5); 
    ycm = cal.getY(gysum/gsum1 + 0.5);   
    
    chull = getConvexHull();    
    //perim = getOuterPerimeterLength();
    perim = getBoundaryLength(BHPolygon.OUTER + BHPolygon.CRACK);
    
    // these match IPLab with very good precision
    m11 = xysum - (xsum*ysum)/ np;
    m20 = x2sum - (xsum*xsum) / np;
    m02 = y2sum - (ysum*ysum) / np;
    mtheta = getAngle(m20, m11, m02);
    
    if (m02 != m20) {
      morient = 0.5 * Math.atan(2*m11/(m20-m02));
    }
      
    computeEllipse(true);
    //see Imaging Book 11.4 Properties of Binary Regions 
    //double a1 = u20 + u02 + Math.sqrt(Math.pow((u20-u02),2.0) + 4.0 * u11*u11);
    //double a2 = u20 + u02 - Math.sqrt(Math.pow((u20-u02),2.0) + 4.0 * u11*u11);
    //mra = cal.getX(Math.sqrt(2.0*a1/u00));
    //mrb = cal.getX(Math.sqrt(2.0*a2/u00));
    //mecc = a1/a2;
    //mecc = Math.sqrt(1.0 - Math.pow(mrb/mra,2.0));
    
    feret = computeFeret(getContourAsXY(BHPolygon.OUTER + BHPolygon.CRACK));
      //many of these can be found in any image processing text.
      // e.g. see table 10.2 in Russ 3 ed.
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
    round = (4.0 * area) / 
      (Math.PI*Math.pow(feret[FERETMAX], 2.0));
  } //computeStats
  
/**
  * Computes the volume of the object by summing the volume of the constituent
  * objects
  */
  public void computeVolume(){
    mvol = 0.0;
    
    BHSingleBlob b = null;
    for (int i = 0; i < vData.size(); i++){
      b = vData.get(i);
      //if (b != null) { mvol += b.computeVolume();}
      if (b != null) {
         b.computeVolume();
         mvol += b.computeVolume();
      }
      
    }
  }//computeVol
  
     
/**
  * Counts the number of non-null singleblobs
  */  
  public int count(){
    int count = 0;
    for (int i = 0; i< vData.size(); i++){
      if (vData.get(i) != null){count++;}
    }
    return count;
  }//count

/**
  * Returns the SingleBlob at index specified
  */
  public BHSingleBlob get(int index){
    return vData.get(index);
  } //get (by index)

/**
  * Returns the blob with the specified label
  */
  public BHSingleBlob getByLabel(int theLabel){
    BHSingleBlob b = null;
    for (int i = 0; i < vData.size(); i++){
      if (vData.get(i).label == theLabel){
        b = vData.get(i);
        break;
      }
    }
    return b;
  }// getByLabel
  
  
 /**
   * Returns the unique labels of the SingleBlobs
   */
   private int[] getBlobLabels(){
      int[] labels = new int[vData.size()];
      int count = 0;
      BHBlob blob;
      for (int i = 0; i < labels.length; i++){
         labels[i] = -1;
         blob = vData.get(i);
         if (blob != null){
            labels[i] = blob.label;
            count++;
         }
      }
      if (count == labels.length){
         return labels;
      }
      
      int[] x = new int[count];
      int j = 0;
      for (int i = 0; i < labels.length; i++){
         if (labels[i] != -1) {
            x[j] = labels[i];
            j++;
         }
      }  
      return x;
   } 
/**
  * Adds a SingleBlob.  If it is the first one entered then this Blob
  * takes on the singleBlob's label
  */  
  public void add(BHSingleBlob b){
    if (vData.size() == 0){
      label = b.label;
      slice = b.slice;
      this.cal = b.cal;
      measureFlag = b.measureFlag;
      
    }
    vData.add(b);
  } //add a BHSingleBlob


/**
  * Adds the contents of one blob to this
  */
  public void add(BHCompoundBlob b){
    if (vData.size()== 0){
     label = b.label;
    }
    for (int i = 0; i < b.count(); i++){
      vData.add(b.get(i)); 
    }     
  }      

/**
  * Adds a point to the SingleBlob that has the same label
  */
  public int addPoint(BHPoint p){
    int ret = 0;
    BHSingleBlob b = getBlobByLabel(p.label);
    if (b != null) {
      ret = b.addPoint(p);
    } 
    return ret;
  }//addPoint
  
  
/**
  * Adds a contour to the SiingleBlob that has the same label
  */
  public void addContour(BHPolygon p){
    BHSingleBlob b = getBlobByLabel(p.label);
    if (b != null) {
      b.addContour(p); 
    }
  }//addContour
  
/**
  * Translate all of the single blobs by dx,dy
  */
  public void translate(int dx, int dy){
   BHSingleBlob sb = null;
   for (int i = 0; i< vData.size(); i++){
      sb = vData.get(i);
      if (sb != null){ sb.translate(dx,dy);}
   } 
   computeStats();
  }//translate
  

/**
  * Labels pixels in image with "label" color
  */
  public void labelContours(ImageProcessor ip, int which){
    labelContours(ip, label, which);
  } 
/**
  * Label contours in image with specified color
  */
  public void labelContours(ImageProcessor ip, int col, int which){
    for (int i = 0; i < vData.size(); i++){
      vData.get(i).labelContours(ip, col, which);
    }
  }//labelContours



/**
  * Labels pixels in image with "label" color
  */
  public void labelPixels(ImageProcessor ip){
      labelPixels(ip, label);
  }

/**
  * Labels pixels in image with specified color
  */
  public void labelPixels(ImageProcessor ip, int col){
    for (int i = 0; i < vData.size(); i++){
      vData.get(i).labelPixels(ip, col);
    }
  }
  
  
  
/**
  * Gets a BHSingleBlob by label
  */
  public BHSingleBlob getBlobByLabel(int theLabel){
    BHSingleBlob blob = null;
    for (int i=0; i< vData.size(); i++){
      if ( vData.get(i).label == theLabel){
        blob = vData.get(i);
        break;
      }
    }
    return blob;
  } //getByLabel  
/**
  * Gets a BHSingleBlob by index
  */
  public BHSingleBlob getBlob(int index){
    if (index >= vData.size()) {return null;}
    return vData.get(index);
  }
/**
  * Returns the aggregate number of points
  */
  public int npoints(){
    int np = 0;
    for (int i = 0; i<vData.size(); i++){
      //  IJ.log("label=" +(vData.get(i)).label + " n = " + (vData.get(i)).points.npoints);
        np += (vData.get(i)).points.npoints;
    }
    return np;
  } //npoints 
  
  

  
/**
  * Returns the number of outer/inner contours that define this blob
  */
  public int countContours(int which){
    int n = 0;
    for (int i = 0; i < vData.size(); i++){
      n += vData.get(i).countContours(which);
    }
    return n;
  }//countContours
  
/**
  * Returns an array of inner/outer contours that define this blob
  * @param which the flag indicating which type of contour to return.  If this
  * flag includes BHPolygon.CHULL then that is returned.
  */
  public BHPolygon[] getContours(int which){
  
     BHPolygon[] p;
     
     if ((which & BHPolygon.CHULL) != 0){
         p = new BHPolygon[1];
         p[0] = getConvexHull(which);
         return p;
     }
     p = new BHPolygon[countContours(which)];
     //IJ.log("N contours = " + p.length);
     int ix = 0;
     for (int i = 0; i< vData.size(); i++){
      BHPolygon[] p2 = vData.get(i).getContours(which);
      if ((p2 != null) && (p2.length > 0)) {
      //    IJ.log("Contour " + i + " n=" + p2.length);
          System.arraycopy(p2, 0, p, ix, p2.length);
          ix += p2.length;
      }    
     } // i-loop
     return p;
  } //getContours 
  

/**
  * Returns the all the points from all the BHSingleBlob objects
  */
  public BHPoint[] getPoints(){
    int np = npoints();
    BHPoint[] p = new BHPoint[np];
    int ix = 0;
    for (int i = 0; i< vData.size(); i++){
     BHSingleBlob b = vData.get(i);
     System.arraycopy(b.points.getArray(), 0, p, ix,b.points.npoints);
     ix += b.points.npoints;
    }//i-loop
    return p;
  }//getPoints
  
 
/**
  * Returns the aggregated contours points as [npoints][2]
  * @param which The type to return such as BHPolygon.OUTER
  */   
  public int[][] getContourAsXY(int which){
   
   if (count() == 0){return null;}
   
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
  * Returns the convex hull of this blob using the OUTER + CRACK polygon
  */  
  public BHPolygon getConvexHull(){ 
   return getConvexHull(BHPolygon.OUTER + BHPolygon.CRACK);
  }
  
/**
  * Returns the convex hull of this blob using the specifed boundary
  */
  public BHPolygon getConvexHull(int which){
   
     // if emtpy then return
    if (vData.size() == 0){ return null;}
    //if only one item - get a copy of its convex hull
    if (vData.size() == 1) {return vData.get(0).getConvexHull(which);}
    //otherwise build one
    int n = 0;
    int ix = 0; 
    BHPolygon[] poly;
    //first get the total number of outer contour points
    for (int i = 0; i < vData.size(); i++){
      poly = (vData.get(i)).getContours(which);
      if (poly != null){
         n += poly[0].npoints;
      }
    } // i-loop
    BHSingleBlob blob = vData.get(0);
    BHPointGroup pg;
    int[] x = new int[n];
    int[] y = new int[n];
    int[] tempx;
    int[] tempy;
    int[] temp; 
    for (int i = 0; i < vData.size(); i++){
      poly = (vData.get(i)).getContours(which);
      pg = poly[0].getSortedPointGroup();
      //tempx = vData.get(i).opixel.get(0).getX();
      //tempy = vData.get(i).opixel.get(0).getY();
      tempx = pg.getX(false);
      tempy = pg.getY(false);
      System.arraycopy(tempx, 0, x, ix, tempx.length);
      System.arraycopy(tempy, 0, y, ix, tempy.length);
      ix += tempx.length;
    }  //i-loop 
    
    //create a new hull using (any) polygon
    BHPolygon hull = vData.get(0).opixel.get(0).getConvexHull(x,y);
    if (hull == null){
      IJ.log("BHCompoundBlob.getConvexHull issue.  Label=" + blob.label + " nblob=" + vData.size());
    } else {
       hull.label = blob.label;
       hull.cal = this.cal;
       //poly.xscale[0] = blob.xscale[0]; poly.xscale[1] = blob.xscale[1];
       //poly.yscale[0] = blob.yscale[0]; poly.yscale[1] = blob.yscale[1];
       //poly.zscale[0] = blob.zscale[0]; poly.zscale[1] = blob.zscale[1];
    }
    return hull;
  }//getConvexHull
  
  
/**
  * Computes the area using the OUTER+CRACK boundary
  */
   public double computeArea(){ 
      return computeArea(BHPolygon.OUTER + BHPolygon.CRACK);
   }
   
/**
  * Returns the summed area enclosed by the inner or outer polygons of each
  * member blob.
  */
  public double computeArea(int which){
    
    double pArea = 0.0;
    for (int i = 0; i< vData.size(); i++){
      pArea += vData.get(i).computeArea(which); 
    }
    return pArea;
  }
  
/**
  * Returns the minimum euclidean distance between this blob on the one provided
  * using the boundary specifed OUTER + CRACK
  */
  public double computeMinimumDistanceTo(BHCompoundBlob b){
   return computeMinimumDistanceTo(b, BHPolygon.OUTER + BHPolygon.CRACK);
  }
/**
  * Returns the minimum euclidean distance between this blob on the one provided
  * using the boundary specifed by which
  */
  public double computeMinimumDistanceTo(BHCompoundBlob b, int which){
    if ((b.countContours(which) == 0) || (countContours(which)==0)){ 
      return Double.NaN;
    }
    double d = Double.MAX_VALUE;
    double temp = 0.0;
    BHPolygon[] b1 = getContours(which);
    BHPolygon[] b2 = b.getContours(which);
    for (int i = 0 ; i < b1.length; i++){
      for (int j = 0; j < b2.length; j++){
        temp = b1[i].computeMinimumDistanceTo(b2[j]);
        if (temp < d) { d = temp;}
      }//j-loop
    }//i-loop
    return d;
  }//computeMinimumDistanceTo
  
/**
  * Returns the outer perimeter length
  */
  public double getBoundaryLength(int which){
    double p = 0.0;
    for (int i = 0; i < vData.size(); i++){ 
      p += (vData.get(i)).getBoundaryLength(which);
    }
    return p;
  } //getOuterPerimeterLength
  
  
/**
  * Returns true if this holds the specified SingleBlob
  */
  public boolean hasLabel(int theLabel){
   boolean ok =  false;
   for (int i = 0; i <vData.size(); i++){
    ok = (vData.get(i).label == theLabel);
    if (ok) break;
   }
   return ok; 
  }
 
/**
  * Post the results to the provided table at the specified row
  */
  public void showInfo(ResultsTable rt, int row){
    
      if (vData.size() > 0) {
        BHPoint pointOne = vData.get(0).points.getStart();
        if ((measureFlag & BASIC) != 0){
          rt.setValue("Blob", row, (double) label);
          rt.setValue("Count", row, vData.size());
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
     	  for (int k = 0; k < vData.size(); k++){
  	   	  BHSingleBlob sblob = (BHSingleBlob) vData.get(k);
  		     if (sblob.solid < 0.4) {ntricky++;}
  		  } //k-loop
  		  rt.setValue("NTricky", row, ntricky);
  		 }//vol
  		 
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
      }
  }
  
  public void showInfo(ResultsTable rt){
    showInfo(rt, rt.getCounter()); 
  }


/**
  * Returns a contour of the equivalent ellipse
  */   
 public PolygonRoi getEllipseRoi(){
  return getEllipseRoi(this);
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
  * Creates a FLoatProcessor of the correct width and height and populated
  * with pixel values - background is 0.  The edge padding defaults to 4.
  */ 
  public ImageProcessor createProcessor(){
    return createProcessor(4);
  }
  
/**
  * Creates a processor of the correct width and height and populated
  * with pixel values - background is 0.  Edge padding with background
  * is as specified.
  * @param pad the number of background pixels to pad around the edges
  */ 
  public ImageProcessor createProcessor(int pad){
    ImageProcessor ip = super.createProcessor(pad);
    if (ip == null) { return null;}
    int x0 =  xrange[0];
    int y0 =  yrange[0];
    BHPoint[] p = getPoints();
        for (int i=0; i< p.length; i++){
           ip.putPixelValue(p[i].x - x0 + pad, p[i].y - y0 + pad, 255.0 );
    }//y-loop
    return ip;
  }//createProcessor
}//BHCompoundBlob