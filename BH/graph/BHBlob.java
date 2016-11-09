import ij.process.*;
import ij.plugin.filter.*;
import ij.*;
import ij.measure.*;
import ij.gui.*;
import java.util.Arrays;

/**
  * This class serves as a base class for BHSingleBlob and BHCompoundBlob.
  * It provides a container for the selection and results of measurments.
  * <br><p>
  * This class also provides a small number of generic methods for the computation
  * of stats.  Because I made the Blob-to-SingleBlob and Blob-to-CompoundBlob
  * an inheritance relationship rather than an interface relationship I sometimes
  * get awkward calling sequences. I often think I'll regret this choice later
  * because all of the java advice is to interface-interface-interface.
  * </p><p>
  * CONCEPTUAL FRAMEWORK:
  * A bi-level 2d image is comprised of foreground and background features.
  * Contiguous groups of pixels are considered SingleBlobs
  * CompoundBlobs are comprised of one or more SingleBlobs
  * </p><p>
  * Contact: Mike Sieracki, msieracki@bigelow.org <br>
  * Author:  Ben Tupper, btupper@bigelow.org <br>
  * History: Developed over the course of 2009 and 2010 using ImageJ 1.43 and Mac OSX 10.5.8 <br>
  *   at Bigelow Laboratory for Ocean Science, www.bigelow.org <p><p>
  * 2010-02-27 BTT switch to using Calibration for managing coordinate conversion
  */


public class BHBlob {

    /** Includes area, perim, count (of constituent blobs) */
    public static final int BASIC = 1; //blobId, area, perim, countOfSingleBlobs
    /** Include locational info: centroid, bounding box, start */
    public static final int LOCATION = 2; // centroid, bounding box, start
    /** Include shape descriptors: feret, circ, solid */
    public static final int SHAPE = 4; //feret, circ, solid, ...
    /** Include volume estimate */
    public static final int VOLUME = 8; //volume estimate
    /** Include slice number */
    public static final int SLICE = 16; //slice number 1,2,3,4,...
    /** Include the ellipse fitting parameters */
    public static final int ELLIPSE = 32; //fit parameters
    /** Include the intensity (gray) info: center of mass, moments, min, max */
    public static final int GRAY = 64; //mass, mean, etc.
    /** Include the Hu and AMI invarient moments */
    public static final int MOMENTS = 128;

    /** This flag is "+"ed with the constants above to determine what is computed/shown.
      Note that some things are computed even if not shown on the table. */
    public int measureFlag = BASIC;


    /** Edge touch flags 0 = none, 1 = left, 2 = top, 4 = right 8 = bottom */
    public static final int NONE = 0, LEFT = 1, TOP = 2, RIGHT = 4, BOTTOM = 8;

    /** Include the perimeters found by CCL perimeters, ImageJ's Wand tool or both?
      * Currently only the CCL perimeter is used.
      */
    public static final int USECCLPERIMETERS = 1;// Use the CCL perimeter
    public static final int USEWANDPERIMETERS = 2; // Use Wand (not good for SVW_VOLUME)
    public static final int USEBOTHPERIMETERS = 4; //I am not sure what this would mean
    public int usePerimeter = USECCLPERIMETERS;

    /** The integer label for this object - called "Blob" in the results table*/
    public int label;
    /** The slice number (1,2,3,...) */
    public int slice;
    /** The side(s) of the image this object touches
      * These are recorded by adding the NONE LEFT RIGHT TOP BOTTOM flags above*/
    public int side;

  /**  pixel scaling and offset info */
    Calibration cal;

    /** The area in scaled units, area is of the pixels, outerArea is area enclosed by outer boundary and
      innerArea is the area enclosed by interior holes, the latetr two are used to compute
      voidFraction*/
    public double area, outerArea, innerArea, voidFraction;
    /* Euler number */
    public int nEuler;
    /** The centroid */
    public double xc, yc; //the centroid
    /** The limits of the object as [min, max] */
    public int[] xrange = new int[2]; //limits in x and y
    public int[] yrange = new int[2];

    /* central moments */
    public double u00, u10, u01, u11, u20, u02, u21, u12, u30, u03;
    /* normalized central moments */
    public double n00, n10, n01, n11, n20, n02, n21, n12, n30, n03;

    public double xsum, ysum, x2sum, y2sum, xysum ; //for moments
    public double x2ysum, xy2sum, x3sum, y3sum;
    public double dx, dy; //not really used
    /** Moments as per Horn 1986 */
    public double m20, m11, m02, mtheta, morient, mecc, mra, mrb; // moments and orientation
    public boolean matchEllipseMoments = true;  //scale the ellipse to match moments of blob(s)?

    public double hferet, vferet; //max width and max height in rotated
    /** The number of elements of the feret array */
    protected int nferet = 7;//number of feret elements
    /** The indices of the various feret measurements */
    public static int FERETMAX = 0, FERETANGLE = 1, FERETMIN = 2, FERETX0 = 3,
      FERETY0 = 4, FERETMEAN = 5, FERETMED = 6;
    /** The feret array */
    public double[] feret; //max, angle, min, x0, y0, mean, median

    /** The Hu invarient moments [H1, ... H7] see Burger and Burge reference*/
    public int nHu = 7;
    public double[] mHu;
    public double[] gHu;

    /** The Affine Moment Invariants [AMI1, ..., AMI4] (See Flusser references)*/
    public int nAMI = 4;
    public double[] mAMI;

    /** Elliptic Fourier Descriptors see Nixon and Aguado reference */
    public int nFD = 5;  //the number of fourier descriptors
    public double[][] afd;  //a coefficients [x][y]
    public double[][] bfd;  //b coefficients [x][y]

    /** mvol is the mechanical "biovolume" estimate" using various methods
        VOL_SVW uses Sieracki, Viles and Webb (1989)
        VOL_ABD simply converts area based diameter to volume using 4/3*pi*(ABD/2)^3
        VOL_NONE skips the calculation and makes mvol=0
    */
    /* Volume estimate in scaled units - historically had "m" prefix because it was derived from mechanical moments */
    public double mvol;
    /** VOL_SKEDM estimate volume by multipling EDM by skeleton.  The result is used as an estimate of the particle radius, r, at that location. Summing disks of radius r and thickness equal to 1 pixel is the volume. VOL_SVW estimates volume with Sieracki Viles and Webb volume-by-slice integeration.  See the SVW reference.  VOL_ABD estimates the volume using 4/3*pi*r^3.  R is determined by backwards solving area=pi*r^2. VOL_NONE means don't compute volume. VOL_AUTO currently uses VOL_SKEDM where solidity < 0.4 and uses VOL_SVW otherwise. */
    public static int VOL_SKEDM = 3, VOL_SVW = 2, VOL_ABD = 1, VOL_NONE = 0, VOL_AUTO = -1;
    public int volMethod = VOL_SKEDM;
    /* Not used */
    public int volMethodUsed = VOL_SKEDM;
    /* The solidity cutoff used when volMethod = VOL_AUTO */
    public double volCutoff = 0.4; //for auto  <cutoff use VOL_SKEDM >= cutoff  use VOL_SVW
    /** Perimeter of the object */
    public double perim; //aggregate perimeter of outer contour(s)
    /** Circularity = area to perimeter ratio (4*pi*area/perim^2) */
    public double circ;
    /** Solidity = outerPerimArea/chullPerimeter */
    public double solid;
    /** Compactness = SQRT(4/pi*area)/maxDiameter */
    public double compact;
    /** Aspect Ratio = maxDiameter/minDiameter */
    public double aspect;
    /** Roundness =  (4*Area)/(pi*maxDiameter) */
    public double round; //
    /** Extent = area/boundingBoxArea */
    public double extent; //
    /** Center of mass - grayscale weighted */
    public double xcm, ycm;
    /** The [min,max] range of grayscale values */
    public double[] grange = new double[2]; //min-max
    /** grayscale mass, mean skewness and kurtosis */
    public double gmass, gmean, gvar, gskew, gkurt;
    /** Grayscale moments */
    public double gsum1, gsum2, gsum3,gsum4,gxsum,gysum; //for gray stats


/**
  * Construct with default calibration, measurements (BASIC) and no label
  */
  public BHBlob(){
    initialize(new Calibration());
  }
  /**
    * Construct with specified Calibration and measurments but no label.
    *  Measurements are provided by summing the BASIC + LOCATION + ...
    */
  public BHBlob(Calibration cal, int measurements){
   measureFlag = measurements;
   initialize(cal);
  }
/**
  * Construct with default calibration, measurements and but label
  */
  public BHBlob(int theLabel){
    label = theLabel;
    initialize(new Calibration());
  }
  /**
    * Construct with and measurmentsand label but default calibration.
    *  Measurements are provided by summing the BASIC + LOCATION + ...
    */
  public BHBlob(int theLabel, int measurements){
    label = theLabel;
    measureFlag = measurements;
    initialize(new Calibration());
  }
  /**
    * Construct with calibration and label but default measurments (BASIC).
    */

  public BHBlob(int theLabel, Calibration cal){
    label = theLabel;
    initialize(cal);
  }

  /**
    * Construct with calibration and label and measurments (BASIC +LOCATION + whatever...).
    */
  public BHBlob(int theLabel, Calibration cal, int measurements){
    label = theLabel;
    measureFlag = measurements;
    initialize(cal);
  }

/**
  * Set the slice number
  */
  public void setSlice(int theSlice){slice = theSlice;}


/**
  * Returns the options for measurements in an array
  */
  public static int[] getMeasurementOptions(){
    int[] opts =  {BHBlob.BASIC, BHBlob.LOCATION, BHBlob.SHAPE,
        BHBlob.VOLUME,  BHBlob.SLICE, BHBlob.ELLIPSE,
        BHBlob.GRAY, BHBlob.MOMENTS};
    return opts;
  }
/**
  * Returns a list of names equivalent to the measurement options
  * @param dialogNames If true then the list of names is augmented in
  * a way suitable for making a GenericDialog
  */
  public static String[] getMeasurmentOptionsNames(boolean dialogNames){
    String[] names = {"Basic", "Location", "Shape", "Volume", "Slice", "Ellipse", "Gray", "Moments"};
    if (dialogNames == true){
      names[0] = names[0] + " info";
      names[1] = names[1] + " info";
      names[2] = names[2] + " info";
      names[3] = names[3] + " estimate";
      names[4] = names[4] + " number";
      names[5] = names[5] + " parameters";
      names[6] = names[6] + "scale info";
      names[7] = names[7] + " (invariant)";
    }
    return names;
  }

/**
  * Initialize - sets everyhting to 0 and empties the containers
  * It does NOT reset the label.
  * NOTE: pixel scaling and offset info  as [units per pixel, offset in unit]. Note
  *  that realLocation = unitsPerPixel * pixelLocation + offsetInUnits.  So we might have
  *  x = xscale[0] * xPixel + xscale[1]
  * @param cal  The calibration class that provides spatial calibration.
  */
  protected void initialize(Calibration cal){
    //scaling info comes from the calibration
    this.cal = cal;
    //xscale[0] = cal.pixelWidth;
    //yscale[0] = cal.pixelHeight;
    //zscale[0] = cal.pixelDepth;
    //xscale[1] = cal.xOrigin;
    //yscale[1] = cal.yOrigin;
    //zscale[1] = cal.zOrigin;
    side = 0;
    slice = 1;

    resetFeaturesToDefaults();
  }

  /*
  * Returns the angle phi given the orinetation theta
  * @theta the angle camputed by calculation of binary moments
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
 * A utility function Returns the range from 0 to count in the array v
 * @param v the array to search
 * @param count the length of values to search within starting from 0
 * @return a 2 element array of [minV, maxV]
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
 * A convenient wrapper around range(v, count) when count is the same
 * as the length of v.
 * @param v the array to search
 * @return a 2 element array of [minV, maxV]
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
  * Returns a contour of the equivalent ellipse
  * see Imaging Book section 11.4
  * @param b A BHSingleBlob from which to compute the ellipse
  * @return A PolygonRoi (see ij.gui.PolygonRoi) with 90 vertices.
  */
 public PolygonRoi getEllipseRoi(BHSingleBlob b){

  int n = 90;
  double phi = getPhi(b.mtheta);
  double sin = Math.sin(phi);
  double cos = Math.cos(phi);
  double[] t = new double[n];
  int[] x = new int[n];
  int[] y = new int[n];
  for (int i = 0; i < n; i++){
     t[i] = Math.toRadians((double) (4*i));
     x[i] = (int) Math.round(b.xc + cos * b.mra * Math.cos(t[i]) - sin * b.mrb * Math.sin(t[i]));
     y[i] = (int) Math.round(b.yc + sin * b.mra * Math.cos(t[i]) + cos * b.mrb * Math.sin(t[i]));
     x[i] = (int) (cal.getX((double) x[i]) + 0.5);
     y[i] = (int) (cal.getY((double) y[i]) + 0.5);
     //x[i] = (int) ((double) x[i] * xscale[0] + xscale[1] + 0.5);
     //y[i] = (int) ((double) y[i] * yscale[0] + yscale[1] + 0.5);
  }

  return new PolygonRoi(x, y, n, Roi.POLYGON);

 }// getEllipseRoi


/**
  * Returns a contour of the equivalent ellipse
  * see Imaging Book section 11.4
  * @param b A BHCompoundBlob from which to compute the ellipse
  * @return A PolygonRoi (see ij.gui.PolygonRoi) with 90 vertices.
  */
 public PolygonRoi getEllipseRoi(BHCompoundBlob b){

  int n = 90;
  double phi = getPhi(b.mtheta);
  double sin = Math.sin(phi);
  double cos = Math.cos(phi);
  double[] t = new double[n];
  int[] x = new int[n];
  int[] y = new int[n];
  //IJ.log("b.mra=" + b.mra + " b.mrb = " + b.mrb);
  for (int i = 0; i < n; i++){
     t[i] = Math.toRadians((double) (4*i));
     x[i] = (int) Math.round(b.xc + cos * b.mra * Math.cos(t[i]) - sin * b.mrb * Math.sin(t[i]));
     y[i] = (int) Math.round(b.yc + sin * b.mra * Math.cos(t[i]) + cos * b.mrb * Math.sin(t[i]));
      //round the values
     x[i] = (int) (cal.getX((double) x[i]) + 0.5);
     y[i] = (int) (cal.getY((double) y[i]) + 0.5);
     //x[i] = (int) ((double) x[i] * xscale[0] + xscale[1] + 0.5);
     //y[i] = (int) ((double) y[i] * yscale[0] + yscale[1] + 0.5);
  }

  return new PolygonRoi(x, y, n, Roi.POLYGON);

 }// getEllipseRoi

/**
  * Computes the 7-element ferets - this is based heavily (almost indentically)
  * on ImageJ's ij.gui,Roi.getFeretValues().  This will return slightly different
  * values since the definition of the polygon describing the boundary
  * differs from ImageJ's Wand tool.
  * Computes [max, angle, min, x0, y0, mean, median] where x0, y0 are the points along
  * the boundary through which the max feret passes at the angle (in radians)
  * @param p the [npoints][2] array of coordinates - typically the
  *   outer contour or contours.
  * @return a 7-element array of [max, angle, min, x0, y0, mean, median]
  */
  public double[] computeFeret(int[][] p){

    double minF=Double.MAX_VALUE, maxF=0.0, angle=0.0, feretX=0.0, feretY=0.0;
    int p1=0, p2=0; //index holders for maxF locations
    double pw = cal.pixelWidth, ph = cal.pixelHeight;
    //double pw = xscale[0], ph = yscale[0];
    double w2 = pw*pw;
    double h2 = ph*ph;
    double d2r = Math.PI/180.0;
    //center of the bounding box
    double centx = (xrange[1]-xrange[0])/2.0;
    double centy = (yrange[1]-yrange[0])/2.0;
    double d = 0.0, dx = 0.0, dy = 0.0, xr = 0.0, yr = 0.0;

    //compute "the" (maximum) feret
    //retain the indices of these two points that define the min/max
    for (int i=0; i<p.length; i++) {
        for (int j=i; j<p.length; j++) {
            dx = p[i][0] - p[j][0];
            dy = p[i][1] - p[j][1];
            d = dx*dx*w2 + dy*dy*h2;
            if (d>maxF) {maxF=d; p1=i; p2=j;}
        }//j-loop
    }//i-loop
    maxF = Math.sqrt(maxF);

    //make arrays for just x and y normalized to centx and centy
    double[] x = new double[p.length];
    double[] y = new double[p.length];
    for (int i=0; i<x.length; i++) {
        x[i] = p[i][0] - centx;
        y[i] = p[i][1] - centy;
    }

    int nAngle = 180;
    double[] f = new double[nAngle+1];
    double fmean = 0.0;
    double fmedian = 0.0;
    //compute the minimum feret
    double sin = 0.0,cos = 0.0;
    double xmin = 0.0,xmax = 0.0,ymin = 0.0,ymax= 0.0;
    double tempMinF = 0.0, width = 0.0, height = 0.0;
    for (int a=0; a<=nAngle; a++) {
      cos = Math.cos(a*d2r);
      sin = Math.sin(a*d2r);
      xmin=Double.MAX_VALUE; ymin=Double.MAX_VALUE;
      xmax=-Double.MAX_VALUE; ymax=-Double.MAX_VALUE;
      for (int i=0; i<p.length; i++) {
          xr = cos*x[i] - sin*y[i];
          yr = sin*x[i] + cos*y[i];
          if (xr<xmin) xmin = xr;
          if (xr>xmax) xmax = xr;
          if (yr<ymin) ymin = yr;
          if (yr>ymax) ymax = yr;
      }//i-lopp
      width = xmax - xmin;
      height = ymax - ymin;
      tempMinF = Math.min(width, height);
      f[a] = tempMinF;
      fmean += f[a];
      minF = Math.min(minF, tempMinF);
    }//a-loop

    fmean = fmean/nAngle;
    //get the median
    Arrays.sort(f);
    int middle = f.length/2;
    if (f.length%2 == 1) {
        fmedian= f[middle];
    } else {
       fmedian = (f[middle-1] + f[middle]) / 2.0;
    }


    //these are the integer versions of the feretX and feretY
    // sorted into "ascending" order
    double x1=p[p1][0];
    double y1=p[p1][1];
    double x2=p[p2][0];
    double y2=p[p2][1];
    if (x1>x2) {
        double tx1=x1, ty1=y1;
        x1=x2; y1=y2; x2=tx1; y2=ty1;
    }
    feretX = x1*pw;
    feretY = y1*ph;
    dx=x2-x1; dy=y1-y2;
    angle = Math.atan2(dy*ph, dx*pw);
    if (angle<0) angle = Math.PI + angle;
    if (pw==ph) {minF *= pw;}
    double[] arr = new double[nferet];
    arr[0] = maxF;
    arr[1] = angle;
    arr[2] = minF;
    arr[3] = feretX;
    arr[4] = feretY;
    arr[5] = fmean;
    arr[6] = fmedian;
    return arr;
  } //computeFeret


/**
  * Creates a processor of the correct width and height and populated
  * with pixel values - background is 0.  The edge padding defaults to 4.
  */
  public ImageProcessor createProcessor(){
    return createProcessor(4);
  }

/**
  * Create a blank processor of the suitable width and height
  * @param pad the number of pixels to pad around the object
  */
  public ImageProcessor createProcessor(int pad){
    int width = xrange[1] - xrange[0] + 1 + (2*pad);
    int height = yrange[1] - yrange[0] + 1 + (2 * pad);
    if ((width <= (2*pad)) || (height <= (2*pad))) { return null;}
    ByteProcessor fip = new ByteProcessor(width, height);
    //IJ.log("w=" + width + " h=" + height + " fip=" + fip);
    return fip;
  } //createProcessor


/**
  * Returns the appropriate touching flag for the specifed width and height
  *
  */
  public int isTouching( int w, int h){
    return isTouching (0,0,w,h);
  }
/**
  * Returns the appropriate touching flag for the specifed width and height
  * @param xs the staring coordinate of the image in pixels (generally 0)
  * @param ys the staring coordinate of the image in pixels (generally 0)
  * @param w the width of the image in pixels
  * @param h the height of the image in pixels
  * @return an integer combination of NONE, LEFT, RIGHT, TOP and BOTTOM indicating
  * where the object encounters the edge of the image
  */
  public int isTouching(int xs, int ys, int w, int h){
    int r = NONE;
    if (xrange[0] <= xs) {r += LEFT;}
    if (xrange[1] >= (xs + w -1)){ r += RIGHT;}
    if (yrange[0] <= ys) {r += TOP;}
    if (yrange[1] >= (ys + h - 1)) { r += BOTTOM;}
    return r;
  }

  /**
    * A method for rounding that avoids rounding to the nearest
    * even number - found in Math.rint(double) - so this required for application
    * of the SVW volume estimate.
    */
   public double rint(double value){
      return Math.floor(value + 0.5);
   }

/**
  * Computes the volume using the product of the skeletonization (SK) and the
  * Euclidean  Distance Map (EDM) for BHSingleBlob class objects.
  * @param blob the BHSingleBlob object
  */
  protected double computeVol_SKEDM(BHSingleBlob blob){

   ByteProcessor edmIp = (ByteProcessor) blob.createProcessor();

   EDM EDMengine = new EDM();
   FloatProcessor edm = EDMengine.makeFloatEDM(edmIp, 0, true);

   ByteProcessor skel = (ByteProcessor) blob.createProcessor();
   skel.invertLut();
   skel.skeletonize();

   float maxV =  (new Double(skel.getMax())).floatValue();
   double v = 0.0;
   mvol = 0.0;
   double extra = 0.25;  //something to play with regarding the width of the edm
   double p3 = cal.getX(1.0) * cal.getY(1.0) * cal.getZ(1.0);
   //double p3 = xscale[0] * xscale[0] * xscale[0];
	 for (int y = 0; y < skel.getHeight(); y++){
	   for (int x = 0; x < skel.getWidth(); x++){
		    if (skel.getPixelValue(x,y) != 0){
		     v = (edm.getPixelValue(x,y) + extra);
		     mvol +=  Math.PI * v * v * p3;  //like a stack of discs
		   } //is there a skelton item?
	   }	//x-loop
	 } 	//y-loop
	 return mvol;
  } //compueVol_SKEDM


 /**
   * A wrapper around computeMoment where x and y are integers
   */
 protected double computeMoment(int x, double x0, int y, double y0, int p, int q, double value){
    return computeMoment((double) x, x0, (double) y, y0, p, q, value);
 }
 /**
   * Computes the p-q order centralized moments of x and y about the centroid x0, y0
   * result = value * (x-x0)^p * (y-y0)^q
   */
 protected double computeMoment(double x, double x0, double y, double y0, int p, int q, double value){
    return value * Math.pow(x-x0, (double) p) * Math.pow(y-y0, (double) q);
 }
 /**
   * Computes normalized centralized moments
   */
 protected double computeNormMoment(double u, double u0, int p, int q){
    return u * Math.pow(1.0/ u0, (double) (p + q + 2) / 2.0);
 }


/**
  * Computes the the moments, centralized moments, normalized centralized moments
  * and Hu moment invariants described in Gonzales and Woods 3d edition and elsewhere
  * NOTE 2009-10-20  There is some duplication between the "old" system and this new system
  * the old could be removed... but be careful about compatability issues.
  * Essentially, the uXX values can replace the xsum, etc values
  * @param p An array of BHPpoints
  */
 protected void computeMoments(BHPoint[] p){

    double dv, dv2;
    //reset these just to make sure
    gsum1 = gsum2 = gsum3 = gsum4 = gxsum = gysum = 0.0;
    xsum = ysum = xysum = x2sum = y2sum = 0.0;
    u00=0.0; u10=0.0; u01=0.0; u11=0.0; u20=0.0; u02=0.0; u21=0.0; u12=0.0; u30=0.0; u03=0.0;
    n20=0.0; n02=0.0; n21=0.0; n12=0.0; n30=0.0; n03=0.0;


    mHu = new double[nHu];
    mAMI = new double[nAMI];

    for (int i = 0; i < p.length; i++){
      xsum += (double) p[i].x;
      ysum += (double) p[i].y;
      x2sum += (double) (p[i].x * p[i].x);
      y2sum += (double) (p[i].y * p[i].y);
      xysum += (double) (p[i].x * p[i].y);
      if (p[i].x < xrange[0]) { xrange[0] = p[i].x;}
      if (p[i].x > xrange[1]) { xrange[1] = p[i].x;}
      if (p[i].y < yrange[0]) { yrange[0] = p[i].y;}
      if (p[i].y > yrange[1]) { yrange[1] = p[i].y;}

      dv = p[i].val; //the value as a double
      if (dv < grange[0]) {grange[0] = dv;}
      if (dv > grange[1]) {grange[1] = dv;}
      dv2 = dv*dv; //its square
      gsum1 += dv; //sum of values
      gsum2 += dv2; //sum of squares
      gsum3 += dv*dv2; //third power
      gsum4 += dv2*dv2;//fourth power
      gxsum += p[i].x*dv;//weighted location
      gysum += p[i].y*dv;//ditto
    }

    //area
    u00 = (double) p.length;
    //central moments
    xc = xsum/u00;
    yc = ysum/u00;
    for (int i = 0; i < p.length; i++){
        //u00 += computeMoment(p[i].x, xc, p[i].y, yc, 0, 0, 1.0);
        u11 += computeMoment(p[i].x, xc, p[i].y, yc, 1, 1, 1.0);
        u10 += computeMoment(p[i].x, xc, p[i].y, yc, 1, 0, 1.0);
        u01 += computeMoment(p[i].x, xc, p[i].y, yc, 0, 1, 1.0);
        u20 += computeMoment(p[i].x, xc, p[i].y, yc, 2, 0, 1.0);
        u02 += computeMoment(p[i].x, xc, p[i].y, yc, 0, 2, 1.0);
        u21 += computeMoment(p[i].x, xc, p[i].y, yc, 2, 1, 1.0);
        u12 += computeMoment(p[i].x, xc, p[i].y, yc, 1, 2, 1.0);
        u30 += computeMoment(p[i].x, xc, p[i].y, yc, 3, 0, 1.0);
        u03 += computeMoment(p[i].x, xc, p[i].y, yc, 0, 3, 1.0);
    } //i loop

    //normalized moments
    n20 = computeNormMoment(u20, u00, 2,0);
    n02 = computeNormMoment(u02, u00, 0,2);
    n11 = computeNormMoment(u11, u00, 1,1);
    n21 = computeNormMoment(u21, u00, 2,1);
    n12 = computeNormMoment(u12, u00, 1,2);
    n30 = computeNormMoment(u30, u00, 3,0);
    n03 = computeNormMoment(u03, u00, 0,3);

    //Hu normalized moments
    mHu[0] = n20 + n02;
    mHu[1] = Math.pow(n20 - n02,2.0) + 4.0 * Math.pow(n11, 2.0);
    mHu[2] = Math.pow(n30 - 3.0 * n12, 2.0) + Math.pow(3.0 * n21 - n03, 2.0);
    mHu[3] = Math.pow(n30 + n12, 2.0) + Math.pow(n21 + n03, 2.0);
    mHu[4] = (n30 - 3.0 * n12)*(n30 + n12)*( Math.pow(n30+n12,2.0) - 3.0 * Math.pow(n21+n03,2.0)) +
             (3.0 * n21 - n03)*(n21 + n03)*(3.0*Math.pow(n30 + n12,2.0) - Math.pow(n21 + n03,2.0));
    mHu[5] = (n20 - n02)*(Math.pow(n30+n12,2.0) - Math.pow(n21 + n03,2.0)) +
             4.0*n11*(n30+n12)*(n21+n03);
    mHu[6] = (3.0*n21-n03)*(n30+n12)*(Math.pow(n30+n12,2.0) - 3.0 * Math.pow(n21+n03,2.0)) +
             (3.0*n12-n03)*(n21+n03)*(3.0*Math.pow(n30+n12,2.0) - Math.pow(n21+n03,2.0));

    //AMI (affine moment invariants)
    mAMI[0] = (u20*u02-u11*u11)/Math.pow(u00,4.0);
    mAMI[1] = (u30*u30*u03*u03 - 6.0*u30*u21*u12*u03 + 4.0*u30*Math.pow(u12,3.0) +
      4.0*u03*Math.pow(u21,3.0) - 3.0*u21*u21*u12*u12)/Math.pow(u00,10.0);
    mAMI[2] = (u20*(u21*u03-u12*u12) - u11*(u30*u03-u21*u12) +
               u02*(u30*u12-u21*u21) )/Math.pow(u00, 7.0);
    mAMI[3] = (Math.pow(u20,3.0)*u03*u30 - 6.0*u20*u20*u11*u12*u03 - 6.0*u20*u20*u02*u21*u03 +
      9.0*u20*u20*u02*u12*u12 + 12.0*u20*u11*u11*u21*u03 + 6.0*u20*u11*u02*u30*u03 -
      18.0*u20*u11*u02*u21*u12 - 8.0*Math.pow(u11,3.0)*u30*u03 - 6.0*u20*u02*u02*u30*u12 +
      9.0*u20*u02*u02*u21*u21 + 12.0*u11*u11*u02*u30*u12 - 6.0*u11*u02*u02*u30*u21 +
      Math.pow(u02,3.0)*u30*u30)/Math.pow(u00, 11.0);

 }//computeMoments


/**
  * Computes the ellipse semi-major and semi-minor axes as well as the
  * eccentricity.  Uses current state of matchEllipseMoments to enable/disable
  * axes scaling of ellipse.
  */
  protected void computeEllipse(){ computeEllipse(matchEllipseMoments); }

/**
  * Computes the ellipse semi-major and semi-minor axes as well as the
  * eccentricity.
  * @param scale if true causes the axes to be normalized so that the central
  * moments of the ellipse match the original blob
  */
  protected void computeEllipse(boolean scale){


      //see Imaging Book 11.4 Properties of Binary Regions
      double a1 = u20 + u02 + Math.sqrt(Math.pow((u20-u02),2.0) + 4.0 * u11*u11);
      double a2 = u20 + u02 - Math.sqrt(Math.pow((u20-u02),2.0) + 4.0 * u11*u11);
      mra = cal.getX(Math.sqrt(2.0*a1/u00));// semi-major
      mrb = cal.getX(Math.sqrt(2.0*a2/u00));//semi-minor
      mecc = a1/a2; //eccentricity

      if (scale) {

         //see the description in ImageJ's EllipseFitter
         double s = Math.sqrt(u00/(Math.PI * mra * mrb));
         mra = mra * s;
         mrb = mrb * s;

      }

  }//computeEllipse

/**
  * Resets the various features to default values
  */
  public void resetFeaturesToDefaults(){
    area = 0.0;
    outerArea = 0.0;
    innerArea = 0.0;
    voidFraction = 0.0;
    nEuler = 0;
    xc = yc = 0.0; //centroid
    xcm = ycm = 0.0; //center of mass (gray)
    xsum = x2sum = 0.0;
    ysum = y2sum = 0.0;
    x2ysum = xy2sum = x3sum = y3sum = 0.0;
    u00 = u10 = u01 = u20 = u02 = u12 = u21 = u30 = u30 = 0.0;
    n00 = n10 = n01 = n20 = n02 = n12 = n21 = n30 = n30 = 0.0;
    xysum = 0.0;
    dx = dy = 0.0;
    m20 = m02 = m11 = mtheta = morient = mecc = mra = mrb = 0.0;
    hferet = vferet = 0.0;
    feret = new double[nferet];
    mHu = new double[nHu];
    gHu = new double[nHu];
    mAMI = new double[nAMI];
    mvol = 0.0;
    perim = 0.0;
    circ = 1.0;
    solid = compact = aspect = round  = extent = 0.0;
    xrange[0] = Integer.MAX_VALUE; xrange[1] = Integer.MIN_VALUE;
    yrange[0] = Integer.MAX_VALUE; yrange[1] = Integer.MIN_VALUE;
    grange[0] = Double.MAX_VALUE; grange[1] = -Double.MAX_VALUE;
    gmass = gmean = gvar = gskew = gkurt = 0.0;
    gsum1 = gsum2 = gsum3 = gsum4 = gxsum = gysum = 0.0;
  }

/**
  * Similar to getBoundingBox, but ALWAYS returns pixel values
  * of [x0, y0, w, h]
  */
  public int[] getBoundingRect(){
    int[] b = new int[4];
    b[0] = xrange[0];
    b[1] = yrange[0];
    b[2] = xrange[1]-xrange[0];
    b[3] = yrange[1]-yrange[0];
    return b;
  }//getBoundingRect

/**
  * Returns the bounding box in pixel units as [x0, y0, w, h]
  */
  public double[] getBoundingBox(){
    return getBoundingBox(true);
  } //getBoundingBox

/**
  * Returns the bounding box as [x0, y0, w, h]
  * @param asPixels If true then return pixel values otherwise scaled values
  */
  public double[] getBoundingBox(boolean asPixels){

      double[] b = new double[4];
//      b[0] = xrange[0];
//      b[1] = yrange[0];
//      b[2] = xrange[1]-xrange[0];
//      b[3] = yrange[1]-yrange[0];

      if (!asPixels) {
         b[0] = cal.getX(xrange[0]);
         b[1] = cal.getY(yrange[0]);
         b[2] = cal.getX(xrange[1])-b[0];
         b[3] = cal.getY(yrange[1])-b[1];
      } else {
         b[0] = xrange[0];
         b[1] = yrange[0];
         b[2] = xrange[1]-xrange[0];
         b[3] = yrange[1]-yrange[0];
      }

      return b;
  }//getBoundingBox

  /**
    * Returns the number of pixels in this blob - in this case 0
    */
  public int npoints(){return 0;}
  /**
    * Returns the convex hull as null
    */
  public BHPolygon getConvexHull(){return null;}
  /**
    * Computes area - returns 0
    */
  public double computeArea (int val) {return 0.0;}


}//BHBlob