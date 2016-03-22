import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;
import ij.measure.*;
import ij.util.Tools;
/**
  * This class serves as a GUI interface to the general BlobHandler class which, 
  * like Image's Analyzer, extracts and measures foreground objects (aka blobs) found
  * in a bi-level image.  Measurements may be redirected to another greyscale image of the same size.
  * Results are shown in the Results window.  The chief advantages of this class 
  * over ParticleAnalyzer is that it is easily extensible, it is easy to work
  * with programmatically, and blobs may be combined to form compound blobs.
  * <p>
  * Users should note that blob perimeters (outer and inner) are defined by the 
  * coordinate of the foreground pixels next to the backgound.  This is different
  * than the boundary found using ImageJ's Wand class - that boundary follows
  * the crack between foreground and background pixels.  This results in small
  * differences between features measured with ParticleAnalyzer and BlobHandler.
  * Future plans include the option to use either the pixel or the crack boundary.
  * The difference is larger for small objects, and for 1 pixel objects some of 
  * BlobHandler's features are Infinity or NaN.
  *  
  * <p>
  * It is not necessary to use the BlobHandler_ class to run BlobHandler programmatically.  But 
  * this class is macro-ready and provides many of the commmon controls including... <br>
  * <ul>
  * <li>Filtering by size, min-max or just min (default = 0-Infinity)
  * <li>Filtering by circularity, min-max or just min (default = 0-1)
  * <li>Combining objects by nearest distance (default = 0, no combinations)
  * <li>Selection of the Volume estimate algorithm (if Volume is selected as a feature below, default = auto)
  * <li>Selection of the grayscale image to redirect measurements to (default = none)
  * <li>Selection of the outputs (Results table cleared is the default)
  * <li>Selection of the features to include in output if Results is checked (see notes below).
  * </ul><p>
  * <center><img src="BlobHandler-gui.png"></center>
  * <p>
  * <p>
  * <b><u>Feature selection</u></b> - many of these have better desciptions in the online ImageJ manual<br>
  * <b>BASIC</b>
  *   <ul>
  *      <li> Blob - the unique identifier of the blob(s), useful when combining blobs
  *      <li> Count - the number of individual blob(s) that make up this compound blob
  *      <li> Area - the area of the blob(s)
  *      <li> Perim - the total perimeter length of the blob(s)
  *   </ul>
  * <b>LOCATION</b>
  *   <ul>
  *      <li> BX, BY, BW and BH  The bounding box description as x0, y0, width and height
  *      <li> XStart, YStart The coordinates of the first pixel encountered in a scan.
  *      <li> Side A flag indicating which side(s) of the image the blob(s) touch(es), where <br>
  *           0 = none, 1 = left, 2 = top, 4 = right 8 = bottom or some additive combination of these
  *   </ul>
  * <b>SHAPE</b>
  *   <ul>
  *      <li> FeretMax The longest dimension found rotating a box at 2-degree increments.
  *      <li> FeretMin the shortest dimension ...
  *      <li> FeretMed the median long feret measurement for the 180 rotations
  *      <li> FeretX, FeretY and Theta use these to reconstruct a line parallel to FeretMax and
  *         passing through a point on the perimeter of the rotated blob(s).
  *      <li> Orient  The orientation of the object derived from mechanical moments.
  *      <li> Circ - circularity - (4*PI*Area)/Perim^2
  *      <li> Round - roundness - (4*Area)/(PI*FeretMax^2)
  *      <li> Solid - solidity - (AreaWithinOuterPerimeter)/(ConvexHullArea)
  *      <li> Compact - compactness - (SQRT(4*Area/PI))/(FeretMax)
  *      <li> Extent - (Area)/(BoundingBoxArea)
  *      <li> Aspect - (FeretMax)/(FeretMin)
  *      <li> VoidFraction - (AreaWithinInnerPerimeter)/(AreaWithinOuterPerimeter)
  *      <li> NEuler - Euler Shape Number - (NumberOuterPerimeters)-(NumberOfInnerPerimeters)
  *   </ul>
  * <b>VOLUME</b>
  *   Volume is useful for those working in the aquatic sciences.
  *   <ul>
  *   <li>We often use a 
  *   method to estimate volume from Sieracki, Viles and Webb (1989, see references). This
  *   method (SVW) has been shown to work well with relatively convex shapes.
  *   <li>To this method we add a simple ABD method which converts Area to an equivalent 
  *   radius (if the pixels were arranged as a circle) then computes volume as 4/3*PI*r^3
  *   <li>And we add a technique that combines the skeletonization with the Euclidian distance map (SKEDM). This
  *   is briefly described in Russ (2007, see refs) but not necessarily in the context of computing volume.
  *   We have found this method works well for objects with solidity lower than 0.4.
  *   <li>Finally, the appropriate method can be found using the AUTO method which uses SKEDM for solidity < 0.4
  *   and SVW in all other cases.
  *   <li> NTricky is reported with volume.  This is the count of blobs in the compound blob that qualified for special 
  *      treatment with the SKEDM algorithm - even if it was not used.
  *   </ul>
  *   
  * <b>ELLIPSE</b>
  *   <ul>
  *      <li> Major - the length of the semi-major axis
  *      <li> Minor - the length of the semi-minor axis
  *      <li> Ecc - eccentricity of the ellipse
  *      <li> Use Orient (above) for the orientation of the long axis
  *   </ul>
  * <b>GREYSCALE INFO</b> most meaningful if redirected to a greyscale image
  *   <ul>
  *      <li> GMin Minimum grey value
  *      <li> GMax Maximum grey value
  *      <li> GMean Mean grey value
  *      <li> GVar Greyscale variance (0 for binary images)
  *      <li> GSkew Greyscale skewness (NaN for binary images)
  *      <li> GKurt Greyscale kurtosis (NaN for binary images)
  *   </ul>
  * <b>MOMENT INVARIANTS</b> We suggest that the references are consulted for details
  *   <ul>
  *    <li> Hu moments are described in Burger and Berge, 2008
  *    <li> Affine Moment Invariants (AMI) are described in Flusser et al, 2009
  *   </ul>     
  * <p>
  * CONCEPTUAL FRAMEWORK:<br> 
  * A bi-level 2d image is comprised of foreground and background features. <br>
  * Contiguous groups of pixels are considered SingleBlobs <br>
  * CompoundBlobs are comprised of one or more SingleBlobs <br>
  * CompoundBlobs are reported in the results table.<br>
  * <p>
  * Contact: Mike Sieracki, msieracki@bigelow.org <br>
  * Author:  Ben Tupper, btupper@bigelow.org <br>
  * History: Developed over the course of 2009 and 2010 using ImageJ 1.43 and Mac OSX 10.5.8 <br>
  *   at Bigelow Laboratory for Ocean Science, www.bigelow.org
  */
  
public class BlobHandler_ implements PlugInFilter {
	private ImagePlus imp;
	private ImageProcessor ip;
	private ImagePlus gimp = null;
	private ImageProcessor gip = null;
	private BlobHandler bh = new BlobHandler();
	private String[] showNames = {"Results table", "Label image", "Clear Results", "Contour image", "Add to Roi Manager"};
	private int RESULTS = 0, LABEL = 1, CLEAR = 2, CONTOUR = 3, ADD = 4;
	private String[] showVolMethods = {"AUTO", "ABD", "SVW", "SKEDM", "None"};
	private int volMethod = BHBlob.VOL_AUTO;
	private boolean[] showThese = {true, false, true, false, true};
	private int measureFlag = BHBlob.BASIC;
	private double[] minMaxArea = {0.0, Double.POSITIVE_INFINITY};
	private double[] minMaxCirc = {0.0, 1.0};
	private boolean filterAfterCombine = false;

	private double combineIfCloserThan = 0.0;
	
	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		return DOES_8G;
	}

	public void run(ImageProcessor ip) {
	 
	  this.ip = ip;
    
	  String options = Macro.getOptions();
    
	  if ((options != null) && (options.length() != 0)){
	     if (parseOptions(options) == false) { return; } 
	  } else {
      if (showDialog() == false) { return; }
	  }
    
    if (gimp != null){
      gip = gimp.getProcessor();
    } else {
      gip = ip;
    }   
    if (process() == false) {
      IJ.showMessage("Process method failed!");
    }
       
	}//run

// the main work engine
  private boolean process(){
    
    Calibration cal = imp.getCalibration();
    //IJ.log("pH = " + cal.pixelHeight);
    //IJ.log("pW = " + cal.pixelWidth);
    
    if (bh.setup(ip, gip, imp.getCalibration(), measureFlag) == false){
        return false;
    }  
    bh.process(); 
     
    if (filterAfterCombine == false){
      bh.blobs.filterByAreaAndCirc(minMaxArea[0], minMaxArea[1], minMaxCirc[0], minMaxCirc[1]);
    }
    
    if ((bh.blobs.count() > 1) && (combineIfCloserThan > 0.0)) {
      int[] labels = bh.blobs.getLabels();
//      IJ.log("Distances:" + combineIfCloserThan);
//      double[][] dd = bh.blobs.computeMinimumDistanceTable();      
      
      int[][] table = bh.blobs.computeMergeTable(combineIfCloserThan);
      
      //IJ.showStatus("Combining nearest neighbors");
      for (int i = 0; i < table.length; i++){
       //IJ.showProgress(i, table.length-1);
       if (table[i].length > 1) {
          bh.blobs.combine(table[i]);
       } //table shows two or more to combine
      } //i-loop 
      //IJ.showProgress(2.0);
    }//combine?
    
    if (filterAfterCombine == true){
      bh.blobs.filterByAreaAndCirc(minMaxArea[0], minMaxArea[1], minMaxCirc[0], minMaxCirc[1]);
    }   
    
    if ((measureFlag & BHBlob.VOLUME) != 0) {
      //first set the volume measurement method
      bh.blobs.setVolumeMethod(volMethod);
      //then compute the volumes
      bh.blobs.computeVolume();
    }
   
  	 int[] labels = bh.blobs.getLabels();

    if (showThese[0] == true){
      ResultsTable rt = ResultsTable.getResultsTable();
    // make a results table to show results
      if (showThese[2] == true){  
  		  //rt = new ResultsTable();
  		  rt.reset();
      } else {
        //rt = ResultsTable.getResultsTable();
      }
  		BHCompoundBlob blob;
  		int row = rt.getCounter();
  		//int[] labels = bh.blobs.getLabels();
  		for (int i = 0; i < labels.length; i++){
  		  blob = bh.blobs.get(labels[i]);
  		  if (blob != null){
    	 	 	rt.incrementCounter(); 
    	 	 	blob.showInfo(rt, row); 
    	 	 	row++;
  		  } 
  		}
  		rt.updateResults();
  		rt.show("Results");
    }
    
    //show labeled image
    if (showThese[1] == true){
      ImageProcessor lip = new ShortProcessor(imp.getWidth(), imp.getHeight());
      bh.blobs.labelPixelsProcessor(lip);
      lip.resetMinAndMax();
      ImagePlus limp = new ImagePlus("Labeled Image", lip);
      limp.show();
     	WindowManager.setTempCurrentImage(limp);
     	
//     	Overlay limpOverlay = limp.getOverlay();
//      if (limpOverlay == null) { limpOverlay = new Overlay();}
//      limpOverlay = bh.blobs.showAllEllipses(limpOverlay);
//      limp.setOverlay(limpOverlay);
//      limp.updateAndDraw();
      
      IJ.run("3-3-2 RGB");
    }

    //show contours
    if (showThese[3] == true){
      ImagePlus theImp = (gimp != null) ? gimp : imp;
      Overlay overlay = theImp.getOverlay();
      if (overlay == null) {overlay = new Overlay();}
      overlay = bh.blobs.showAllContours(BHPolygon.OUTER + BHPolygon.CRACK, overlay);
      overlay = bh.blobs.showAllContours(BHPolygon.INNER + BHPolygon.CRACK, overlay);
      theImp.setOverlay(overlay);
      theImp.updateAndDraw();
   
      //ImageProcessor cip = new ShortProcessor(imp.getWidth(), imp.getHeight());
      //bh.blobs.labelContoursProcessor(cip, BHPolygon.OUTER);
      //bh.blobs.labelContoursProcessor(cip, BHPolygon.INNER);      
      //cip.resetMinAndMax();
      //ImagePlus cimp = new ImagePlus("Contoured Image", cip);
      //cimp.show();
      //WindowManager.setTempCurrentImage(cimp); 
      //IJ.run("3-3-2 RGB");

    } 
    return true;   
  }//process
  
  /**
    * This returns a list of the names of the image windows plus "None"
    */
  private String[] getImageList(){
    int n = WindowManager.getImageCount();
    if (n == 0){ return null;}
    
    int[] ids = WindowManager.getIDList();
    String[] out = new String[n+1];
    out[0] = "None";
    ImagePlus theImp = null;
    for (int i = 0; i < n; i++){
      theImp = WindowManager.getImage(ids[i]);
      out[i+1] = theImp.getTitle();
    }//i-loop
    return out; 
  }//getImageList
  
  /**
    * Returns the measureFlag as a series of booleans
    */
    private boolean[] getMeasureBoolean(){
       int[] opt = BHBlob.getMeasurementOptions();
       boolean[] b = new boolean[opt.length];
       for (int i = 0; i< opt.length; i++){
         b[i] = (measureFlag & opt[i]) != 0;
       }
      return b;
    }
    
    
/**
  * Parses the options passed from a macro - an example...
  * "area=0.0-Infinity circularity=0.0-1.0 filter combine=0.000 method=AUTO redirect=None results label clear contour basic location shape volume slice ellipse grayscale"
  */
  private boolean parseOptions(String options){
    boolean ok = true;
    for (int i = 0; i < showThese.length; i++){showThese[i] = false;}
    measureFlag = 0;
    String[] s = options.split(" ");
    String[] t;    
    for (int i = 0; i < s.length; i++){
      t = s[i].split("=");
      if (t[0].equalsIgnoreCase("area")){
        minMaxArea = parseRange(t[1], Double.POSITIVE_INFINITY);
      } else if (t[0].equalsIgnoreCase("circularity")){
        minMaxCirc = parseRange(t[1], 1.0);
      } else if (t[0].equalsIgnoreCase("combine")){
        combineIfCloserThan = Tools.parseDouble(t[1]);
      } else if (t[0].equalsIgnoreCase("method")){
         if (t[1].equalsIgnoreCase("Auto")){
            volMethod = BHBlob.VOL_AUTO;
         } else if (t[1].equalsIgnoreCase("SVW")){
            volMethod = BHBlob.VOL_SVW;
         } else if (t[1].equalsIgnoreCase("ABD")){
            volMethod = BHBlob.VOL_ABD;
         } else if (t[1].equalsIgnoreCase("SKEDM")){
            volMethod = BHBlob.VOL_SKEDM;
         } else {
            volMethod = BHBlob.VOL_NONE;
         }        
      } else if (t[0].equalsIgnoreCase("redirect")){
        if (!t[0].equalsIgnoreCase("none")){gimp = WindowManager.getImage(t[1]);}       
      } else if (t[0].equalsIgnoreCase("filter")){
        filterAfterCombine = true;
      } else if (t[0].equalsIgnoreCase("results")){
        showThese[RESULTS] = true;
      } else if (t[0].equalsIgnoreCase("clear")){
        showThese[CLEAR] = true;
      } else if (t[0].equalsIgnoreCase("label")){
        showThese[LABEL] = true;
      } else if (t[0].equalsIgnoreCase("contour")){
        showThese[CONTOUR] = true;
      }  else if (t[0].equalsIgnoreCase("add")){
        showThese[ADD] = true;      
      } else if (t[0].equalsIgnoreCase("basic")){
        measureFlag += BHBlob.BASIC;
      } else if (t[0].equalsIgnoreCase("location")){
        measureFlag += BHBlob.LOCATION;
      } else if (t[0].equalsIgnoreCase("shape")){
        measureFlag += BHBlob.SHAPE;
      } else if (t[0].equalsIgnoreCase("volume")){
        measureFlag += BHBlob.VOLUME;
      } else if (t[0].equalsIgnoreCase("gray")){  //gray and grayscale are the same
        measureFlag += BHBlob.GRAY;
      } else if (t[0].equalsIgnoreCase("grayscale")){
        measureFlag += BHBlob.GRAY;
      } else if (t[0].equalsIgnoreCase("slice")){
        measureFlag += BHBlob.SLICE;
      } else if (t[0].equalsIgnoreCase("ellipse")){
        measureFlag += BHBlob.ELLIPSE;
      } else if (t[0].equalsIgnoreCase("moments")){
        measureFlag += BHBlob.MOMENTS;
      } else if (t[0].equalsIgnoreCase("add")){
        showThese[ADD] = true;      
      } else {
        //parameter not known
      }
    }     
    
     
    return ok;
  }
  
  
  /**
    * Parses a string like "0-Infinity" into a double array
    */
  private double[] parseRange(String arg, double upperDefault){
    String[] t2 = arg.split("-");
    double[] d = new double[2];
    d[0] = Tools.parseDouble(t2[0]);
    if (Double.isNaN(d[0])) {d[0] = 0.0;} //make sure it is a real number
    d[1] = (t2.length >= 2) ? Tools.parseDouble(t2[1]) : upperDefault;
    if (Double.isNaN(d[1])) {d[1] = upperDefault;} //make sure it is a real number
    
    return d;
  }
  /**
    * Returns true of the user accepts and false if the user cancels.
    * Borrows heavily from ImageJ's ParticleAnalyzer.showDialog()
    */
  private boolean showDialog(){
    String[] names = getImageList();
    GenericDialog gd = new GenericDialog("Blob Handler Setup");
    gd.addStringField("Area range:" , minMaxArea[0] + "-" + minMaxArea[1]);
    gd.addStringField("Circularity range:" , minMaxCirc[0] + "-" + minMaxCirc[1]);
    gd.addCheckbox("Filter after combining blobs", filterAfterCombine );
    gd.addNumericField("Combine if closer than:", combineIfCloserThan, 3);
    gd.addChoice("Method for volume", showVolMethods, showVolMethods[0]);
    boolean[] b = getMeasureBoolean();
    int[] bvalues = BHBlob.getMeasurementOptions();
    String[] bnames = BHBlob.getMeasurmentOptionsNames(true);
    gd.addChoice("Redirect to:", names, names[0]);
    
    gd.addMessage("Select the results to show");
    gd.addCheckboxGroup(2,2, showNames, showThese);
    
    gd.addMessage("Select the measurements to include");
    int nrow = b.length/2 + (b.length % 2) ;
    gd.addCheckboxGroup(nrow, 2, bnames, b);
    
    gd.showDialog();
    
    if (gd.wasCanceled() == true){return false;}
    
   String redirectName = "";
   String text = "";
   String[] text2 = {"",""};
   
   //Area range
   text2 = Tools.split(gd.getNextString(),"-");
   minMaxArea[0] = Tools.parseDouble(text2[0]);
   if (Double.isNaN(minMaxArea[0])) { minMaxArea[0] = 0.0;}
   minMaxArea[1] = (text2.length==2) ? Tools.parseDouble(text2[1]) : Double.POSITIVE_INFINITY;
 
   text2 = Tools.split(gd.getNextString(),"-");
   minMaxCirc[0] = Tools.parseDouble(text2[0]);
   if (Double.isNaN(minMaxCirc[0])) { minMaxCirc[0] = 0.0;}
   minMaxCirc[1] = (text2.length==2) ? Tools.parseDouble(text2[1]) : 1.0;   
  
   filterAfterCombine = gd.getNextBoolean();
   
   combineIfCloserThan = gd.getNextNumber();
   if (Double.isNaN(combineIfCloserThan) == true) {combineIfCloserThan = 0.0;}
   
   String theName = gd.getNextChoice();
   if (theName.equalsIgnoreCase("Auto")){
      volMethod = BHBlob.VOL_AUTO;
   } else if (theName.equalsIgnoreCase("SVW")){
      volMethod = BHBlob.VOL_SVW;
   } else if (theName.equalsIgnoreCase("ABD")){
      volMethod = BHBlob.VOL_ABD;
   } else if (theName.equalsIgnoreCase("SKEDM")){
      volMethod = BHBlob.VOL_SKEDM;
   } else {
      volMethod = BHBlob.VOL_NONE;
   }
   //IJ.log("volMethod = " + volMethod);
   redirectName = gd.getNextChoice();
   gimp = WindowManager.getImage(redirectName);
   
   for (int i = 0; i< showThese.length; i++){
    showThese[i] = gd.getNextBoolean();
   }
   
   int theNewMeasure = 0;
   for (int i = 0; i < b.length; i++){
     if (gd.getNextBoolean() == true) { theNewMeasure += bvalues[i];}
   }
   measureFlag = theNewMeasure;
   
   return true;  
  }//showDialog

}
