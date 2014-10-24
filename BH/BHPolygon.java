import ij.*;
import ij.gui.*;
import ij.process.*;
import ij.measure.*;
import java.awt.Color;
import java.awt.geom.*;
/**
  * This is a simple container wrapper around BHPointGroup specifically for
  * contour points, convex hulls and ellipses.  Contour points are retained for the 
  * pixel locations (upper left of perimeter pixels) AND the crack locations
  * (interpixel like ImageJ's Wand). 
  * <p>
  * Two items of special interest are the convex hull computation
  * which is exactly like ImageJ's and the minimum distance computation for 
  * computing the distance between two polygons. 
  *
  * <p>
  * Contact: Mike Sieracki, msieracki@bigelow.org <br>
  * Author:  Ben Tupper, btupper@bigelow.org <br>
  * History: Developed over the course of 2009 and 2010 using ImageJ 1.43 and Mac OSX 10.5.8 <br>
  *   at Bigelow Laboratory for Ocean Science, www.bigelow.org <br>
  * 2010-03-05 BTT slightly modified the descriptors to include CRACK and changed <br>
  *  ocont and icont to opixel and ipixel to make it clear;
  *
  */


public class BHPolygon extends BHPointGroup{
  
   /** Polygon type is coded by OUTER, INNER, CHULL or ELLIPSE */
   public static final int OUTER = 1; //boundary may be pixel or crack
   public static final int INNER = 2; //boundary may be pixel or crack
   public static final int CRACK = 4;
   public static final int CHULL = 8; 
   public static final int ELLIPSE = 16;
   
   /** The default type is OUTER */
   public int type = OUTER;
   
   /** This is the blob label associated with this contour */
   public int label = 0;
   public Calibration cal;

 /*
  * Construct as empty container
  */
 public BHPolygon(){
  setScale();
 }
 
 /*
  * Construct as empty container
  */
 public BHPolygon(Calibration cal){
  setScale(cal);
 }
 
 /*
  * Construct as empty container
  */
 public BHPolygon(int theType,  int theLabel){
  type = theType;
  label = theLabel;
  setScale();
 }
 
 /*
  * Construct as empty container
  */
 public BHPolygon(int theType, int theLabel, Calibration cal){
  type = theType;
  label = theLabel;
  setScale(cal);
 }
  
 /*
  * Construct with the first seed point
  */
  
 public BHPolygon(int theType, int theLabel, int x, int y){
  type = theType;
  label = theLabel;
  setScale();
  add(x,y);
 }

 /*
  * Construct with the first seed point
  */
  
 public BHPolygon(int theType, int theLabel, int x, int y, Calibration cal){
  type = theType;
  label = theLabel;
  setScale(cal);
  add(x,y);
 }
  
 public BHPolygon(int[] x, int[] y, int theType, int theLabel){
    type = theType;
    label = theLabel;
    setScale();
    for (int i = 0; i < x.length; i++){add(x[i],y[i]);}
 } 

 public BHPolygon(int[] x, int[] y, int theType, int theLabel, Calibration cal){
    type = theType;
    label = theLabel;
    setScale(cal);
    for (int i = 0; i < x.length; i++){ add(x[i],y[i]); }
 }
 
 public BHPolygon(int theType, int theLabel, BHPoint thePoint){
    type = theType;
    label = theLabel;
    setScale();
    add(thePoint);
 }
 
  public BHPolygon(int theType, int theLabel, BHPoint thePoint, Calibration cal){
    type = theType;
    label = theLabel;
    setScale(cal);
    add(thePoint);
 }
 
 public BHPolygon(int theType, int theLabel, BHPoint[] thePoints){
    type = theType;
    label = theLabel;
    setScale();
    add(thePoints);
 } 

 public BHPolygon(int theType, int theLabel, BHPoint[] thePoints, Calibration cal){
    type = theType;
    label = theLabel;
    setScale(cal);
    add(thePoints);
 } 
 
 
/**
  * Returns TRUE if this is an outer perimeter
  */
  public boolean isOuter(){
      return isOuter(type);
  }
  
  public static boolean isOuter(int typeCode){
     return (typeCode & OUTER) != 0;
  }

/***
    Returns  TRUE if this is an inner contour
    */
   public boolean isInner(){
      return isInner(type);
  }
  
  public static boolean isInner(int typeCode){
     return (typeCode & INNER) != 0;
  }    
  
/**
  * Returns true if this is a crack boundary
  */
   public boolean isCrack(){
      return isCrack(type);
  }
  
  public static boolean isCrack(int typeCode){
     return (typeCode & CRACK) != 0;
  }
    
/***
   * Returns TRUE if this is a convex hull
   */
   public boolean isCHull(){
      return isCHull(type);
  }
  
  public static boolean isCHull(int typeCode){
     return (typeCode & CHULL) != 0;
  }   

/***
   * Returns TRUE if this is an ellipse
   */
   public boolean isEllipse(){
      return isEllipse(type);
  }
  
  public static boolean isEllipse(int typeCode){
     return (typeCode & CHULL) != 0;
  }   
  

/**
  * Returns the color for this object
  */
  public Color getColor(){
   return getColor(type);
  }
/**
  * Returns the Color for the specified type of type of object<br>
  * OUTER = green (pixel or crack)<br>
  * INNER = blue (pixel or crack)<br>
  * CHULL = yellow<br>
  * ELLIPSE = orange<br>
  */
  public static Color getColor(int which){
  
   Color color = null;
   switch(which){
      case BHPolygon.OUTER: color = Color.GREEN; break;
      case BHPolygon.INNER: color = Color.BLUE; break;
      case BHPolygon.CHULL: color = Color.YELLOW; break;
      case BHPolygon.ELLIPSE: color = Color.ORANGE; break;
      case (BHPolygon.OUTER + BHPolygon.CRACK): color = Color.GREEN; break;
      case (BHPolygon.INNER + BHPolygon.CRACK): color = Color.BLUE; break;
      default: color = Color.CYAN;
   } 
  
  return color;
  }
/**
  * Sets the spatial scaling with default calibration
  */
  public void setScale(){
    setScale(new Calibration());
  }

/**
  * Sets the spatial scaling derived from the input Calibration
  */
  public void setScale(Calibration cal){
   this.cal = cal;
  }
    
/**
  * Adds a point - assigns the point to this polygon
  */
  public void add(int x, int y){
    BHPoint point = new BHPoint(x, y, label);
    add(point); 
  } //add


/**
  * Returns the length just like computed by ImageJ's
  * PolygonRoi.getTracedPerimeter
  */
  public double getTracedLength(){
      
      if (npoints == 1) { return 0.0 ;}
      int[][] xy= getXY(false);  
      int sumdx = 0;
      int sumdy = 0;
      int nCorners = 0;
      double dx1 = xy[0][0] - xy[npoints-1][0];
      double dy1 = xy[0][1] - xy[npoints-1][1];
      int side1 = (int) Math.abs(dx1) + (int) Math.abs(dy1); //one of these is 0
      boolean corner = false;
      int nexti, dx2, dy2, side2;
      for (int i=0; i<npoints; i++) {
          nexti = i+1;
          if (nexti==npoints)
            nexti = 0;
          dx2 = xy[nexti][0] - xy[i][0];
          dy2 = xy[nexti][1] - xy[i][1];
          sumdx += Math.abs(dx1);
          sumdy += Math.abs(dy1);
          side2 = (int) Math.abs(dx2) + (int) Math.abs(dy2);
          if (side1>1 || !corner) {
            corner = true;
            nCorners++;
          } else
            corner = false;
          dx1 = dx2;
          dy1 = dy2;
          side1 = side2;
      }
      double w = cal.pixelWidth;
      double h = cal.pixelHeight;
      
      return sumdx*w+sumdy*h-(nCorners*((w+h)-Math.sqrt(w*w+h*h)));
  }
  
/**
  * Returns the length of this polygon - this is the naive method
  * - to match ImageJ's use getTracedLength
  */
  public double getLength(){
   if (npoints == 1) { return 0.0 ;}
   double[][] xy= getXYDouble(true);
   double len = 0.0;
   double dx = 0.0;
   double dy = 0.0;
   for (int i = 0; i < (xy.length-1); i++){
     dx = (xy[i+1][0] - xy[i][0]);
     dy = (xy[i+1][1] - xy[i][1]);
     len += Math.sqrt(dx*dx + dy*dy);
   }
   return len;
  }//getLength  
  
/**
  * Computes the area enclosed by a contour - assumes a simple polygon (non-self-intersecting)
  * See http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/
  * and ImageJ's Roi.getArea()
  * area = 0.5*ABS(sum(x[i]*y[i+1]-x[i+1]*y[i])) for i = 0...npoints
  */
  public double computeArea(){
    double area = 0.0;
    //if (npoints == 1) { return xscale[0] * yscale[0];}
    if (npoints > 2){
      double[][] xy = getXYDouble(true); //slight of hand - has first point attached to end
      for (int i = 0; i < npoints; i++){ //so we can use npoints+1
        area += xy[i][0]*xy[i+1][1] - xy[i+1][0]*xy[i][1];
      }
      area = 0.5 * Math.abs(area);
    } else {
      area = cal.pixelWidth * cal.pixelHeight;
      //area = xscale[0] * yscale[0];
    }
    return area;
  }//computeArea
  
/**
  * rotates the points about the specified centroid by the specified angle
  * See Horn 1986.  Instead of rotating this one, it simply returns 
  * the rotated coords as [nvertices][2] - this type of contour only handles
  * integer type coords so it isn't suitable for rotation storage
  */
  public double[][] getRotated(double xc, double yc, double angle){
    double[][] xy = getXYDouble(false);
    double sin = Math.sin(angle);
    double cos = Math.cos(angle);
    double dx = 0.0;
    double dy = 0.0;
    for (int i = 0; i < xy.length; i++){
      dx = xy[i][0] - xc;
      dy = xy[i][1] - yc;
      xy[i][0] = xc + (dx * cos) + (dy * sin);
      xy[i][1] = yc - (dx * sin) + (dy * cos);
    }
    
    return xy;
  }
  
  
/**
  * Returns a PolygonRoi see ij.gui.PolygonRoi
  */
  public PolygonRoi getPolygon(){
   PolygonRoi poly = new PolygonRoi(getX(false), getY(false), npoints, Roi.POLYGON); 
   poly.setStrokeColor(getColor());
   return poly;
  } //getPolygon

 /*
  * Gets the points in a n x 2 array [x][y] order with the option
  * to close the path (first and last point the same) as double
  */ 
 public double[][] getXYDouble(boolean closeIt){
    if (npoints == 0) {return null;}
    int n = npoints;
    if (closeIt == true) { n++; }
    double[][] xy = new double[n][2];
    for (int i = 0; i< n; i++){
      xy[i][0] = cal.getX(points[i % npoints].x);
      xy[i][1] = cal.getY(points[i % npoints].y);  
      //xy[i][0] = (double) points[i % npoints].x * xscale[0] + xscale[1];
      //xy[i][1] = (double) points[i % npoints].y * yscale[0] + yscale[1];      
    }
    return xy;
 } // getXY

  
/**
  * Returns the convex hull of this contour
  */
   public BHPolygon getConvexHull(){
    return getConvexHull(this);
   }//getConvexHull
   
/**
  * Returns the convex hull of the specified contour
  * This is adapted (only slightly) from Wayne Rasband's 
  * gift-wrapping implementation in ij.gui.PolygonRoi class
  */
   public BHPolygon getConvexHull(BHPolygon c) {
      return getConvexHull(c.getX(false), c.getY(false));
    }//getConvexHull
 
 
 /**
  * Returns the convex hull of the specified contour
  * This is adapted (only slightly) from Wayne Rasband's 
  * gift-wrapping implementation in ij.gui.PolygonRoi class
  */
 public BHPolygon getConvexHull(int[] xCoordinates, int[] yCoordinates){
      //printPoints();
      int n = xCoordinates.length;
      if (yCoordinates.length != n) {
         IJ.log("BHPolygon.getConvexHull inputs x and y must be same length");
         return null;
      }
      if (n == 1) {return new BHPolygon(xCoordinates, yCoordinates, CHULL, label);}
      
        
//        for (int i = 0; i < xCoordinates.length; i++){
//         IJ.log("xC,yC " + label + ": " + i + " = [" + xCoordinates[i] + ", " + yCoordinates[i] + "]");
//        }

      int xbase = 0;
      int ybase = 0;
      
      int[] xx = new int[n];
      int[] yy = new int[n];
      int n2 = 0;
      int smallestY = Integer.MAX_VALUE;
      int x, y;
      for (int i=0; i<n; i++) {
          y = yCoordinates[i];
          if (y<smallestY)
          smallestY = y;
      }
      int smallestX = Integer.MAX_VALUE;
      int p1 = 0;
      for (int i=0; i<n; i++) {
          x = xCoordinates[i];
          y = yCoordinates[i];
          if (y==smallestY && x<smallestX) {
              smallestX = x;
              p1 = i;
          }
      }
      int pstart = p1;
      int x1, y1, x2, y2, x3, y3, p2, p3;
      int determinate;
      int count = 0;
      do {
          x1 = xCoordinates[p1];
          y1 = yCoordinates[p1];
          p2 = p1+1; if (p2==n) p2=0;
          x2 = xCoordinates[p2];
          y2 = yCoordinates[p2];
          p3 = p2+1; if (p3==n) p3=0;
          do {
              x3 = xCoordinates[p3];
              y3 = yCoordinates[p3];
              determinate = x1*(y2-y3)-y1*(x2-x3)+(y3*x2-y2*x3);
              if (determinate>0)
                  {x2=x3; y2=y3; p2=p3;}
              p3 += 1;
              if (p3==n) p3 = 0;
          } while (p3!=p1);
          if (n2<n) { 
              xx[n2] = xbase + x1;
              yy[n2] = ybase + y1;
              n2++;
          } else {
              count++;
              if (count>10) return null;
          }
          p1 = p2;
      } while (p1!=pstart);

      //IJ.log("n, n2 " + n + ", " + n2 + ", " + (n2 < n));
      if (n2 < n) {
        int[] tx = new int[n2];
        int[] ty = new int[n2];
        System.arraycopy(xx, 0, tx, 0, n2);
        System.arraycopy(yy, 0, ty, 0, n2);
        xx = tx;
        yy = ty;
        
//        for (int i = 0; i < xx.length; i++){
//         IJ.log("CHULL " + label + ": " + i + " = [" + xx[i] + ", " + yy[i] + "]");
//        }
      }
      return new BHPolygon(xx, yy, CHULL, label);
 }//getConvexHull


// /** ORIGINAL
//  * Returns the convex hull of the specified contour
//  * This is adapted (only slightly) from Wayne Rasband's 
//  * gift-wrapping implementation in ij.gui.PolygonRoi class
//  */
// public BHPolygon getConvexHull(int[] inX, int[] inY){
//    
//      int n = inX.length;
//      if (inY.length != n) {
//         IJ.log("BHPolygon.getConvexHull inputs x and y must be same length");
//         return null;
//      }
//      if (n == 1) {return new BHPolygon(inX, inY, CHULL, label);}
//      
//      int[] r = getBounds();//[xmin,xmax, ymin,ymax]
//      IJ.log("r=" + r[0]+ ", "+ r[1]+ ", " + r[2]+ ", " + r[3]);
//      for (int i=0;i<inX.length;i++){
//         IJ.log(" " + inX[i] + ", " + inY[i]);}
//      
//      int nmax = n/2;
//      if (nmax < 10) {nmax = 10;}
//      
//      int[] xx = new int[n];
//      int[] yy = new int[n];
//      int n2 = 0;
//      int smallestY = Integer.MAX_VALUE;
//      int x = 0;
//      int y = 0;
//      for (int i=0; i<n; i++) {
//          y = inY[i];
//          if (y<smallestY)
//          smallestY = y;
//      }
//      int smallestX = Integer.MAX_VALUE;
//      int p1 = 0;
//      for (int i=0; i<n; i++) {
//          x = inX[i];
//          y = inY[i];
//          if (y==smallestY && x<smallestX) {
//              smallestX = x;
//              p1 = i;
//          }
//      }
//      
//      
//      int pstart = p1;
//      int x1, y1, x2, y2, x3, y3, p2, p3;
//      int determinate;
//      int count = 0;
//      
//      IJ.log("p1, x,y = " + p1 + ", "+ x + ", " + y);
//      
//      do {
//          x1 = inX[p1];
//          y1 = inY[p1];
//          p2 = p1+1; if (p2==n) p2=0;
//          x2 = inX[p2];
//          y2 = inY[p2];
//          p3 = p2+1; if (p3==n) p3=0;
//          do {
//              x3 = inX[p3];
//              y3 = inY[p3];
//              determinate = x1*(y2-y3)-y1*(x2-x3)+(y3*x2-y2*x3);
//              if (determinate>0)
//                  {x2=x3; y2=y3; p2=p3;}
//              p3 += 1;
//              if (p3==n) p3 = 0;
//          } while (p3!= p1);
//          if (n2<n) { 
//              IJ.log("n2, x1, y1 = " + n2 + ", " + x1 + ", " + y1);
//              if ( (x1 == inX[pstart]) && (y1 == inY[pstart])){
//               p2 = p1;
//              } else {
//                 xx[n2] = x1;
//                 yy[n2] = y1;
//                 n2++;
//              }
//          } else {
//              count++;
//              if (count>nmax) {
//               IJ.log("BHPolygon.getConvexHull:"+label+" n2 >= n more than "+nmax+" times, n=" + n);
//               return null;
//            }
//          }
//          p1 = p2;
//      } while (p1!=pstart);
//      
//      if (n2 < n) {
//        int[] tx = new int[n2];
//        int[] ty = new int[n2];
//        System.arraycopy(xx, 0, tx, 0, n2);
//        System.arraycopy(yy, 0, ty, 0, n2);
//        xx = tx;
//        yy = ty;
//      }
//      return new BHPolygon(xx, yy, CHULL, label);
// }//getConvexHull

/**
  * Labels the pixels in the IP with the specifed color and offsets from the 
  * origin
  */
  public void labelContour(ImageProcessor ip, int col, int xoff, int yoff){
      if (npoints > 0){
        for (int i = 0; i < npoints; i++){
          ip.putPixel(points[i].x + xoff, points[i].y + yoff, col); 
        }
      }
  } //labelContour
  
/**
  * Labels the pixels in the IP with the specifed color
  */
  public void labelContour(ImageProcessor ip, int col){
    labelContour(ip, col, 0, 0);
  } //labelContour
  
/**
  * Labels the pixels in the IP with the label color
  */
  public void labelContour(ImageProcessor ip){
    labelContour(ip, label);
  } //labelContour

/**
  * Prints the polygon to the log
  */
 public void printPoints(){
  String typename =  (type == OUTER) ? "outer" : "inner";
  IJ.log("[X,Y] coordinates for "+ typename + " contour for label = " + label);
  for (int i = 0; i< npoints; i++){
    IJ.log(points[i].x + ", " + points[i].y);
  } 
 }


/***
   * Compute the shortest line connecting two object
   **/
   public Line computeShortestLineTo(BHPolygon poly){
   
      double d = Double.MAX_VALUE;
      double temp = 0.0;
      double tempx = 0.0;
      double tempy = 0.0;
      double[][] p1 = getXYDouble(false);
      double[][] p2 = poly.getXYDouble(false);
      int ip1 = 0;
      int ip2 = 0;
      
      for (int i = 0; i < npoints; i++){
         for (int j = 0; j < poly.npoints; j++){
            tempx = (p1[i][0] - p2[j][0]);
            tempy = (p1[i][1] - p2[j][1]);
            temp = tempx*tempx + tempy*tempy;
            if (temp < d) { 
               d = temp;
               ip1 = i;
               ip2 = j;
            }
         } //j-loop
      } // i-loop

      Line line = new Line(p1[ip1][0], p1[ip1][1], p2[ip2][0],p2[ip2][1]);
      return line; 
   } 
/**
  * Computes the minimum distance between this polygon and the 
  * provided polygon - using pythagorean quadratic 
  *  sqrt(dx^2 + dy^2)
  */
  public double computeMinimumDistanceTo(BHPolygon poly){
    double d = Double.MAX_VALUE;
    double temp = 0.0;
    double tempx = 0.0;
    double tempy = 0.0;
    double[][] p1 = getXYDouble(false);
    double[][] p2 = poly.getXYDouble(false);
    
    for (int i = 0; i < npoints; i++){
      for (int j = 0; j < poly.npoints; j++){
        tempx = (p1[i][0] - p2[j][0]);
        tempy = (p1[i][1] - p2[j][1]);
        temp = tempx*tempx + tempy*tempy;
        if (temp < d) { d = temp;}
      } //j-loop
    } // i-loop
    
    return Math.sqrt(d);
  }// minimumDistanceBetween
    
/**
  * Duplicates the properties of this polygon and returns a new reference
  */
  public BHPolygon duplicate(){
    BHPoint[] newpoints = new BHPoint[npoints];
    for (int i = 0; i < npoints; i++){
      newpoints[i] = points[i].duplicate();
    }
    BHPolygon newp = new BHPolygon(type, label, newpoints, cal);
    return newp;
  }//duplicate
  
  
/***
   * Returns a java.awt.geom.GeneralPath
   */
   public GeneralPath getGeneralPath(){
   
      GeneralPath path = new GeneralPath();
      double[][] xy = getXYDouble(true);
      for (int i = 0; i < xy.length; i++){
         path.lineTo(xy[i][0], xy[i][1]);
      }
      return path;
   }
}