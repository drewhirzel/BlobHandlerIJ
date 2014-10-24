import ij.*;
import java.util.Arrays;
/**
  * BHPointGroup is a simple container for an ensemble of BHPoint objects.
  * This is intended to act like Vector, but simplify the 
  * getting as an array - <Vector>.toArray() seems like overkill
  * besides - this now makes the collection of points similar to
  * the contours.  Original label associations can be found in the 
  * label attribute of each point.
  *
  * <p>
  * Contact: Mike Sieracki, msieracki@bigelow.org <br>
  * Author:  Ben Tupper, btupper@bigelow.org <br>
  * History: Developed over the course of 2009 and 2010 using ImageJ 1.43 and Mac OSX 10.5.8 <br>
  *   at Bigelow Laboratory for Ocean Science, www.bigelow.org
  */



public class BHPointGroup {
  /** the number of points */
  public int npoints = 0;
  /** the array of points */
  public BHPoint[] points;
  /** The size of the array increase - changed internally as needed */
  private int chunkSize = 100;
  
/**
  * Generic constructor
  */
  public BHPointGroup(){
   initialize(); 
  } 
/**
  * Constructor with the first point to add
  */
  public BHPointGroup(BHPoint point){
    initialize();
    add(point);
  }
  
/**
  * Initialize the points
  */
  private void initialize(){
      points = new BHPoint[chunkSize];
      npoints = 0;
  }//initialize

/**
  * Adds the points from one group to this group
  */
  public int add(BHPointGroup otherPoints){
    for (int i = 0; i<otherPoints.npoints; i++){
        add(otherPoints.get(i));
    }
    return npoints;
  }
  
/**
  * Adds an array of points
  */
  public int add(BHPoint[] newPoints){
    for (int i = 0; i < newPoints.length; i++){
      add(newPoints[i]);
    } 
    return npoints;
  }
/**
  * Adds a point
  */
  public int add(BHPoint point){
    if (npoints == points.length){
      BHPoint[] temp = new BHPoint[chunkSize*2];
      System.arraycopy(points, 0, temp, 0, npoints);
      points = temp;
      chunkSize *= 2;
    }
    points[npoints] = point;
    npoints++;
    return npoints;
  }  
  
  //trims the size of the point array to npoints
  public int trim(){
    if (npoints < points.length){
      BHPoint[] temp = new BHPoint[npoints];
      System.arraycopy(points, 0, temp, 0, npoints);
      points = temp;
    }
    npoints = points.length;
    return npoints;
  }//trim

/**
  * Returns the specified point or null if the index is out of bounds
  * @param index the zero-based index of the point to retrieve
  * @return the BHPoint requested or null if the index is bad
  */
  public BHPoint get(int index){
   
   if ((index < 0) || (index >= npoints)) {
      return null;
   } else {
      return points[index];
   }
    
  }//get
  
/**
  * Returns the entire array of points - the array is trimmed first
  */
  public BHPoint[] getArray(){
    trim();
    return points;
  }


/**
  * Returns the points as an array 
  */
  public BHPoint[] getPointsAsArray(){
    BHPoint[] arr =  new BHPoint[npoints];
    System.arraycopy(points, 0, arr, 0, npoints);
   return arr;    
  }
    
/**
  * Returns the upper-left point (which might not be the first in the list)
  */
  public BHPoint getStart(){
    if (npoints == 0) {return null;}
    
    int ix = 0;
    for (int i = 0; i < npoints; i++){
      if ((points[i].x < points[ix].x) && (points[i].y < points[ix].y)){
        ix = i;
      }
    }
    return points[ix];
  }
  
/**
  * Returns the entire array of points - a wrapper around getArray()
  */
  public BHPoint[] get(){
    return getArray();
  }
/**
  * Returns the size of the collection
  */
  public int size(){
    return npoints;
  }
/**
  * Clears all the data
  */
  public void clear(){
    initialize();
  }
  
  
/**
  * Translate the points by [dx,dy]
  * @param dx the amount to translate in x
  * @param dy the amount tp translate in y
  */
  public void translate(int dx, int dy){
    if (npoints > 0) {
      for (int i = 0; i<npoints; i++){
        points[i].translate(dx, dy);
      }
    } 
  }
  
  
/**
  * Computes the bounding rectangular coordinates
  * @return a 4 element array of [xmin,xmax, ymin,ymax]
  */
  public int[] getBounds(){
   
   int[] b = {Integer.MAX_VALUE, Integer.MIN_VALUE, Integer.MAX_VALUE, Integer.MIN_VALUE};
   for (int i = 0; i < npoints; i++){
    if (points[i].x < b[0]) { b[0] = points[i].x;}
    if (points[i].x > b[1]) { b[1] = points[i].x;}
    if (points[i].y < b[2]) { b[2] = points[i].y;}
    if (points[i].y > b[3]) { b[3] = points[i].y;}    
   } 
   return b;
  }//getBounds
/**
  * print the points
  */
  public void printPoints(){
   if (npoints != 0){
     IJ.log("i = [x,y, label, value]");
     for (int i = 0; i < npoints; i++){
      IJ.log("i=" + i + ": " + points[i].x + ", " + points[i].y + ", " + points[i].label + ", " + points[i].val);
     }
   }
  }
  
/**
  * Duplicates the properties of this group and returns a new reference
  */
  public BHPointGroup duplicate(){
    BHPointGroup newpg = new BHPointGroup();
    for (int i = 0; i < npoints; i++){
      newpg.add(points[i].duplicate());
    }
    return newpg; 
  }
  
  
 /**
   * sorts the points into x then by y and no returns a sorted copy
   */
   public BHPoint[] getSortedPoints() {    
      BHPoint[] pg = duplicate().getArray();
      Arrays.sort(pg);
      return pg;
   }

 /**
   * sorts the points into x then by y and no returns a sorted copy
   */
   public BHPointGroup getSortedPointGroup() {    
      BHPoint[] p = getSortedPoints();
      BHPointGroup pg = new BHPointGroup();
      pg.add(p);
      return pg;
   } 
 
 
/*
 * Gets the points in a n x 2 array [x][y] order
 */
 public int[][] getXY(){
   return getXY(false);
 } //getXY
  
 /**
  * Returns the x points
  */
  public int[] getX(){
//    int[] x = new int[npoints];
//    for (int i = 0; i < npoints; i++){x[i] = points[i].x;}
    return getX(true);
  }
  
/**
  * Returns the y points
  */
  public int[] getY(){
//    int[] y = new int[npoints];
//    for (int i = 0; i < npoints; i++){y[i] = points[i].y;}
    return getY(true);
  }  
  

/**
  * gets the X points with optional closure
  */
  public int[] getX(boolean closeIt){

    if (npoints == 0) {return null;}
    int n = npoints;
    if (closeIt == true) { n++; }  
    int[] x = new int[n];
    
    for (int i = 0; i < n; i++){
        x[i] = points[i % npoints].x;
    }
    return x;  
  }

/**
  * gets the Y points with optional closure
  */
  public int[] getY(boolean closeIt){

    if (npoints == 0) {return null;}
    int n = npoints;
    if (closeIt == true) { n++; }  
    int[] y = new int[n];
    
    for (int i = 0; i < n; i++){
        y[i] = points[i % npoints].y;
    }
    return y;  
  }



 
 /*
  * Gets the points in a n x 2 array [x][y] order with the option
  * to close the path (first and last point the same)
  */ 
 public int[][] getXY(boolean closeIt){
    if (npoints == 0) {return null;}
    int n = npoints;
    if (closeIt == true) { n++; }
    int[][] xy = new int[n][2];
    
    for (int i = 0; i < n; i++){
        xy[i][0] = points[i % npoints].x;
        xy[i][1] = points[i % npoints].y;
    }
    return xy;
 } // getXY


 /*
  * Gets the points in a n x 2 array [x][y] order with the option
  * to close the path (first and last point the same) as double.
  * No scaling is performed
  */ 
 public double[][] getXYDoubleRaw(boolean closeIt){
    if (npoints == 0) {return null;}
    int n = npoints;
    if (closeIt == true) { n++; }
    double[][] xy = new double[n][2];
    for (int i = 0; i< n; i++){
      xy[i][0] = (double) points[i % npoints].x;
      xy[i][1] = (double) points[i % npoints].y;      
    }
    return xy;
 } // getXY
 
// /*
//  * Gets the points in a n x 2 array [x][y] order with the option
//  * to close the path (first and last point the same) as double
//  */ 
// public double[][] getXYDouble(boolean closeIt){
//    if (npoints == 0) {return null;}
//    int n = npoints;
//    if (closeIt == true) { n++; }
//    double[][] xy = new double[n][2];
//    for (int i = 0; i< n; i++){
//      xy[i][0] = (double) points[i % npoints].x * xscale[0] + xscale[1];
//      xy[i][1] = (double) points[i % npoints].y * yscale[0] + yscale[1];      
//    }
//    return xy;
// } // getXY
 
/**
  * Returns a two element array of [xc, yc]
  */
  public double[] getCentroid(){
   
   int xsum=0;
   int ysum=0;
   int zsum=0;
   for (int i = 0; i < npoints; i++){
      xsum += points[i].x;
      ysum += points[i].y;
      zsum += points[i].z;
   } 
    
   double[] c = new double[3];
   if (npoints > 0){
     c[0] = xsum / (double) npoints;
     c[1] = ysum / (double) npoints;
     c[3] = zsum / (double) npoints;
   } 
   return c; 
  } //getCentroid
  
  
/**
  * Returns the centroid distance r(t) = ([x(t) – xc]^2+ [y(t) - yc]^2)^1/2
  * for each point on the polygon
  */
  public double[] getCentroidDistance( double xc, double yc){
   
    double[] r = new double[npoints]; 
    
    for (int i = 0; i < npoints; i++){
      r[i] = Math.sqrt(Math.pow(points[i].x - xc,2.0) + Math.pow(points[i].y-yc,2.0));
    }
    return r;
  }//getCentroidDistance
  public double[] getCentroidDistance(){
    double[] xy = getCentroid();
    return getCentroidDistance(xy[0], xy[1]);
  }

  
///**
//  * Returns a 4 integer array of [xmin, xmax, ymin, ymax] of the bounds
//  * of this polygon
//  */
//  public int[] getBounds(){
//       
//    int[] r = {points[0].x, points[0].x, points[0].y, points[0].y};  
//    for (int i = 0; i < npoints; i++){
//      if (points[i].x < r[0]) { r[0] = points[i].x;}
//      if (points[i].x > r[1]) { r[1] = points[i].x;}
//      if (points[i].y < r[2]) { r[2] = points[i].y;}
//      if (points[i].y > r[3]) { r[3] = points[i].y;}
//    } //i loop
//    return r;
//  }//getBounds

///**
//  * rotates the points about the specified centroid by the specified angle
//  * See Horn 1986.  Instead of rotating this one, it simply returns 
//  * the rotated coords as [nvertices][2] - this type of contour only handles
//  * integer type coords so it isn't suitable for rotation storage
//  */
//  public double[][] getRotated(double xc, double yc, double angle){
//    double[][] xy = getXYDouble(false);
//    double sin = Math.sin(angle);
//    double cos = Math.cos(angle);
//    double dx = 0.0;
//    double dy = 0.0;
//    for (int i = 0; i < xy.length; i++){
//      dx = xy[i][0] - xc;
//      dy = xy[i][1] - yc;
//      xy[i][0] = xc + (dx * cos) + (dy * sin);
//      xy[i][1] = yc - (dx * sin) + (dy * cos);
//    }
//    
//    return xy;
//  }
  
/**
  * rotates the points about the specified centroid by the specified angle
  * See Horn 1986.  Instead of rotating this polygon, it simply returns 
  * the rotated coords as [nvertices][2] - this type of contour only handles
  * integer type coords even though the returned coords are double
  * Coordinates are untransformed by scale - thus they are "raw" -
  * be  sure that xc and yc inputs are also untransformed
  * @param xc  The unscaled x coordinate of the centroid
  * @param yc The unscaled y coordinate of the centroid
  * @param angle the angle in radians about which to rotate
  */
  public double[][] getRotatedRaw(double xc, double yc, double angle){
    double[][] xy = getXYDoubleRaw(false);
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
  }//getRotatedRaw
      
}//BHPointGroup defintion