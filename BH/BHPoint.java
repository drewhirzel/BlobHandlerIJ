
/**
  * A generic point handler for BHBlob - points are handled in pixel units.
  * A key behavior of this class is that it implements Comparable which
  * allows for rapid sorting based upon x and y (and possibly z later)
  *
  * Contact: Mike Sieracki, msieracki@bigelow.org <br>
  * Author:  Ben Tupper, btupper@bigelow.org <br>
  * History: Developed over the course of 2009 and 2010 using ImageJ 1.43 and Mac OSX 10.5.8 <br>
  *   at Bigelow Laboratory for Ocean Science, www.bigelow.org
  */
  public class BHPoint implements Comparable {
   /** The pixel's x coordinate in pixel units */
   public int x;
   /** The pixel's y coordinate in pixel units */
   public int y;
   public int z; //unused
   /** The pixel's grouping label */
   public int label; //is this ever used?
   /** The pixel's grayscale value */
   public double val;// grayscale value
   
   public BHPoint(int theX, int theY){
    x = theX;
    y = theY;
   } //constructor
    
   public BHPoint(int theX, int theY, int theLabel){
    x = theX;
    y = theY;
    label = theLabel;
   } //constructor    

   public BHPoint(int theX, int theY, int theLabel, double theValue){
    x = theX;
    y = theY;
    label = theLabel;
    val = theValue;
   } //constructor    
    
   public BHPoint(int theX, int theY, int theLabel, float theValue){
    x = theX;
    y = theY;
    label = theLabel;
    val = (double) theValue;
   } //constructor  
   
           
/**
  * Move the point by [dx, dy]
  */
   public void translate(int dx, int dy){
     x += dx;
     y += dy;
   }
   
/**
  * Duplicates its values and returns a new reference
  */
  public BHPoint duplicate(){
     return new BHPoint(x,y,label, val);      
  }
  
  
 /**
   * Compares this point to the provided one
   * First by x, then by y 
   *  -1 means comes before
   *  0 means equal
   *  1 comes after
   */
   public int compareTo(Object obj) throws ClassCastException {
      BHPoint p = (BHPoint) obj;
      
      int dx = x - p.x;
      int dy = y - p.y;
      
      if (dx != 0) {
         return dx;
      } else if (dy != 0) { 
         return dy;
      } 
      
      return 0;
   } //compareTo
   
  }// Point class