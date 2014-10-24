public class BHMoments {

   
   /** raw moments */
   public double m00, m10, m01, m11, m20, m02, m21, m12, m30, m03;
   /** centralized moments */
   public double u00, u10, u01, u11, u20, u02, u21, u12, u30, u03;
   /** normalized central moments */
   public double n00, n10, n01, n11, n20, n02, n21, n12, n30, n03;
   public int n;
   public double xc, yc;
   public double xsum, ysum, x2sum, y2sum, xysum ;
   public double x2ysum, xy2sum, x3sum, y3sum; 


   public BHMoments(){}
   public BHMoments(BHPoint[] p){
      compute(p)
   }

   public void zero(){
      n = 0;
      m00=0.0; m10=0.0; m01=0.0; m11=0.0; m20=0.0; m02=0.0; m21=0.0; m12=0.0; m30=0.0; m03=0.0;
      u00=0.0; u10=0.0; u01=0.0; u11=0.0; u20=0.0; u02=0.0; u21=0.0; u12=0.0; u30=0.0; u03=0.0; 
      n00=0.0; n10=0.0; n01=0.0; n11=0.0; n20=0.0; n02=0.0; n21=0.0; n12=0.0; n30=0.0; n03=0.0;
      xc=0.0; yc=0.0;
      xsum=0.0; ysum=0.0; x2sum=0.0; y2sum=0.0; xysum=0.0;
      x2ysum=0.0; xy2sum=0.0; x3sum=0.0; y3sum=0.0;
   }
   
   public void compute(BHPoint[] p){
      zero();
      computeSums(p);
      computeMoments();
   }
   public void computeSums(BHPoint[] p){
   
      n = p.length();
      for (int i = 0; i < n; i++){
      
         
      
      }
   
   
   }
}