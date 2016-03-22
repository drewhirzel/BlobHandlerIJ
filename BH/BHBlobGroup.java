import ij.*;
import ij.gui.*;
import java.util.Vector;
import ij.process.*;
import java.awt.*;
import java.util.*;
import graph.*;

/*
 * An Vector based container of BHCompoundBlob(s) which, in turn, are be comprised of
 * one or more BHSingleBlob(s).
 * <p> 
 * We use a Vector here because sometimes items are added to a location beyond
 * the "end" and so the intervening items must be 
 * created on the fly as null.  Blobs are uniquely numbered
 * with a label, which is used here as the indexed location
 * within the vector.
 * <p>
 * This class also has useful methods for filtering constituent objects based upon 
 * a number of features such as area, circularity,...  It would be simple to add others.
 * <p>
 * This class also manages the "combine" methods.
 * <p>
 * Contact: Mike Sieracki, msieracki@bigelow.org <br>
 * Author:  Ben Tupper, btupper@bigelow.org <br>
 * History: Developed over the course of 2009 and 2010 using ImageJ 1.43 and Mac OSX 10.5.8 <br>
 *   at Bigelow Laboratory for Ocean Science, www.bigelow.org
 */

public class BHBlobGroup {

/** For not very good reasons this container holds a copy of the measurment
  * flags used by the Blobs it contains 
  */   
  public int measureFlags = BHBlob.BASIC;
  
/**
  * The container
  */
  private Vector <BHCompoundBlob> vData = new Vector<BHCompoundBlob>();
  
/**
  * Generic constructor - label = 0, type = 1
  */
  public BHBlobGroup(){
  }
 
/**
  * Constructor with specified blob and location
  */
  public BHBlobGroup(int index, BHCompoundBlob b){
      add(index, b);
  }
/**
  * Constructor with first item to append
  */  
  public BHBlobGroup(BHCompoundBlob b){
    add(b.label, b);
  }
  
/**
  * Clears the data
  */
  public void clear(){
    vData.clear();
  }
  
/**
  * Get as array
  */
  public BHCompoundBlob[] toArray(){
    //String[] sv = ( String[] )v.toArray( new String[ v.size() ] );
    return (BHCompoundBlob[]) vData.toArray(new BHCompoundBlob[vData.size()]);
  }
/**
  * Gets a Blob by label
  */
  public BHCompoundBlob get(int label){
    //IJ.log("Get Blob " + label + " is there? " + hasBlob(label));
    BHCompoundBlob blob = null;
    for (int i = 0; i < vData.size(); i++){
      blob = vData.get(i);
      if ((blob != null) && (blob.label == label)){ break;}
    }
    return blob;
  } 

/**
  * Returns the unique labels
  */
  public int[] getLabels(){
   int count = 0;
   BHCompoundBlob blob;
   int[] x = new int[vData.size()];
   for (int i = 0; i < vData.size(); i++){
       blob = vData.get(i);
       if (blob != null){
          x[count] = blob.label;
          count++;
       }
   }
   int[] labels = new int[count];
   System.arraycopy(x, 0, labels, 0, count);
   return labels;
  }
  
/**
  * Returns a index LUT for labels.  So that the index location maybe 
  * quickly determined for a set of labels
  */
  public int[] getLabelsLUT(){
   
   int[] labels = getLabels();
   int[] lut = new int[vData.size()];
   for (int i = 0; i< lut.length; i++){lut[i] = -1;}
   for (int i = 0; i < labels.length; i++){lut[labels[i]] = i;}
    return lut;
  }//getLabelsLUT
/**
  * Count returned
  */
  public int countBlobs(){
    return count();
  }//countBlobs
 
/**
  * Returns the count of valid blobs
  */
  public int count(){
    int count = 0;
    for (int i = 0; i < vData.size(); i++){
      if (vData.get(i)!= null) { count++;}
    }
    return count;
  } 
/**
  * Adds a BHCompoundBlob to the specifed location
  */
  public void add(int index, BHCompoundBlob b){
    if (vData.size() < index) {
      vData.setSize(index);
    }
    vData.set(index, b);
  }//add
  
/**
  * Adds a BHCompoundBlob to the end
  */
  public void add(BHCompoundBlob b){
    vData.add(b);
  }
/**
  * Adds a SingleBlob to the labeled Blob at the specifed index
  */  
  public void add(int index, BHSingleBlob sb){
    
    if (vData.size() <= index) {
      vData.setSize(index+1);
      BHCompoundBlob blob = new BHCompoundBlob(sb);
      vData.set(index, blob);
      return;
    }
    
    if (hasBlob(index)){
      vData.get(index).add(sb);
      return; 
    }
    
    vData.set(index, new BHCompoundBlob(sb));
 } //add
 
/**
  * Places the SingleBlob at it's label index
  */
  public void add(BHSingleBlob sb){
    add(sb.label, sb);
  }
  

/**
  * Appends an array of BHCompoundBlobs
  */
  public void add(BHCompoundBlob[] blobs){
    for (int i = 0; i < blobs.length; i++){
      add(blobs[i].label, blobs[i]);
    }
  }// add an array of contours
 
 
/**
  * Removes an element specified by label - returns the element that
  * was removed and replaces its location with null.  Note that this
  * is different than the "remove" methods for Vector which do not 
  * perform the replacement.
  * @param label The label of the blob to remove
  * @return the BHCompoundBlob element that was removed (and replaced by null);
  */
  public BHCompoundBlob removeBlob(int label){
    int ix = getLabelIndex(label);
    BHCompoundBlob b = null;
    if (ix >= 0){
        b = vData.get(ix);
        vData.setElementAt(null, ix);
    }
    return b;
  }//removeBlob
  

/**
  * Moves a BHCompoundBlob from idxA to idxB - places a null in the original location idxA
  */
  public void moveBlob(int idxA, int idxB){
    BHCompoundBlob blobA = removeBlob(idxA);
    vData.setElementAt(blobA, idxB);  
  }  //moveBlob
  
/**
  * indicates if the specified blob exists
  */ 
  public boolean hasBlob(int label){
   if (label >= vData.size()){ return false;}
   int idx = getLabelIndex(label);
   if (idx < 0){return false;}
   BHCompoundBlob b =  vData.get(label);
   return (b != null);
  }

/**
  * Gets the zero-based index for the specifed label.
  * @param theLabel the label to find the location of.
  * @return a zero-based index of the label's location, if not found then -1 is returned
  */
  public int getLabelIndex(int theLabel){
    int idx = -1;
    BHCompoundBlob b = null;
    for (int i = 0; i < vData.size(); i++){
      b = vData.get(i);
      if ((b !=null) && (b.label == theLabel)){
        idx = i;
        break;
      }
    }
    return idx;
  }
  
/**
  * Merge (via combine(one, two)) the list of blobs into the first listed
  * in the argument.  For example...  merge({3,6,11})
  * will add blobs 6 and 11 to 3.
  * @param labelList the list of label indices to combine
  */
  public void combine(int[] labelList){
    if (labelList.length > 1) {
      for (int i = 1; i < labelList.length; i++){   
        combine(labelList[0], labelList[i]);
      }
    }
  }
/**
  * Combines the second to the first and voids the position of the second
  * @param one the label of the first blob
  * @param two the label of the second blob
  */
  public void combine(int one, int two){
   
    if (!hasBlob(one)) { return ;} 
    if (!hasBlob(two)) { return ;}
    
    BHCompoundBlob b1 = get(one);
    BHCompoundBlob b2 = get(two);
    
    b1.add(b2);
    
    vData.setElementAt(null, two);
    b1.computeStats();
    
  }
  
/**
  * Adds a point to the blob at "label" - if the blob doesn't exists it is created
  */
  public int addPoint(int label, BHPoint point){
     if (label >= vData.size()){
        BHSingleBlob b = new BHSingleBlob(label);
        add(b);
        return b.addPoint(point);
     } else {
       BHCompoundBlob b = (BHCompoundBlob) vData.get(label);
       if (b == null){
         b = new BHCompoundBlob(label);
       }
      return b.addPoint(point); 
     }   
  }//addPoint
  
/**
  * Adds a point - if the blob doesn't exists it is created
  */
  public int addPoint(int label, int x, int y){
     return addPoint(label, new BHPoint(x,y, label));
  }
  
/**
  * Adds a point - if the blob doesn't exists it is created
  */
  public int addPoint(int label, int x, int y, double val){
     return addPoint(label, new BHPoint(x,y, label, val));
  }

/**
  * Adds a point - if the blob doesn't exists it is created
  */
  public int addPoint(int label, int x, int y, float val){
     return addPoint(label, new BHPoint(x,y, label, (double) val));
  }
    
/**
  * Adds a contour - the contour already contains its own inside/outside
  * flag so the blob will store it properly.
  * If the blob at [label] is does not exists it will be created.
  */ 
  public void addContour(BHPolygon c){
    //IJ.log("BlobGroup adding contour to " + c.label);
    BHCompoundBlob b = vData.get(c.label);
    if (b == null) {return ;}
    b.addContour(c); 
  }
/**
  * Moves all members by dx, dy
  */
  public void translate(int dx, int dy){
  
      if (count() == 0) { return;}
      BHCompoundBlob b;
      int[] labels = getLabels();
      //IJ.showStatus("Translating blobs");
      for (int i = 0; i < labels.length; i++){
        //IJ.showProgress(i, labels.length-1);
        b = get(labels[i]);
        b.translate(dx,dy);
      }
      //IJ.showStatus(" ");
      //IJ.showProgress(2.0);
    //
    //if (vData.size() != 0) {
    //  BHCompoundBlob[] b = toArray();
    //  for (int i = 0; i<vData.size(); i++){
    //    if (b[i] != null) {b[i].translate(dx,dy);}
    //  }
    //}
  }//translate

/**
  * Sets the side flag for each compound object based upon the BHBlob.NONE, BHBlob.LEFT,...
  * @param w the width of the original image in pixles
  * @param h the height of the original image in pixels
  */
  public void setSideFlag(int w, int h){
    setSideFlag(0, 0, w, h);
  }   
/**
  * Sets the side flag for each compound object based upon the BHBlob.NONE, BHBlob.LEFT,...
  * @param xs the origin of the original image in pixels
  * @param ys the origin of the original image in pixels
  * @param w the width of the original image in pixles
  * @param h the height of the original image in pixels
  */
  public void setSideFlag(int xs, int ys, int w, int h){
    
    if (count() == 0) { return;}
    BHCompoundBlob b;
    int[] labels = getLabels();
    //IJ.showStatus("Computing side touches");
    for (int i = 0; i < labels.length; i++){
      //IJ.showProgress(i, labels.length-1);
      b = get(labels[i]);
      b.side = b.isTouching(xs, ys, w, h);
    }
    //IJ.showStatus(" ");
    //IJ.showProgress(2.0);
    
  }//setSideFlag  
  
/**
  * Computes the minimum distances separating the outside
  * contour(s) of the each blob - returns a double[][] that is only 1/2 populated
  * (since distance from A to B is the same as B to A). The returned is just
  * like a distance table on a road map.  The diagonal is always 0.
  * NOTE - the returned array is [nblob][nblob] not [size][size]thus where blobs
  * don't exists (null).  Use getLabels() to retrieve the labels associated
  * with each row/column
  * 
  */  
  public double[][] computeMinimumDistanceTable(){
    int nblobs = count();
    if (nblobs == 0) {return null;}
    if (nblobs == 1) {return new double[1][1];}
    double[][] dist = new double[nblobs][nblobs];
    int[] labels = getLabels(); //this will be nblobs long
    BHCompoundBlob leftHandBlob;
    //IJ.showStatus("Computing nearest neighbor distances for "+nblobs+" blobs");
    for (int y = 0; y < nblobs; y++){
      //IJ.showProgress(y, nblobs-1);
      leftHandBlob = get(labels[y]);
      for (int x = (y+1); x < nblobs; x++){
        dist[y][x] = leftHandBlob.computeMinimumDistanceTo(get(labels[x]));
        dist[x][y] = dist[y][x];//mirror image
      } //x-loop
    }// y-loop
    //IJ.showStatus(" ");
    //IJ.showProgress(2.0);
    return dist;
  }//computeMinimumDistanceTable


 public int[][] computeMergeTable(double mergeIfCloserThan){
    return computeMergeTable6(mergeIfCloserThan);
 }
 
 public int[][] computeMergeTable6(double mergeIfCloserThan){
    
    if (count() == 0){
      return null;
    }
    //get the square matrix of Euclidean distances
    double[][] dist = computeMinimumDistanceTable();
    int n = dist.length;
    // make an index table
    int[][] adj = new int[n][n];
    int[] labels = getLabels();//this is dist.length long
    int[] lut = getLabelsLUT();    
    
     
    //first make a look-up table based upon labels
    //rather than distance
    for (int y = 0; y < adj.length; y++){
      for (int x = 0; x < adj[0].length; x++){ // from 0
        if ((dist[y][x] != 0) && (dist[y][x] < mergeIfCloserThan)){
          adj[y][x] = labels[x];
          adj[x][y] = labels[y];//symmetric
        } 
      }// x-loop
    }//y-loop

   //create an instance of the Graph maker
	Graph	g = new graph.Graph();
	//extract the sub graphs
	int nSubs = g.extractSubgraphs(adj);
	//
	return g.getAdjacencyList();    
 }
 
// public int[][] computeMergeTable5(double mergeIfCloserThan){
//    
//    if (count() == 0){
//      return null;
//    }
//    //get the square matrix of Euclidean distances
//    double[][] dist = computeMinimumDistanceTable();
//    int n = dist.length;
//    int[] nPerRow = new int[dist.length];
//    // make an index table
//    int[][] idx = new int[n][n];
//    int[] labels = getLabels();//this is dist.length long
//    int[] lut = getLabelsLUT();
//    int[] start = new int[n];
//    
//    
//     
//    //first make a look-up table based upon labels
//    //rather than distance
//    for (int y = 0; y < idx.length; y++){
//      start[y] = labels[y];
//      //ptr[y] = labels[y];  //the pointer
//      //idx[y][0] = labels[y];  //self-label the first entry
//      for (int x = 0; x < idx[0].length; x++){ // from 0
//      //for (int x = y+1; x < idx[0].length; x++){ // from y+1
//        if (dist[y][x] == 0){
//         //idx[y][x] = labels[y];
//        } else if ((dist[y][x] != 0) && (dist[y][x] < mergeIfCloserThan)){
//          idx[y][x] = labels[x];
//          idx[x][y] = labels[y];//labels[x];//symmetric
//          if (idx[y][x] < start[y]) {start[y] = idx[y][x];} 
//        } 
//      }// x-loop
//      //if (start[y] == n){start[y]= labels[y];} //means you are alone
//    }//y-loop
//    
//    
//    
//    IJ.log("mergeIfCloserThan = " + mergeIfCloserThan);
//    IJ.log("N =" + n);
//    String fmt = "%3d";
//   
//    String hdr = String.format("%9s",":");
//    for (int i = 0; i < n; i++){ hdr = hdr + String.format(fmt,labels[i]) + " ";}
//
//    //make the following true to spit out the raw connectivity table
//    if (true) {
//       IJ.log("idx table in bh.computeMergeTable() BEFORE");
//       IJ.log(hdr);
//       for (int y = 0; y < idx.length; y++){
//         String s = String.format(fmt, labels[y]) +" " + String.format(fmt, start[y]) + ": ";
//         for (int x = 0; x < idx[0].length; x++){
//         s = s + String.format(fmt, idx[y][x]) + " ";
//         }
//         IJ.log(s);
//      }
//    }
// 
//   //this is the list of possible nodes
//   BHNode[] nodeArr = new BHNode[n];
//   for (int i = 0; i < n; i++){ nodeArr[i] = new BHNode(labels[i]);}
//   //this is the container for the subgraphs
//   Vector<Vector<BHNode>> graphs = new Vector<Vector<BHNode>>();
//   
//   Vector<BHNode> subG = new Vector<BHNode>();
//   for (int y = 0; y < n; y++){
//      for (int x = 0; x < n; x ++){  
//         if ((idx[y][x] != 0) && (subG.size() == 0)){
//            graphs.add(subG);
//            BHNode node = nodeArr[idx[y][x]];
//            node.visited = true;
//            subG.add(node);
//         } 
//         if ((idx[y][x] != 0) && (subG.size() != 0)){
//            DFS(idx, x, y, nodeArr, subG);
//         }
//         
//         
//      }
//   }
//
//   int[][] ret = new int[graphs.size()][];
//   for (int i = 0; i < ret.length; i++){
//      subG = (Vector<BHNode>) graphs.get(i);
//      BHNode[] r = (BHNode[]) subG.toArray(new BHNode[0]); 
//      ret[i] = new int[r.length];
//      for (int j = 0; j < r.length; j++){ ret[i][j] = r[i].label;}
//   }
//   return ret;
//  }
//
// 
//  public void DFS(int[][] idx, int x0, int y0, BHNode[] nodeArr, Vector<BHNode> subG){
//  
//      Stack<BHNode> s = new Stack<BHNode>();  
//      int label = idx[y0][x0];
//      BHNode theNode = nodeArr[label-1];
//      //theNode.visited = true;
//      s.push(theNode);
//      while (!s.isEmpty()){
//         BHNode n = s.pop();
//         if (n.visited == false){
//            n.visited = true;
//            subG.add(n);
//            for (int x=0; x<idx[n.label-1].length; x++){
//               //if (nodeArr[x].visted == false){
//                  s.push(nodeArr[x]);
//               //}//if not visited
//            }//step along the row
//         } //if not visited
//         
//      
//      }//while
//  
//  
//  
//  
//  }//DFS
// 
// 
// 
// public int[][] computeMergeTable4(double mergeIfCloserThan){
//    
//    if (count() == 0){
//      return null;
//    }
//    //get the square matrix of Euclidean distances
//    double[][] dist = computeMinimumDistanceTable();
//    int n = dist.length;
//    int[] nPerRow = new int[dist.length];
//    // make an index table
//    int[][] idx = new int[n][n];
//    int[] labels = getLabels();//this is dist.length long
//    int[] lut = getLabelsLUT();
//    int[] start = new int[n];
//    
//    //first make a look-up table based upon labels
//    //rather than distance
//    for (int y = 0; y < idx.length; y++){
//      start[y] = labels[y];
//      //ptr[y] = labels[y];  //the pointer
//      //idx[y][0] = labels[y];  //self-label the first entry
//      for (int x = 0; x < idx[0].length; x++){ // from 0
//      //for (int x = y+1; x < idx[0].length; x++){ // from y+1
//        if (dist[y][x] == 0){
//         idx[y][x] = labels[y];
//        } else if ((dist[y][x] < mergeIfCloserThan)){
//          idx[y][x] = labels[x];
//          idx[x][y] = labels[y];//labels[x];//symmetric
//          if (idx[y][x] < start[y]) {start[y] = idx[y][x];} 
//        } 
//      }// x-loop
//      //if (start[y] == n){start[y]= labels[y];} //means you are alone
//    }//y-loop
//    
//    IJ.log("mergeIfCloserThan = " + mergeIfCloserThan);
//    IJ.log("N =" + n);
//    String fmt = "%3d";
//   
//    String hdr = String.format("%9s",":");
//    for (int i = 0; i < n; i++){ hdr = hdr + String.format(fmt,labels[i]) + " ";}
//
//    //make the follwoing true to spit out the raw connectivity table
//    if (true) {
//       IJ.log("idx table in bh.computeMergeTable() BEFORE");
//       IJ.log(hdr);
//       for (int y = 0; y < idx.length; y++){
//         String s = String.format(fmt, labels[y]) +" " + String.format(fmt, start[y]) + ": ";
//         for (int x = 0; x < idx[0].length; x++){
//         s = s + String.format(fmt, idx[y][x]) + " ";
//         }
//         IJ.log(s);
//      }
//    }
//  
//  
//  
//   
//  
//  
//   // now we surf start and idx from the bottom upward
//   for (int y = idx.length-1; y > 0; y--){
//   
//      if (start[y] != labels[y]) {
//         //if start is not equal to the index of the row
//         //then there is some lower valued label that we can transfer these
//         // items to
//         //for (int x = 0; x < idx[y].length; x++){
//         //   if (idx[y][x] != 0){
//               //if (idx[y][x] != start[y]){
//                  //if it is non-zero AND not the start value then move it up vertically
//                  //to where start points to
//                  int ptr = y;
//                  int prevPtr = y;
//                  int theLabel = labels[y];
//                  int theStart = start[y];
//                  while (theStart != theLabel){
//                     theStart = start[ptr];
//                     theLabel = labels[ptr];
//                     ptr = start[ptr]-1;
//                     if (theStart != theLabel){
//                        copyRow(idx[prevPtr], idx[ptr], false);
//                        prevPtr = ptr;
//                     }
//                  }
//                  
//                  //copyRow(idx[y], idx[ptr], false);
//                  for (int x = 0; x < idx[y].length;x++){
//                     idx[y][x] = 0;
//                  }     
//                  //idx[ptr][x] = idx[y][x];
//                  
//               //} else {
//                  //otherwise, if non-zero AND it is the start value
//                  //insert the label in the [start][label] location which
//                  //is somewhere above this row
//                //  idx[start[y]-1][labels[y]-1] = labels[y];
//               //}   
//            //} //if idx[y][x] != 0
//            
//            //if the start value is less than the label value for this row
//               //then zero out the element 
////               if (start[y] < labels[y]) {
////                  //idx[y][x] = 0-idx[y][x];
////                  idx[y][x] = 0;
////               }
//         //} //x-loop
//         //start[y] = 0-start[y];
//      } else {
//         //hmmmmmm   what does it mean to be here if we do nothing   
//   
//      }
//   
//   } //y-loop
//  
//    //make the follwoing true to spit out the raw connectivity table
//    if (true) {
//       IJ.log("idx table in bh.computeMergeTable() AFTER");
//       IJ.log(hdr);
//       for (int y = 0; y < idx.length; y++){
//         String s = String.format(fmt, labels[y]) +" " + String.format(fmt, start[y]) + ": ";
//         for (int x = 0; x < idx[0].length; x++){
//         s = s + String.format(fmt, idx[y][x]) + " ";
//         }
//         IJ.log(s);
//      }
//    }
//     
//   int[] connCount = new int[n];
//   int cblobCount = 0;
//   for (int y = 0; y < idx.length; y++){
//      for (int x = 0; x < idx[y].length; x++){
//         if (idx[y][x] != 0) {connCount[y]++;}
//      }
//      if (connCount[y] != 0){ cblobCount++;}
//   }
//   
//   int[][] ret = new int[cblobCount][];
//   
//   int j = 0;
//   for (int y = 0; y < idx.length; y++){
//      if (connCount[y] > 0){
//         ret[j] = new int[connCount[y]];
//         int k = 0;
//         for (int x = 0; x < idx[y].length; x++){
//            if (idx[y][x] != 0){
//               ret[j][k] = idx[y][x];
//               k++;
//            }
//         }
//         j++;
//      }
//   }
//   
//   return ret;
//  }
//
//
// public int[][] computeMergeTable3(double mergeIfCloserThan){
//    
//    if (count() == 0){
//      return null;
//    }
//    //get the square matrix of Euclidean distances
//    double[][] dist = computeMinimumDistanceTable();
//    int n = dist.length;
//    int[] nPerRow = new int[dist.length];
//    // make an index table
//    int[][] idx = new int[n][n];
//    int[] labels = getLabels();//this is dist.length long
//    int[] lut = getLabelsLUT();
//    int[] start = new int[n];
//    
//    //first make a look-up table based upon labels
//    //rather than distance
//    for (int y = 0; y < idx.length; y++){
//      start[y] = labels[y];
//      //ptr[y] = labels[y];  //the pointer
//      //idx[y][0] = labels[y];  //self-label the first entry
//      for (int x = 0; x < idx[0].length; x++){ // from 0
//      //for (int x = y+1; x < idx[0].length; x++){ // from y+1
//        if (dist[y][x] == 0){
//         idx[y][x] = labels[y];
//        } else if ((dist[y][x] < mergeIfCloserThan)){
//          idx[y][x] = labels[x];
//          idx[x][y] = labels[y];//labels[x];//symmetric
//          if (idx[y][x] < start[y]) {start[y] = idx[y][x];} 
//        } 
//      }// x-loop
//      //if (start[y] == n){start[y]= labels[y];} //means you are alone
//    }//y-loop
//    
//    IJ.log("mergeIfCloserThan = " + mergeIfCloserThan);
//    IJ.log("N =" + n);
//    String fmt = "%3d";
//   
//    String hdr = String.format("%9s",":");
//    for (int i = 0; i < n; i++){ hdr = hdr + String.format(fmt,labels[i]) + " ";}
//
//    //make the follwoing true to spit out the raw connectivity table
//    if (true) {
//       IJ.log("idx table in bh.computeMergeTable() BEFORE");
//       IJ.log(hdr);
//       for (int y = 0; y < idx.length; y++){
//         String s = String.format(fmt, labels[y]) +" " + String.format(fmt, start[y]) + ": ";
//         for (int x = 0; x < idx[0].length; x++){
//         s = s + String.format(fmt, idx[y][x]) + " ";
//         }
//         IJ.log(s);
//      }
//    }
//  
//  
//   // now we surf start and idx from the bottom upward
//   for (int y = idx.length-1; y > 0; y--){
//   
//      if (start[y] != labels[y]) {
//         //if start is not equal to the index of the row
//         //then there is some lower valued label that we can transfer these
//         // items to
//         for (int x = 0; x < idx[y].length; x++){
//            if (idx[y][x] != 0){
//               //if (idx[y][x] != start[y]){
//                  //if it is non-zero AND not the start value then move it up vertically
//                  //to where start points to
//                  int ptr = y;
//                  int theLabel = labels[y];
//                  int theStart = start[y];
//                  while (theStart != theLabel){
//                     theStart = start[ptr];
//                     theLabel = labels[ptr];
//                     ptr = start[ptr]-1;
//                  }
//                       
//                  idx[ptr][x] = idx[y][x];
//                  
//               //} else {
//                  //otherwise, if non-zero AND it is the start value
//                  //insert the label in the [start][label] location which
//                  //is somewhere above this row
//                //  idx[start[y]-1][labels[y]-1] = labels[y];
//               //}   
//            } //if idx[y][x] != 0
//            
//            //if the start value is less than the label value for this row
//               //then zero out the element 
//               if (start[y] < labels[y]) {
//                  //idx[y][x] = 0-idx[y][x];
//                  idx[y][x] = 0;
//               }
//         } //x-loop
//         //start[y] = 0-start[y];
//      } else {
//         //hmmmmmm   what does it mean to be here if we do nothing   
//   
//      }
//   
//   } //y-loop
//  
//    //make the follwoing true to spit out the raw connectivity table
//    if (true) {
//       IJ.log("idx table in bh.computeMergeTable() AFTER");
//       IJ.log(hdr);
//       for (int y = 0; y < idx.length; y++){
//         String s = String.format(fmt, labels[y]) +" " + String.format(fmt, start[y]) + ": ";
//         for (int x = 0; x < idx[0].length; x++){
//         s = s + String.format(fmt, idx[y][x]) + " ";
//         }
//         IJ.log(s);
//      }
//    }
//     
//   int[] connCount = new int[n];
//   int cblobCount = 0;
//   for (int y = 0; y < idx.length; y++){
//      for (int x = 0; x < idx[y].length; x++){
//         if (idx[y][x] != 0) {connCount[y]++;}
//      }
//      if (connCount[y] != 0){ cblobCount++;}
//   }
//   
//   int[][] ret = new int[cblobCount][];
//   
//   int j = 0;
//   for (int y = 0; y < idx.length; y++){
//      if (connCount[y] > 0){
//         ret[j] = new int[connCount[y]];
//         int k = 0;
//         for (int x = 0; x < idx[y].length; x++){
//            if (idx[y][x] != 0){
//               ret[j][k] = idx[y][x];
//               k++;
//            }
//         }
//         j++;
//      }
//   }
//   
//   return ret;
//  }
//
//
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
//  public int[][] computeMergeTable2(double mergeIfCloserThan){
//    
//    if (count() == 0){
//      return null;
//    }
//    //get the square matrix of Euclidean distances
//    double[][] dist = computeMinimumDistanceTable();
//    int n = dist.length;
//    int[] nPerRow = new int[dist.length];
//    // make an index table
//    int[][] idx = new int[n][n];
//    int[] labels = getLabels();//this is dist.length long
//    int[] lut = getLabelsLUT();
// 
////    for (int i=0; i < labels.length; i++){
////      IJ.log("Labels[" + i + "]=" + labels[i]);
////    }
//    Vector<BHNodeVector> allv = new Vector<BHNodeVector>();
//    
//    //first make a look-up table based upon labels
//    //rather than distance
//    for (int y = 0; y < idx.length; y++){
//      //ptr[y] = labels[y];  //the pointer
//      //idx[y][0] = labels[y];  //self-label the first entry
//      //for (int x = 0; x < idx[0].length; x++){ // from 0
//      for (int x = y+1; x < idx[0].length; x++){ // from y+1
//        if ((dist[y][x] > 0.0) && (dist[y][x] < mergeIfCloserThan)){
//          idx[y][x] = labels[x];
//          idx[x][y] = labels[y];//labels[x];//symmetric
//        }
//      }
//    }
//    
//    IJ.log("mergeIfCloserThan = " + mergeIfCloserThan);
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
//  
//    //start a new nodevector and populate it with the first node
//    BHNodeVector nv = new BHNodeVector();
//   
//    for (int y = 0; y < idx.length; y++){
//      //starting a new row.  If the value in the first position is > 0
//      //then we haven't been here yet  - start a new nodeVector
//      if (idx[y][0] > 0){
//         nv = new BHNodeVector();
//         nv.add(new BHNode(idx[y][0]));
//         idx[y][0] = -1; //reassign to -1
//         allv.add(nv);
//      }   
//      for (int x = 0; x < idx[0].length; x++){
//         if (idx[y][x] > 0){
//            chaseConnectivity(idx, idx[y][x]-1, nv);
//         }
//      }
//    }//y-loop
//    
//    
//   
//   if (true){
//      IJ.log("idx table after initial pass");
//      for (int i=0; i < allv.size(); i++){
//         ((BHNodeVector) allv.get(i)).log();
//      }
//   }
//   
//   int[][] ret = new int[allv.size()][];
//   int[] count = new int[allv.size()];
//   for (int i = 0; i < allv.size(); i++){
//      ret[i] = ((BHNodeVector)allv.get(i)).getLabels();
//   }
//   
//   
//   
//   
//   return ret;
//  }
//
//
//  public void chaseConnectivity(int[][] idx, int y, BHNodeVector nv){
//   for (int x = 0; x < idx[y].length; x++){
//      if (idx[y][x] > 0){
//         //if the node label is not zero and it hasn't been added before
//         //then chase it
//         if (nv.add(new BHNode(idx[y][x]))){
//            chaseConnectivity(idx, idx[y][x]-1, nv);
//            //idx[y][x] = -1;
//         }
//      } //if not zero
//   }//i loop along row
//  
//  }//chaseConnectivty
//  
//  public int[][] computeMergeTable1(double mergeIfCloserThan){
//    
//    if (count() == 0){
//      return null;
//    }
//    //get the square matrix of Euclidean distances
//    double[][] dist = computeMinimumDistanceTable();
//    int n = dist.length;
//    int[] nPerRow = new int[dist.length];
//    // make an index table
//    int[][] idx = new int[n][n];
//    int[] labels = getLabels();//this is dist.length long
//    int[] lut = getLabelsLUT();
//    // this is a pointer of label associations
//    //int[] ptr = new int[idx.length]; 
//       
//    //first make a look-up table based upon labels
//    //rather than distance
//    for (int y = 0; y < idx.length; y++){
//      //ptr[y] = labels[y];  //the pointer
//      idx[y][0] = labels[y];  //self-label the first entry
//      //for (int x = 0; x < idx[0].length; x++){ // from 0
//      for (int x = y+1; x < idx[0].length; x++){ // from y+1
//        if ((dist[y][x] > 0.0) && (dist[y][x] < mergeIfCloserThan)){
//          idx[y][x] = labels[x];
//          //idx[x][y] = lut[labels[y]];//labels[x];//symmetric
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
//    // idx holds the adjacency labels
//    // so along row 1 the objects directly adjacent are listed
//    // to reduce replication the matrix is only populated on the upper right
//    // now we scan through each row-portion (from diagonal toward right)
//    for (int y = 0; y < idx.length; y++){
//      for (int x = (y+1); x < idx[y].length; x++){
//         //if the occupied entry is not empty
//        if (idx[y][x] > 0){ 
//          //get a copy of the lut value - this is a label
//          int lutidx = lut[idx[y][x]];
//          if (idx[lutidx][0] == idx[y][x]){
//            //we are here if it is copy of itself
//            idx[lutidx][0] = idx[y][0];
//            //ptr[lutidx] = idx[y][0];
//          } else {
//            //we are here 
//            idx[y][0] = idx[lutidx][0];
//            //ptr[y] = idx[lutidx][0];
//          }
//        } //if
//      } //x-loop
//     nPerRow[lut[idx[y][0]]]++; //increment the bin count
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
//      if (nPerRow[y] > 0) { connCount++;}
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
//      if (nPerRow[y] > 0){
//        ret2[nconn] = new int[nPerRow[y]];
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

   
   private void copyRow(int[] src, int[] dst){
      copyRow(src, dst, true);
   }
   private void copyRow(int[] src, int[] dst, boolean zero){
      for (int i = 0; i < src.length; i++){
         dst[i] = src[i];
         if (zero) {src[i] = 0;}
      }
   }


/**
  * 
  */
/**
  * Labels the image processor to make it a contour mask
  */
  public void labelContoursProcessor(ImageProcessor ip, int inOrOut){
    int[] labels = getLabels();
    BHCompoundBlob blob;
    if (labels.length == 0) { return;}
    for (int i = 0; i < labels.length; i++){
      blob = get(labels[i]);
      blob.labelContours(ip, blob.label, inOrOut);
    }    
  }//labelPixels  
/**
  * Labels the image processor to make it a count mask
  */
  public void labelPixelsProcessor(ImageProcessor ip){
    int[] labels = getLabels();
    if (labels.length == 0) { return;}
    for (int i = 0; i < labels.length; i++){
      get(labels[i]).labelPixels(ip);
    }    
  }//labelPixels
  
/**
  * Labels the image processor to make it a count mask
  */
  public void labelPixelsProcessor(ImageProcessor ip, int col){
    int[] labels = getLabels();
    if (labels.length == 0) { return;}
    for (int i = 0; i < labels.length; i++){
      get(labels[i]).labelPixels(ip, col);
    }    
  }//labelPixels
   
/**
  * Removes items that are smaller / bigger than the area specified
  * @param minArea the inclusive minimum area or 0.0, Double.NEGATIVE_INFINITY or 
  *   Double.NaN if no lower limit.
  * @param maxArea the inclusive maximum area or Double.NaN 
  *   Double.MAX_VALUE, or Double.POSITIVE_INFINITY if no upper limit
  */
  public void filterByArea(double minArea, double maxArea){
    BHCompoundBlob b;
    double[] r = checkRange(minArea, maxArea);
    boolean removeMe = false;
    int[] labels = getLabels();
    for (int i = 0; i < labels.length; i++){
      removeMe = false;
      b = vData.get(labels[i]);
      if (b != null){
        if ((b.area < r[0]) || (b.area > r[1])) { removeMe = true;}
      }
      if (removeMe == true) { removeBlob(b.label);}
    }//i-loop
    
  }//filterByArea

/**
  * Removes items that are smaller / bigger than the Circ specified
  * @param minCirc the inclusive minimum Circ or 0.0, Double.NEGATIVE_INFINITY or 
  *   Double.NaN if no lower limit.
  * @param maxCirc the inclusive maximum Circ or Double.NaN 
  *   Double.MAX_VALUE, or Double.POSITIVE_INFINITY if no upper limit
  */
  public void filterByCirc(double minCirc, double maxCirc){
    BHCompoundBlob b = null;
    double[] r = checkRange(minCirc, maxCirc);
    boolean removeMe = false;
    int[] labels = getLabels();
    for (int i = 0; i < labels.length; i++){
      b = vData.get(labels[i]);
      removeMe = false;
      if (b != null){
        if (b.circ < r[0]){ removeMe = true;}
        if (b.circ > r[1]){ removeMe = true;}         
      }
      if (removeMe == true) { removeBlob(b.label);}
    }//i-loop
    
  }//filterByCirc  

/**
  * Filter by area and circularity
  * @param minArea the inclusive minimum area or 0.0, Double.NEGATIVE_INFINITY or 
  *   Double.NaN if no lower limit.
  * @param maxArea the inclusive maximum area or Double.NaN 
  *   Double.NaN if no lower limit.
  * @param minCirc the inclusive minimum Circ or 0.0, Double.NEGATIVE_INFINITY or 
  *   Double.NaN if no lower limit.
  * @param maxCirc the inclusive maximum Circ or Double.NaN 
  *   Double.MAX_VALUE, or Double.POSITIVE_INFINITY if no upper limit
  */
  public void filterByAreaAndCirc(double minArea, double maxArea,double minCirc, double maxCirc){
    filterByArea(minArea, maxArea);
    filterByCirc(minCirc, maxCirc);
  }//filterByAreaAndCirc
  

/**
  * Filter items with MaxESD outside of the specified range
  * @param min The inclusive shortest MaxESD 
  * @param max the inclusive longest MaxESD
  */
  public void filterByMaxESD(double min, double max){
    double[] r = checkRange(min,max);
    BHCompoundBlob b = null;
    boolean removeMe = false;
    int[] labels = getLabels();
    for (int i = 0; i < labels.length; i++){
      b = vData.get(labels[i]);
      removeMe = false;
      if (b != null){
        if (b.feret[b.FERETMAX] < r[0]){ removeMe = true;}
        if (b.feret[b.FERETMAX] > r[1]){ removeMe = true;}         
      }
      if (removeMe == true) { removeBlob(b.label);}
    }//i-loop 
  }  
 
/**
  * Filter items with MaxESD outside of the specified range
  * @param min The inclusive shortest MaxESD 
  * @param max the inclusive longest MaxESD
  */
  public void filterByMinESD(double min, double max){
    double[] r = checkRange(min,max);
    BHCompoundBlob b = null;
    boolean removeMe = false;
    int[] labels = getLabels();
    for (int i = 0; i < labels.length; i++){
      b = vData.get(labels[i]);
      removeMe = false;
      if (b != null){
        if (b.feret[b.FERETMIN] < r[0]){ removeMe = true;}
        if (b.feret[b.FERETMIN] > r[1]){ removeMe = true;}         
      }
      if (removeMe == true) { removeBlob(b.label);}
    }//i-loop 
  }  

 
/**
  * Filter items with MeanESD outside of the specified range
  * @param min The inclusive shortest MeanESD 
  * @param max the inclusive longest MeanESD
  */
  public void filterByMeanESD(double min, double max){
    double[] r = checkRange(min,max);
    BHCompoundBlob b = null;
    boolean removeMe = false;
    int[] labels = getLabels();
    for (int i = 0; i < labels.length; i++){
      b = vData.get(labels[i]);
      removeMe = false;
      if (b != null){
        if (b.feret[b.FERETMEAN] < r[0]){ removeMe = true;}
        if (b.feret[b.FERETMEAN] > r[1]){ removeMe = true;}         
      }
      if (removeMe == true) { removeBlob(b.label);}
    }//i-loop 
  }  
         
/**
  *
  * Removes items that are smaller / bigger than the ABD specified
  * @param minABD the inclusive minimum area or 0.0, Double.NEGATIVE_INFINITY or 
  *   Double.NaN if no lower limit.
  * @param maxABD the inclusive maximum area or Double.NaN 
  *   Double.MAX_VALUE, or Double.POSITIVE_INFINITY if no upper limit
  */
  public void filterByABD(double minABD, double maxABD){
    double[] r = checkRange(minABD, maxABD);
    double minArea = Math.PI * (r[0]/2.0)*(r[0]/2.0);
    double maxArea = Math.PI * (r[1]/2.0)*(r[1]/2.0);
    filterByArea(minArea, maxArea);
  }//filterByABD
  
/**
  * Checks the input values a, b and returns A and B if a or B are NaN or infinite
  * @param a the lower range value 
  * @param b the upper range value
  * @param A the replacement for a if it is NaN or infinite
  * @param B the replacement for b if it is NaN or infinite
  * @return a two element range value
  */
  public double[] checkRange(double a, double b, double A, double B){
    double[] r = {a, b};
    if (Double.isNaN(a) || Double.isInfinite(a)){r[0] = A;}
    if (Double.isNaN(b) || Double.isInfinite(b)){r[1] = B;}
    
    return r;
  }//checkRange
  
/**
  * Checks the input values a, b and returns 0.0 or Double.MAX_VALUE for either 
  * if they are NaN or infinite
  */    
  public double[] checkRange(double a, double b){
    return checkRange(a, b, 0.0, Double.MAX_VALUE);
  }
    
/**
  * Remove if touching the specified side, caller can specify up to 4
  * sides to test - or less than four using BHBlob.NONE as the pass through value.
  * @param w the width of the image
  * @param h the height of the image
  * @param edge1 a bit-style flag of the first side not to touch.
  * @param edge2 a bit-style flag of the second side not to touch.
  * @param edge3 a bit-style flag of the third side not to touch.
  * @param edge4 a bit-style flag of the fourth side not to touch.
  * <p>
  * KEY: 
  * 0 = BHBlob.NONE, 1 = LEFT, 2 = TOP, 4 = RIGHT, 8 = BOTTOM or some combination of these.
  * <p>
  * For example...
  * group.removeIfTouchingEdge(640, 480, BHBlob.RIGHT, BHBlob.BOTTOM, BHBlob.NONE, BHBlob.NONE);
  */
  public void removeIfTouchingEdge(int w, int h, int edge1, int edge2, int edge3, int edge4){
    if (countBlobs() <= 1) { return;}
    int[] labels = getLabels();
    BHCompoundBlob b;
    int isTouching;
    for (int i = 0; i < labels.length; i++){
      b = get(labels[i]);
      isTouching = b.isTouching(w,h);
      if ( ((isTouching & edge1) != 0) || 
        ((isTouching & edge2) != 0) || 
        ((isTouching & edge3) != 0) ||
        ((isTouching & edge4) != 0)) {
          //IJ.log("  Remove label=" + labels[i] + " flag=" + isTouching);
          removeBlob(labels[i]);
      }
    }
  }
    
/**
  * Eliminates all particles but the one with the largest area
  * In the event of a tie, the first is retained
  */  
  public void removeAllButLargest(){
    if (countBlobs() <= 1) { return;}
    int[] labels = getLabels();
    double area = -Double.MAX_VALUE;
    int ix = 0;
    BHCompoundBlob b = null;
    //identify the one with the biggest area
    for (int i = 0; i < labels.length; i++){
       b = get(labels[i]);
       if (b.area > area) {
          area = b.area;
          ix = i;
       }
    }//i-loop
    //now remove all but the one at ix
    for (int i = 0; i<labels.length; i++){
      if (i != ix){ removeBlob(labels[i]);}
    }
   }//removeAllButLargest

/**
  * Calls the computeVolume method for specified BHCompoundBlob
  * @param theLabel the label identfier
  */
  public void computeVolume(int theLabel){
   BHCompoundBlob b = get(theLabel);
   if (b != null){ b.computeVolume();} 
  } //computeVolume
/**
  * Calls the computeVolume method for each constituent BHCompoundBlob
  */
  public void computeVolume(){
    int[] labels = getLabels();
    for (int i = 0; i < labels.length; i++){
      computeVolume(labels[i]);
    }
  }//computeVolume


/**
  * Sets the volume computation method for the specified blob
  * @param theLabel the label ID of the blob
  * @param theMethod See BHBlob for possible values like VOL_SVW, VOL_AUTO, ...
  */
  public void setVolumeMethod(int theLabel, int theMethod){
   BHCompoundBlob b = get(theLabel);
   if (b != null){ b.setVolumeMethod(theMethod);} 
  } //computeVolume

/**
  * Sets the volume computation method for each compound blob
  * @param method See BHBlob for possible values like VOL_SVW, VOL_AUTO, ...
  */
  public void setVolumeMethod(int method){
    int[] labels = getLabels();
    for (int i = 0; i < labels.length; i++){
      setVolumeMethod(labels[i], method);
    }    
  }
  
  public Overlay showAllContours(int which, ImagePlus imp){
    Overlay overlay = new Overlay();
    return showAllContours(which,  overlay);
 }
  
  public Overlay showAllContours(int which, Overlay overlay){
    BHCompoundBlob b = null;
    BHPolygon[] poly;
    int[] labels = getLabels();
    for (int i = 0; i < labels.length; i++){
         b = get(labels[i]);
         poly = b.getContours(which);
         if ((poly != null) && (poly.length > 0) ){
            for (int j = 0; j < poly.length; j++) {overlay.add(poly[j].getPolygon());}
         }
    }     
    return overlay;
  }
  
  public Overlay showAllEllipses(Overlay overlay){
    BHCompoundBlob b = null;
    PolygonRoi poly;
    int[] labels = getLabels();
    for (int i = 0; i < labels.length; i++){
         b = get(labels[i]);
         poly = b.getEllipseRoi(b);
         overlay.add(poly);
    }
    return overlay;
  }
  
}//BHBlobGroup class