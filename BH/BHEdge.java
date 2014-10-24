//see http://algowiki.net/wiki/index.php/Edge

public class BHEdge implements Comparable<BHEdge> {
   
   final BHNode from, to;
   final double weight = 0.0;
   
   public BHEdge(final Node argFrom, final Node argTo, final double argWeight){
       from = argFrom;
       to = argTo;
       weight = argWeight;
   }
   
   public int compareTo(final BHEdge argEdge){
       return weight - argEdge.weight;
   }
}
