<html>
<head>
<title>BlobHandler</title>
</head>
<body bgcolor = "#ffffff">

<font color="#224488" face="Helvetica, Arial" size=+1>
<b>BlobHandler for Particle Analysis</b>
</font>

<p>

<blockquote>

<table cellspacing=5>

<tr valign=top><td ><b>Author:</b></td><td >
Ben Tupper (btupper at bigelow.org) and Mike Sieracki (msieracki at bigelow.org) <br>
Bigelow Laboratory for Ocean Science - www.bigelow.org
</td></tr>

<tr valign=top><td ><b>History:</b></td><td >
2010-02-18: First version<br>
</td></tr>

<tr valign=top><td ><b>Source:</b></td><td >

<A HREF="download/BlobHandler.zip">BlobHandler.zip</A>
</td></tr>

<tr valign=top><td ><b>Installation:</b></td><td >
Download <A HREF="download/BlobHandler.zip">BlobHandler.zip</A> to the plugins folder and unzip. Restart ImageJ or use the menu selection Help > Update Menus.
</td></tr>

<tr valign=top><td ><b>References:</b></td><td >
See 'references.rtf' for a list of resources we have relied upon.  Three were indispensible: 
   (1) the ImageJ API, examples and source code and (2) Digital Image Processing An Algorithmic Introduction using Java, Burger and Berge, 
   and (3) the many gems and thoughtful contributions found on ImageJ's mailing list archives.
</td></tr>

<tr valign=top><td ><b>Description:</b></td><td >
   BlobHandler is quite similar to ImageJ's built-in ParticleAnalyzer.  The basic 
   work flow is extraction and measurement of foreground features (blobs) from bi-level images.
   The BlobHandler_ class is the user interface for selecting various outputs and making
   simple filtering of blobs based upon area and circularity.  This user interface class
   is macro-friendly (see included BlobHandler_Demo_Macro.ijm) and easily modified.  
   <p>
   The workhorse class that handles blob detection and extraction is BlobHandler which can be extended.  This class
   can be used by plugin programmers without the use of the user interface class (BlobHandler_).
   See the BlobHandler_Demo (included) for examples on using BlobHandler programmatically.
   <p>
   While it is not necessary to use the BlobHandler_ class to run BlobHandler programmatically, the
   user interface class, like the demo, demonstrates how to use BlobHandler programatically. Under the hood it
   shows how to select features, filter and select outputs.  The user interface includes...<br>
   <ul>
   <li>Filtering by area, min-max or just min (default = 0-Infinity)
   <li>Filtering by circularity, min-max or just min (default = 0-1)
   <li>Filter after combining blobs by area and circularity (default = false, filter before combining)
   <li>Combining objects by nearest proximity (default = 0, no combinations)
   <li>Selection of the Volume estimate algorithm (if Volume is selected as a feature below, default = auto)
   <li>Selection of the grayscale image to redirect measurements to (default = none)
   <li>Selection of the outputs (Results table cleared is the default)
   <li>Selection of the features to include in output if Results is checked (see notes below).
   </ul><p><p>
   <center><img src="blobhandler-gui.png"></center>
   <p>
   <p>
   <center><b><u>Feature selection</u></b> <br> (Many of these have better desciptions in the online ImageJ manual)<br><br></center>
   <b>BASIC</b> features
     <ul>
        <li> Blob - the unique identifier of the blob(s), useful when combining blobs
        <li> Count - the number of individual blob(s) that make up this compound blob
        <li> Area - the area of the blob(s)
        <li> Perim - the total perimeter length of the blob(s)
     </ul>
   <b>LOCATION</b> information relative to image
     <ul>
        <li> BX, BY, BW and BH  The bounding box description as x0, y0, width and height
        <li> XC and YC  The centroid of the blob(s)
        <li> XStart, YStart The coordinates of the first pixel encountered in a scan.
        <li> Side A flag indicating which side(s) of the image the blob(s) touch(es), where <br>
             0 = none, 1 = left, 2 = top, 4 = right 8 = bottom or some additive combination of these
     </ul>
   <b>SHAPE</b> features
     <ul>
        <li> FeretMax The longest dimension found rotating a box at 2-degree increments around the blob(s) (ala <i>rotating caliper</i>).
        <li> FeretMin the shortest dimension ...
        <li> FeretMed the median long feret measurement for the 180 rotations
        <li> FeretX, FeretY and Theta use these to reconstruct a line parallel to FeretMax and
           passing through a point on the perimeter of the rotated blob(s).
        <li> Orient  The orientation of the object derived from mechanical moments.
        <li> Circ - circularity - (4*PI*Area)/(Perim^2)
        <li> Round - roundness - (4*Area)/(PI*FeretMax^2)
        <li> Solid - solidity - (AreaWithinOuterPerimeter)/(ConvexHullArea)
        <li> Compact - compactness - SQRT(4*Area/PI)/(FeretMax)
        <li> Extent - (Area)/(BoundingBoxArea)
        <li> Aspect - (FeretMax)/(FeretMin)
        <li> VoidFraction - (AreaWithinInnerPerimeter)/(AreaWithinOuterPerimeter)
        <li> NEuler - Euler Shape Number - (NumberOuterPerimeters)-(NumberInnerPerimeters)
     </ul>
   <b>VOLUME</b>
     Volume estimate is useful for those working with cells in the aquatic sciences.
     <ul>
     <li><b>SVW</b> We often use a method to estimate volume from Sieracki, Viles and Webb (1989, see references). This
     method (SVW) has been shown to work well with relatively convex shapes.
     <li><b>ABD</b>To this method we add a simple ABD method which converts Area to an equivalent 
     radius (if the pixels were arranged as a circle) then computes volume as 4/3*PI*r^3
     <li><b>SKEDM</b> And we add a technique that combines the skeletonization with the Euclidian distance map (SKEDM). This
     is briefly described in Russ (2007, see refs) but not necessarily in the context of computing volume.
     We have found this method works well for objects with solidity lower than 0.4.
     <li> <b>AUTO</b> Finally, the appropriate method can be found using the AUTO method which uses SKEDM for solidity < 0.4
     and SVW in all other cases.
     <li><b>NTricky</b> is reported with volume.  This is the count of single blobs in the compound blob that qualified for special 
        treatment with the SKEDM algorithm - even if that algorithm was not selected.
     </ul>
     
   <b>ELLIPSE</b> fitting parameters
     <ul>
        <li> Major - the length of the semi-major axis
        <li> Minor - the length of the semi-minor axis
        <li> Ecc - eccentricity of the ellipse
        <li> Use Orient (above) for the orientation of the long axis
     </ul>
   <b>GREYSCALE</b> information - most meaningful if redirected to a greyscale image
     <ul>
        <li> GMin Minimum grey value
        <li> GMax Maximum grey value
        <li> GMean Mean grey value
        <li> GVar Greyscale variance (0 for binary images)
        <li> GSkew Greyscale skewness (NaN for binary images)
        <li> GKurt Greyscale kurtosis (NaN for binary images)
     </ul>
   <b>MOMENT INVARIANTS</b> for binary blobs(s). We suggest that the references are consulted for details
     <ul>
      <li> Hu moments are described in Burger and Burge, 2008
      <li> Affine Moment Invariants (AMI) are described in Flusser et al, 2009
     </ul>   
   <p>
   <u><b>About Perimeters</b></u><br>
   Users should note that blob perimeters (outer and inner) are defined by the 
   coordinate of the cracks separating foreground pixels from backgound.  This is 
   the boundary found using ImageJ's Wand class. There are various methods for 
   computing the perimter based upon these coordinates, and  we have selected a
   simple-minded approach of summing the lengths between vertices.  See some discussion
   about an alternative method <a href = http://www.dentistry.bham.ac.uk/landinig/software/software.html>here</a>
   <p>
   <b><u>Underlying Conceptual Framework</b></u><br> 
   <ul>
   <li>A bi-level 2d image is comprised of foreground (blobs) and background features. 
   <li>Contiguous groups of foreground pixels are considered SingleBlobs.
   <li>SingleBlobs have one outer perimeter and zero or more inner perimeters (aka 'contours').
   <li>CompoundBlobs are comprised of one or more SingleBlobs.
   <li>CompoundBlobs have as many outer perimeters are the number of SingleBlobs it contains.
   <li>CompoundBlobs have zero or more inner perimeters.
   <li>Features of CompoundBlobs are reportable in the results table.
   </ul>

</td></tr>


</table>

</blockquote>
<p><b>|<a href="index.html">Plugins</a> | <a href="../index.html">Home</a> |</b>

</body>
</html>

