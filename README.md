## BlobHandler

BlobHandler is a set of particle analysis tools for working with foreground objects (blobs) in 2d images.   

* A *Blob* is a set of contiguous foreground pixels.
* Two or more *Blob* objects may be associated to form a *Compound Blob*

BlobHandler is quite similar to ImageJ's built-in ParticleAnalyzer.  The basic 
work flow is extraction and measurement of foreground features (blobs) from bi-level images.
The BlobHandler_ class is the user interface for selecting various outputs and making
simple filtering of blobs based upon area and circularity.  This user interface class
is macro-friendly (see included BlobHandler_Demo_Macro.ijm) and easily modified.  

The workhorse class that handles blob detection and extraction is BlobHandler which can be extended.  This class
can be used by plugin programmers without the use of the user interface class (BlobHandler_).
See the BlobHandler_Demo (included) for examples on using BlobHandler programmatically.

While it is not necessary to use the BlobHandler_ class to run BlobHandler programmatically, the
 user interface class, like the demo, demonstrates how to use BlobHandler programatically. Under the hood it
 shows how to select features, filter and select outputs.  The user interface includes...

   + Filtering by area, min-max or just min (default = 0-Infinity)
   
   + Filtering by circularity, min-max or just min (default = 0-1)
   
   + Filter after combining blobs by area and circularity (default = false, filter before combining)
   
   + Combining objects by nearest proximity (default = 0, no combinations)
   
   + Selection of the Volume estimate algorithm (if Volume is selected as a feature below, default = auto)
   
   + Selection of the grayscale image to redirect measurements to (default = none)
   
   + Selection of the outputs (Results table cleared is the default)
   
   + Selection of the features to include in output if Results is checked

![BlobHandler GUI]()


### Installation

Clone the repository to an ImageJ plugins directory. 

### Usage

The repository includes compiled code (.jar and .class) to ease the user experience.  You should be able to simply restart ImageJ (or use "Help > Refresh Menus").  The navigate to "Plugins > BlobHandler > BlobHandler".




### References ###

##### REF1: For SVW (Sieracki, Viles and Webb) biovolume description

Biovolume calculation based upon Sieracki,Viles and Webb, 1989, "Algorithm to estimate cell biovolume using image analyzed microscopy", Cytometry 10:551-557.

##### REF2: for SVW biovolume source code

pmbiovol.c used for Sieracki, Viles and Webb (1989).  Contact Mike Sieracki at Bigelow Laboratory for Ocean Science

##### REF3: general guide and template

see ImageJ's Roi, ParticleAnalyzer, ResultsTable, Analyzer

##### REF4: *Excellent* for contour and labeling, Moments, Hu Moments, etc.

Wilhelm BURGER, Mark J. BURGE, Digital Image Processing An Algorithmic Introduction using Java, ISBN: 978-1-84628-379-6, Springer 2008

##### REF5: For contouring and labeling
F. Chang, C.-J. Chen, and C.-J. Lu, A linear-time component-labeling algorithm using contour tracing technique, Computer Vision and Image Understanding, vol. 93, no. 2, pp. 206-220, 2004. 

##### REF6: general reference, moments
B.K.P.Horn, Robot Vision ,  ISBN-10:0-262-08159-8, The MIT Press, 1986.

##### REF7: For AMI (affine moment invariants)
Flusser J., Suk T.: Pattern Recognition by Affine Moment Invariants. Pattern Recognition , 26 (1993), 1, 167-174. 
Flusser, J., Suk, T.: Affine Moment Invariants: A new tools for character recognition. Patter Recognition Letters, 15 (1994), 433-436.
http://zoi.utia.cas.cz/moment_invariants  

This book is available as of Jan 2010 - it adds more moments to the AMI collection but they are not included at this time.
Jan Flusser, Tomáš Suk, and Barbara Zitová, Moments and Moment Invariants in Pattern Recognition, Wiley & Sons Ltd., 2009 (312 pp., ISBN 978-0-470-69987-4)

##### REF8:  Elliptical Fourier Descriptors
Mark Nixon and Alberto Aguado, Feature Extraction and Image Processing, 2nd Ed, Academic Press 2008, ISBN 978-0-1237-2538-7 

##### REF9: general reference and the SKEDM (SKeletonization * EDM)  method of estimating the width of strand-like particles
John Russ, The Image Processing Handbook, 5th ed. CRC Press, 2007, ISBN 0-8493-7254-2

##### REF10: Paul Bourke's excellent geometry references...  
http://paulbourke.net/geometry/

##### REF11: Connected components (aka graphs and subgraphs)
http://www.kirupa.com/developer/actionscript/depth_breadth_search.htm
see also REF4 Chapter 11 for similar searching used in FloodFilling

##### REF12: ImageJ API, source code and examples
Rasband, W.S., ImageJ, U. S. National Institutes of Health, Bethesda, Maryland, USA, http://imagej.nih.gov/ij/, 1997-2015. 

### Who do I talk to? ####
* Ben Tupper btupper@bigelow.org
* Nicole Poulton npoulton@bigelow.org
* Michael Sieracki msieracki@bigelow.org