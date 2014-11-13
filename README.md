# README #

BlobHandler is a set of particle analysis tools for working with foreground objects (blobs) in 2d images.   

* A *Blob* is a set of contiguous foreground pixels.
* Two or more *Blob* objects may be associated to form a *Compound Blob*


### References ###

REF1: For SVW (Sieracki, Viles and Webb) biovolume descriotion

Biovolume calculation based upon Sieracki,Viles and Webb, 1989, "Algorithm to estimate cell biovolume using image analyzed microscopy", Cytometry 10:551-557.

REF2: for SVW biovolume source code

pmbiovol.c used for Sieracki, Viles and Webb (1989).  Contact Mike Sieracki at Bigelow Laboratory for Ocean Science

REF3: general guide and template

see ImageJ's Roi, ParticleAnalyzer, ResultsTable, Analyzer

REF4: *Excellent* for contour and labeling, Moments, Hu Moments, etc.

Wilhelm BURGER, Mark J. BURGE, Digital Image Processing An Algorithmic Introduction using Java, ISBN: 978-1-84628-379-6, Springer 2008

REF5: For contouring and labeling
F. Chang, C.-J. Chen, and C.-J. Lu, A linear-time component-labeling algorithm using contour tracing technique, Computer Vision and Image Understanding, vol. 93, no. 2, pp. 206-220, 2004. 

REF6: general reference, moments
B.K.P.Horn, Robot Vision ,  ISBN-10:0-262-08159-8, The MIT Press, 1986.


REF7: For AMI (affine moment invariants)
Flusser J., Suk T.: Pattern Recognition by Affine Moment Invariants. Pattern Recognition , 26 (1993), 1, 167-174. 
Flusser, J., Suk, T.: Affine Moment Invariants: A new tools for character recognition. Patter Recognition Letters, 15 (1994), 433-436.
http://zoi.utia.cas.cz/moment_invariants

This book is available as of Jan 2010 - it adds more moments to the AMI collection but they are not included at this time.
Jan Flusser, Tomáš Suk, and Barbara Zitová, Moments and Moment Invariants in Pattern Recognition, Wiley & Sons Ltd., 2009 (312 pp., ISBN 978-0-470-69987-4)

REF8:  Elliptical Fourier Descriptors
Mark Nixon and Alberto Aguado, Feature Extraction and Image Processing, 2nd Ed, Academic Press 2008, ISBN 978-0-1237-2538-7 

REF9: general reference and the SKEDM (SKeletonization * EDM)  method of estimating the width of strand-like particles
John Russ, The Image Processing Handbook, 5th ed. CRC Press, 2007, ISBN 0-8493-7254-2

REF10: Paul Bourke's excellent geometry references...  
http://local.wasp.uwa.edu.au/~pbourke/geometry/

REF11: Connected components (aka graphs and subgraphs)
http://www.kirupa.com/developer/actionscript/depth_breadth_search.htm
see also REF4 Chapter 11 for similar searching used in FloodFilling


### Who do I talk to? ####
* Ben Tupper btupper@beigelow.org
* Nicole Poulton npoulton@bigelow.org
* Michael Sieracki msieracki@bigelow.org