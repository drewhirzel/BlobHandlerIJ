A BlobHandler.jar is provided in the repos, but if you wish to make your own so that other 
plugins can have easy access to the BlobHandler class and its derivatives.

(1) Remove all the classes in the plugins/BlobHandler/BH directory

rm *.class

(2) Compile from ImageJ the BlobHandler.java located in plugins/BlobHandler/BH directory

(3) Navigate to the BH directory and Build a jar for BlobHandler and then remove the classes

cd <ImageJ-path>/plugins/BlobHandler/BH
jar cvfM ../BlobHandler.jar *.class graph/*.class
rm *.class graph/*.class

(4) Compile BlobHandler_.java and BlobHandler_Demo.java separately.


Now you'll have in the BlobHandler directory...

   About_BlobHandler.class
   BlobHandler_.class
   BlobHandler_Demo.class
   BlobHandler.jar

plus the 
   java source code for About_BlobHandler, BlobHandler_, BlobHandler_Demo
   macro example (BlobHandler_Demo_Macro.ijm)
   the api folder
   the example image (test.tif)
   the references documentation (references.rtf)
   the index.html (introduction page and one image blobhandler-gui.png)
   a text file of miscellaneous notes (notes.txt)
   this document (make-jar-doc.txt)



(Outdated) Here are some notes on creating a jar for the BlobHandler and BH classes followed 
by some brief notes on making the api documentation.

//-------
Use the following to compile the classes into a jar
//-------

$ jar cvfM BlobHandler_.jar BHBlob.class BHSingleBlob.class BHCompoundBlob.class BHBlobGroup.class BHPoint.class BHPointGroup.class BHPolygon.class BlobHandler.class BlobHandler_.class
$ rm *.class



//----------
Make the documentation
//----------

 javadoc -d /Applications/ImageJ/plugins/BlobHandler/api /Applications/ImageJ/plugins/BlobHandler/BH/*.java
 