text = " \n*** BlobHandler Demo Macro  *** \n \n" + 
   "This macro demonstrates use of the BlobHandler where \n " + 
   " blobs whose perimeters are closer than 12 pixels apart \n " +
   " are combined.  The 'Labeled' and 'Contoured' images show the \n " + 
   " grouping by color. Note the 'Count' value in the results table. \n " +
   " 'Blob' in the table refers to the identifier of the first \n " + 
   " blob in the groupings \n"
print(text);
run("Blobs (25K)");
run("Duplicate...", "title=blobs-mask.gif");
setThreshold(125, 255);
run("Convert to Mask");
run("BlobHandler ", "area=0.0-Infinity circularity=0.0-1.0 combine=12 redirect=blobs.gif results clear label contour basic location shape volume ellipse grayscale moments");
