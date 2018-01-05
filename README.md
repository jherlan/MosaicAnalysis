# MosaicAnalysis
This packages contains several functions useful for the analysis of 2D orthomosaics and orthoprojections:

## 1 - FromPictoRdata (Transform a PNG into a matrix in R)
This function allows you to transform an annotated PNG file into a matrix with every cell value corresponding to a specific type of organism.

## 2 - IndPatchSize (Patch size distribution extractor)
This function allows you to extract the size and type of all the patches present in an image.

## 3 - IndPatchSizePoly (Patch size distribution extractor)
This function allows you to extract the size and type of all the patches present in an image using the Python library GDAL. WARNING: You need the Python library GDAL to use this function. Go to: https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/ for more information on how to proceed.

## 4 - Pw_An_VPC (Power Analysis for Point Counts)
This function allows you to proceed with a power analysis of Point Counts using an annotated image.

## 5 - Pw_An_size_dis (Power Analysis for Size Distribution)
Using an annotated image, this function allows you to proceed with a power analysis for the 5 most commons metric used to describe a size distribution: mean, geometrical mean, kurtosis, skewness, coefficient of variation.

## 6 - Buffer_Analysis (Patch Buffer Analysis)
This function allows you to measure the relative composition of the buffer of a distribution of patch. WARNING: You need the Python library GDAL to use this function. Go to: https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/ for more information on how to proceed.

## 7 - PatchDistance (Distance between patch)
This function allows you to extract the distance between patches as well as their sizes, excluding colonies which are too close from the edges (i.e. wether directly touching an edge or for which an edge is closer than the identified closest neighbor). WARNING: You need the Python library GDAL to use this function. Go to: https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/ for more information on how to proceed.
