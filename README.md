# MosaicAnalysis
This packages contains several functions useful for the analysis of 2D othomosaics and orthoprojections:

1 - Buffer_Analysis (Patch Buffer Analysis) 
This function allows you to measure what is within the buffer of a given size of a distribution of contiguous patch. WARNING: You need the Python library GDAL to use this function. Go to: \code{https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/} for more informations on how to proceed.

2 - Pw_An_VPC (Power Analysis for Point Counts)
This function allows you to proceed with a power analysis of Point Counts on a annotated image.

3 - Pw_An_size_dis (Power Analysis for Size Distribution)
This function allows you to proceed with a power analysis for the 5 most commons metric used to describe a size distribution: mean, geometrical mean, kurtosis, skewness, coefficient of variation.

4 - IndPatchSize (Patch size distribution extractor)
This function allows you to extract the size and type of all the patches present in an image

5 - FromPictoRdata (Transform a PNG into a matrix in R)
This function allows you to transform an annotated PNG file into a matrixwith every cell value corresponding to a a specific type of organism
