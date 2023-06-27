# RSAworkshop
 practice code and data for the workshop "Representational Similarity Analyses of fMRI data" 2023 at University Hamburg for members of the Research Training Group "Emotional Learning and Memory"

 ## ./Niftitools 
 includes some functions that are used for working with nifti-images (Copyright (c) 2014, Jimmy Shen)

 ## ./ROIbased
 includes scripts for computing and extracting model-based RSAs from regions of interest (ROIs)

 ## ./ROIs
 example regions of interest, mostly based on the Harvard-Oxford-Atlas and resized to the example beta images

 ## ./SLbased
 includes scripts for computing searchlight-based RSAs and exporting the results in nifti-images

 ## ./data
 example data that can be used for trying out the scripts

 ## ./modelbased
 includes scripts for comparing model-RSAs with conceptual models

 ## ./preparingRSAs
 example scripts for setting up single-trial and regular, condition-wise first-level analyses using SPM batch-scripts
 script for transforming single-trial trial beta-images into t-images 
 example script for transforming ROIs from MNI space into native space (in case an RSA in native space is planned)

 ## analysis.R
 example R-script for analyzing ROI-based RSA-results using linear mixed models and ANOVAs
