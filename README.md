# DuraClone-gating
 Automated gating of Basic and Bcell panel of DuraClone. This repository has the automated gating pipeline for 2 of the DuraClone panels.


Files can downloaded through the flowRepository link decribed on the manuscript. Make sure both LMD/FCS files and attachments will be downloaded in a same folder. Everything should be in  EXACTLY 1 folder.

OrganizeFiles.R needs to be run after downloading the files so that files are organized and in correct formatting for the automated analysis. Make sure you change the the path of where the files are, and where you want the files to be moved after they are organizd.

Make sure the metadata is downloaded from flowRepository, as it is needed for the step above.

The automated pipeline is based on flowDensity automated gating, and the manual gating strategy. For a detailed information check the manuscript. 

The latest version of flowDensity can be downloaded from devel branch of BioConductor, here:
https://master.bioconductor.org/packages/devel/bioc/html/flowDensity.html.

The automated pipeline is based on different datasets and also cytometer; in order to run the pipeline, only Run_automated.R needs to be modified and run.

Before running the code, make sure path are updated to correspond to a right folder in your directory ( i.e., input.path, res.path, and path of sourced codes). 
dataset, type and panel need to be set.


Look for "###This needs to be modified" within the R codes to update the paths, so that they correspod to correct paths in your local directory.
