# DuraClone-gating
 Automated gating of Basic and Bcell panel of DuraClone. This repository has the automated gating pipeline for 2 of the DuraClone panels.


Files can downloaded through the flowRepository Link decribed on the manuscript. Make sure both  LMD/FCS files and attachments will be downloaded in a same folder. Everything should be in 1 folder.
OrganizeFiles.R needs to be run after downloading the files so that files are organized and in correct formatting for the automated analysis. 
Make sure the metadata is downloaded from flowRepository, as it is needed for the step above.

The automated pipeline is based on flowDensity automated gating, and the manual gating strategy. For a detailed information check the manuscript.
The automated pipeline is based on different datasets and also cytometer.
To run the pipeline, only Run_automated.R needs to be modified and run.
Before running the code, make sure path are updated to correspond to a right folder in your directory. 
dataset, type and panel need to be set.
