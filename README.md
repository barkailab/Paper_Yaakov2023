# Paper_Yaakov2023
Thank you for your interest in our research!
This repo contains the scripts and matfiles to recapitulate our analysis.
The name of each scripts indicates the corresponding Figure in the paper.
Each script is divided into subsections corresponding to individual subfigures.

The mat files contain the variables necessary to run the script and are loaded in the corresponding sections in the script.
However, in the scripts the "load command" often refers to absolute paths. it is necessary to adjust these when running the code - but should be straightforward.
Also the naming/meaning of the variables should be intuitive but please feel free to contact us with specific questions.

Moreover, due to size constrains the actual data is not part of this repo. instead of the data files are empty variables with the right name.
These variables need to be ppopulated with the processed data from the corresponding GEO dataset (GSE231810).
The affected matfiles are:

chapRelocGithub.mat
  normProfile contains the procesed data in the samples of sTable. (the first dimension, rows, correpsond to all bases in the yeast genome, the second dimension columns are the different samples, according to sTable).

dataFs-singleGithub.mat
  data contians the processed data of the samples in smeta. (see above)
  normProfile containts the mean data of the  Chaperone Chec seq profiles indicated in chapNames
  
dataFScombGithub.mat
  data (same as data in dataFS-singleGithub.mat)

TFdataGithub.mat
  normProfile (same as normProfile in chapRelocGithub.mat)

Please reach out to us if you have any further questions.
Enjoy analysing the data.
