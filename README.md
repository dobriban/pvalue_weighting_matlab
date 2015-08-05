#Pvalue Weighting Package for Matlab

This MATLAB package contains open source implementations of several p-value weighting methods, including Spjotvoll, exponential and Bayes weights. These are methods for improving power in multiple testing via the use of prior informaton. 

Please see the readme.pdf file for more information, including examples of using the code. For the details, consult the paper *Optimal Multiple Testing Under a Gaussian Prior on the Effect Sizes* by Dobriban, Fortney, Kim, Owen:  http://arxiv.org/abs/1504.02935

In addition, this package has all the code to reproduce the computational and data analysis results of the above paper.

This package is provided only for convenience of Matlab users. The R package **pweight** (https://github.com/edgardobriban/pweight) is more actively maintained and contains more functionality.

# Pipeline for Pvalue Weighting of GWAS data

To reproduce the data analysis results of the paper *Optimal Multiple Testing Under a Gaussian Prior on the Effect Sizes* by Dobriban, Fortney, Kim, Owen:  http://arxiv.org/abs/1504.02935,
We provide the data analysis pipeline.
This pipeline can be used to efficiently and reproducibly perform standard data analyses with p-value weighting of GWAS data.  It has the following features:

**High-level functions that perform a standard set of analyses:** The user specifies the prior and current data sets, and the methods described in our paper are automatically performed.

**Customizable parameters:** The user can can change many aspects of the analysis by specifying appropriate parameters. For instance, it is possible to change the parameters of the weighting schemes.

**Reproducibility:** The results of analyses are recorded in a standard format and directory structure, which, if coupled with a disciplined workflow, enable reproducible analysis.

**Extensibility:**  It is possible to add new data sets, weighting schemes, and perform custom workflows. Some of this will require writing new code.

### How do I get set up? ###

* Pull the repository to a folder on a local machine or a cluster
* To inspect the results of existing analyses, go to the folder Data Analysis/[AnalysisName]/Results and browse text and image files.
* To reproduce existing analyses, you first need to download and process the original data files. Due to data access policies these could not be included in the package.
  1. To download data files, go to Data/Raw/[Data Set Name] and consult the description.txt file for the description web link for the appropriate data sets.
  2. Download and unpack each data set into its own Data/Raw/[Data Set Name] folder
  3. Run the MATLAB scripts in the Data/Raw/[Data Set Name]/Code folder to process the raw data into MATLAB files. They will be deposited in Data/Raw/.
  4. Repeat the steps above for each data set.
* To reproduce or change existing analyses go to the folder Data Analysis/[AnalysisName] and run or modify analysis.m
* Add new analyses (for instance on newly specified pairs of GWAS data sets) by creating a new folder Analysis/[AnalysisName]
* Add new data sets by creating a new folder Data/Raw/[DataName]/ where you can copy the data for a first processing. Using the existing data sets as a template, process the new data into a MATLAB format, which will then be saved in the folder Data/Processed.
* Write new core analysis code by adding scripts to the Code/ folder.
* The change directory statement in analysis.m is configured to the authors' local machine and should be modified accordingly.
* Requires MATLAB. It was tested on R2014a.


### Who do I talk to? ###

* Repo owner: Edgar Dobriban. dobriban@stanford.edu
