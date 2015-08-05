#Pvalue Weighting Package for Matlab

This MATLAB package contains open source implementations of several p-value weighting methods, including Spjotvoll, exponential and Bayes weights. These are methods for improving power in multiple testing via the use of prior informaton. 

Please see the readme.pdf file for more information, including examples of using the code. For the details, consult the paper *Optimal Multiple Testing Under a Gaussian Prior on the Effect Sizes* by Dobriban, Fortney, Kim, and Owen:  http://arxiv.org/abs/1504.02935

In addition, this package has all the code to reproduce the computational and data analysis results of the above paper.

The R package **pweight** (https://github.com/edgardobriban/pweight) is more actively maintained and contains more functionality.

## Pipeline for Pvalue Weighting of GWAS data

To reproduce the data analysis results of the paper *Optimal Multiple Testing Under a Gaussian Prior on the Effect Sizes* by Dobriban, Fortney, Kim, and Owen:  http://arxiv.org/abs/1504.02935,
we provide the data analysis pipeline used there.

##Dependencies and Licensing
* Requires MATLAB. It was tested on R2014a/b.
* License: GPL-3
* Parts of this package depend on MATLAB packages written by third-party users. These are available from the MATLAB Central archive. They are also included in the current package, in /Code/Helper Code folder. Each of them has a BSD License, which explicitly allows re-distribution in the present form, and they are copyright of their respective owners. The packages are: displaytable, fdr_bh, figtitle, kde2d, savetightfigure. 

### Who do I talk to? ###

* Repo owner: Edgar Dobriban. dobriban@stanford.edu
