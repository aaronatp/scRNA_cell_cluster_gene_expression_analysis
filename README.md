# Tools for analyzing scRNA-seq counts data
This repository contains tools that I have found helpful for analyzing scRNA-seq data generated with Seurat. In the top of the 'seurat_object_to_csv.Rmd' file, I have included some code for porting the data used here from Seurat objects to csv files. I hope that the tools here will be useful and helpful for other people analyzing scRNA-seq data. Please let me know if you think there is any more functionality that would be useful to have in this sort of package.

If you don't know what Anaconda is or if you don't have it installed on your computer, have a look at this (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html - this may be confusing. If you need further help, try to find a video on Youtube on installing Anaconda or shoot me a message). If you don't know what virtual environments are, have a look at this (https://stackoverflow.com/questions/41972261/what-is-a-virtualenv-and-why-should-i-use-one). Once you have Anaconda installed and you know what virtual environments are, come back and keep reading. The configuration of my virtual environment can be found in the 'environment.yaml' file. 

To load configure your virtual environment using these settings: first, go to the terminal app on your computer and type "git clone clone this GitHub repository.  run '''conda env create --name gene_expression_analysis -f /path/to/this/environment.yaml''' from the command line.
