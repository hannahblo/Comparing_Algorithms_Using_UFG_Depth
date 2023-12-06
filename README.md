# Comparing Machine Learning Algorithms by Union-Free Generic Depth

This repository contains R-code and data sets corresponding to the "Comparing Machine Learning Algorithms by Union-Free Generic Depth" article written by H.Blocher, G. Schollmeyer, M. Nalenz and C. Jansen. We apply the ufg depth to two application examples considering the performance of classifiers.

The structure of the repository is as follows:
- Folder OpenML: contains the file "main_file.R" that exectutes the application on the data sets provided by the OpenMl repository (see: [https://www.openml.org/](https://www.openml.org/), accessed: 05.12.2023)
- Folder UCI_data: contains the file "main_file.R" that exectutes the application on the data sets provided by the UCI repository (Informations about the data see Jansen et al: [https://jmlr.org/papers/v24/22-0902.html](https://jmlr.org/papers/v24/22-0902.html), accessed: 05.12.2023). This folder contains also the R-file validate_classification.R that gives some information on the performance computation.

# Set up
- First, clone the repository please. 
- Second, please install all necessary R-packages listed in the beginning of the main_file.R files. Please note the following two points
  - For the computation of the linear programs, we used the R interface of gurobi optimizer, see [here](https://www.gurobi.com/) (accessed: 08.02.2023). This is a commercial
solver that offers a free academic licenses which can be found [here](https://www.gurobi.com/features/academic-named-user-license/) (accessed: 08.02.2023). To install this package, please follow the instructions there. A documentation can be found [here](https://www.gurobi.com/wp-content/plugins/hd_documentations/documentation/9.0/refman.pdf) (page 643ff) (accessed: 08.02.2023).
  - Some of the packages need to be installed using Bioconductor.
  - The oofos and ddandrda packages are still under development and need to be downloaded from their respective Github pages.
