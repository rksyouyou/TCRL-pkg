# tcrl: Association Analysis of T-Cell Receptor Repertoire and Clinical Phenotypes

[![License](https://img.shields.io/badge/license-LGPL--2.0-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.0.html)

## Overview

T cell receptors (TCRs) are essential for adaptive immunity, and advances in genome technology now enable analysis of the TCR repertoire at the individual sequence level. Studying TCR repertoire associations with clinical phenotypes can provide insights into immune-mediated diseases, yet methods for such analyses remain underdeveloped. We introduce **TCR-L**, a statistical tool for evaluating associations between the TCR repertoire and disease outcomes. Built on a mixed-effect model, **TCR-L** accounts for both **fixed effects**, representing explicitly extracted TCR sequence features, and **random effects**, capturing hidden sequence characteristics. Independent statistical tests assess these effects separately before combining their results.

Simulation studies confirm that **TCR-L** effectively controls type I error and outperforms methods considering only fixed or random effects ([https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04690-2](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04690-2)). Application to real-world data from a skin cutaneous melanoma study reveals a significant association between TCR repertoire and patient survival. By incorporating both observable and latent TCR features, **TCR-L** offers a robust approach for investigating TCR repertoire associations with disease outcomes.

The **TCR-L** has been implemented into the **tcrl** package. The **tcrl** package provides a suite of functions for analyzing associations between the T-cell receptor (TCR) repertoire and clinical phenotypes. It is designed to help researchers evaluate the diversity and composition of the TCR repertoire and its relationship to clinical outcomes, such as disease states or treatment responses. The package includes methods for both continuous and binary phenotypes, incorporating homology matrices to account for correlations between subjects' TCR repertoires.

### Features:
- Functions for calculating amino acid frequencies in TCR repertoires.
- Score-based tests for association between TCR repertoire features and clinical phenotypes.
- Support for both continuous and binary outcomes.
- Tools for calculating TCR repertoire homology between subjects using substitution matrices such as BLOSUM and PAM.

## Installation

You can install the tcrl package by downloading the tcrl_0.1.0.tar.gz file to your local system. Once downloaded, install it from your local file path using the following command in R:


```{r}
install.packages("path/to/tcrl_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Usage

### Example 1: Amino Acid Frequency Calculation

```{r}
# Load example dataset
data("example.data")
attach(example.data)
str(example.data)
 List of 4
  $ TCRdat:'data.frame': 118 obs. of  3 variables:
   ..$ Subject.ID: chr [1:118] "1" "1" "1" "1" ...
   ..$ AAseq     : chr [1:118] "CASSNSTKTVFF" "CASSWKRCLLFF" "CSARASAMHVIMFF" "CSARFLRVFF" ...
   ..$ Abundance : chr [1:118] "1" "1" "5" "5" ...
  $ Y     : num [1:20, 1] 0.301 1.226 -0.327 0.155 1.144 ...
   ..- attr(*, "dimnames")=List of 2
   .. ..$ : NULL
   .. ..$ : chr "Y_b"
  $ X     : num [1:20, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
  $ W     : num [1:20, 1] 1.8 -4.5 -3.5 -3.5 2.5 -3.5 -3.5 -0.4 -3.2 4.5 ...
   ..- attr(*, "dimnames")=List of 2
   .. ..$ : chr [1:20] "A" "R" "N" "D" ...
   .. ..$ : NULL
```


`example.data` is a simulated dataset designed to illustrate the analysis of T-cell receptor (TCR) repertoires in relation to clinical and biological outcomes. It includes TCR sequence information, a continuous response variable, and covariates for 20 patients, making it a useful resource for demonstrating function usage in association studies.

The dataset is structured as a list with multiple components. `TCRdat` is a data frame containing three key columns: `Subject.ID`, which identifies individual patients; `AAseq`, representing amino acid sequences in the TCR repertoire; and `Abundance`, a numeric vector indicating the frequency of each sequence per patient. The dataset also includes `Y`, a numeric vector representing continuous response variables, such as clinical or biological measurements associated with each patient. Additionally, `X` is a matrix of covariates where each row corresponds to a patient and each column represents different variables, such as age, gender, or other relevant clinical factors. Lastly, `W` is a numeric vector providing Kyte and Doolittle hydrophobicity scores for each amino acid sequence.

This dataset is intended to facilitate the evaluation of associations between TCR diversity and clinical outcomes. By incorporating TCR sequences, their abundance, and patient-specific covariates, `example.data` provides a structured framework for exploring immunological and bioinformatics analyses.

```{r}
# Compute amino acid frequency
freq_matrix <- AAfreq(TCRdat[,1], TCRdat[,2], as.numeric(TCRdat[,3]))

# View the result
print(freq_matrix)

detach(example.data)
```

### Example 2: TCR Repertoire Homology Calculation

```{r}
# Load example dataset
data("example.data")
attach(example.data)

# Compute the homology matrix using BLOSUM62
S <- seqhom(TCRdat[,1], TCRdat[,2], as.numeric(TCRdat[,3]), "BLOSUM62")

# View the homology matrix
print(S)

detach(example.data)
```

### Example 3: Score Test for Continuous Outcomes

```{r}
# Load example dataset
data("example.data")
attach(example.data)

# Extract TCR repertoire features
fR <- AAfreq(TCRdat[,1], TCRdat[,2], as.numeric(TCRdat[,3]))

# Compute the homology matrix
S <- seqhom(TCRdat[,1], TCRdat[,2], as.numeric(TCRdat[,3]), 'BLOSUM62')

# Perform score test for continuous outcome
TCRL_cont(Y, X, fR, W, S)

detach(example.data)
```

## Citation

If you use the **tcrl** package in your research, please cite:

Liu, M., Goo, J., Liu, Y., Sun, W., Wu, M.C., Hsu, L. and He, Q., 2022. *TCR-L: an analysis tool for evaluating the association between the T-cell receptor repertoire and clinical phenotypes*. BMC Bioinformatics, 23(1), p.152.

## License

This package is licensed under the LGPL-2.0 License. See the [LICENSE](LICENSE) file for details.

## Authors

- **Meiling Liu** - [mliu@fredhutch.org](mailto:mliu@fredhutch.org)
