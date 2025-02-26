# tcrl: Association Analysis of T-Cell Receptor Repertoire and Clinical Phenotypes

[![License](https://img.shields.io/badge/license-LGPL--2.0-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.0.html)

## Overview

T cell receptors (TCRs) are essential for adaptive immunity, and advances in genome technology now enable analysis of the TCR repertoire at the individual sequence level. Studying TCR repertoire associations with clinical phenotypes can provide insights into immune-mediated diseases, yet methods for such analyses remain underdeveloped. We introduce **TCR-L**, a statistical tool for evaluating associations between the TCR repertoire and disease outcomes. Built on a mixed-effect model, **TCR-L** accounts for both **fixed effects**, representing explicitly extracted TCR sequence features, and **random effects**, capturing hidden sequence characteristics. Independent statistical tests assess these effects separately before combining their results.

Simulation studies confirm that **TCR-L** effectively controls type I error and outperforms methods considering only fixed or random effects ([https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04690-2](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04690-2)). Application to real-world data from a skin cutaneous melanoma study reveals a significant association between TCR repertoire and patient survival. By incorporating both observable and latent TCR features, **TCR-L** offers a robust approach for investigating TCR repertoire associations with disease outcomes.

The **tcrl** package provides a suite of functions for analyzing associations between the T-cell receptor (TCR) repertoire and clinical phenotypes. It is designed to help researchers evaluate the diversity and composition of the TCR repertoire and its relationship to clinical outcomes, such as disease states or treatment responses. The package includes methods for both continuous and binary phenotypes, incorporating homology matrices to account for correlations between subjects' TCR repertoires.

### Features:
- Functions for calculating amino acid frequencies in TCR repertoires.
- Score-based tests for association between TCR repertoire features and clinical phenotypes.
- Support for both continuous and binary outcomes.
- Tools for calculating TCR repertoire homology between subjects using substitution matrices such as BLOSUM and PAM.

## Installation

You can install the development version of **tcrl** from GitHub using the following commands:

```{r}
# Install devtools if necessary
install.packages("devtools")

# Install tcrl from GitHub
devtools::install_github("yourusername/tcrl")
```

## Usage

### Example 1: Amino Acid Frequency Calculation

```{r}
# Load example dataset
data("example.data")
attach(example.data)

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
