# tcrl: Association Analysis of T-Cell Receptor Repertoire and Clinical Phenotypes

[![License](https://img.shields.io/badge/license-LGPL--2.0-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.0.html)

## Overview

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
