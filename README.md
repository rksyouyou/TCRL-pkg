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

detach(example.data)

# View the result
head(freq_matrix[,1:3])
            [,1]        [,2]       [,3]
 [1,] 0.11646586 0.080321285 0.01606426
 [2,] 0.08227848 0.008438819 0.05485232
 [3,] 0.07883817 0.029045643 0.01244813
 [4,] 0.10243902 0.039024390 0.04390244
 [5,] 0.09810127 0.069620253 0.03164557
 [6,] 0.09523810 0.000000000 0.04761905
```
The `AAfreq` function computes amino acid frequencies within the TCR repertoire across multiple subjects. This function is essential for understanding the distribution and abundance of amino acids, which plays a crucial role in evaluating immune diversity and response. By inputting subject-specific TCR sequence data, users can generate a frequency matrix for further analysis, making it particularly useful in immunology and bioinformatics studies that aim to characterize TCR diversity and its associations with clinical phenotypes.

The function returns a **numeric matrix**, where each **row** represents a subject, each **column** corresponds to an amino acid. The values indicate the frequency of each amino acid in the subject’s TCR repertoire. This output matrix can be used for further analyses, including evaluating amino acid diversity patterns across a study cohort or investigating potential links between TCR repertoire composition and clinical outcomes.


### Example 2: TCR Repertoire Homology Calculation

```{r}
# Load example dataset
data("example.data")
attach(example.data)

# Compute the homology matrix using BLOSUM62
S <- seqhom(TCRdat[,1], TCRdat[,2], as.numeric(TCRdat[,3]), "BLOSUM62")

# View the homology matrix
print(S[1:3,1:3])                                                                                                                                                                         
           [,1]      [,2]      [,3]
 [1,] 1.0000000 0.2360752 0.2830869
 [2,] 0.2360752 1.0000000 0.3423294
 [3,] 0.2830869 0.3423294 1.0000000

detach(example.data)
```

The `seqhom` function quantifies the homology between TCR repertoires across different subjects based on amino acid sequence similarities. It utilizes substitution matrices, such as BLOSUM and PAM, to evaluate sequence similarity, providing a numerical measure of TCR repertoire homology. The function returns a symmetric homology matrix (S), where each entry S[i,j] represents the computed homology score between subjects i and j. Higher values indicate greater similarity between TCR repertoires. This method provides a powerful framework for analyzing immune repertoire similarities and investigating shared immune responses across individuals.


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

fixed.effect  rand.effect Overall.pval
   0.3738187    0.2862339    0.3461362

detach(example.data)
```

The **TCRL_cont** function evaluates the association between T-cell receptor (TCR) repertoire features and a **continuous clinical phenotype**, incorporating both **fixed and random effects** to account for variability across subjects. In this analysis, the **fixed effect p-value** (0.3738) suggests that TCR features do not have a consistent association with the phenotype across all subjects. Similarly, the **random effect p-value** (0.2862) indicates that allowing the effect to vary between individuals does not reveal a statistically significant relationship. The **overall p-value** (0.3461), obtained by combining the fixed and random effects using **Fisher’s procedure**, further supports this conclusion, as it is well above the standard significance threshold of 0.05. These results suggest that, in this example dataset, there is no strong evidence of an association between TCR repertoire features and the continuous clinical phenotype. 

## Citation

If you use the **tcrl** package in your research, please cite:

Liu, M., Goo, J., Liu, Y., Sun, W., Wu, M.C., Hsu, L. and He, Q., 2022. *TCR-L: an analysis tool for evaluating the association between the T-cell receptor repertoire and clinical phenotypes*. BMC Bioinformatics, 23(1), p.152.

## License

This package is licensed under the LGPL-2.0 License. See the [LICENSE](LICENSE) file for details.

## Authors

- **Meiling Liu** - [mliu@fredhutch.org](mailto:mliu@fredhutch.org)
