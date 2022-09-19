![logo](https://user-images.githubusercontent.com/57667417/191082345-57fed066-37e9-4a8a-a65a-c9562d0625a4.png)

# Conditional frequency distribution

The (m,n)-mer R package was created to summarize biological data into numerical characteristics. It reads a FASTA file and generates a table describing the conditional frequency distribution of the selected (m,n)-mer in the sequences. This output is combined with class information to generate a feature matrix for classification. 

(m,n)-mers are an alternative case of k-mers (Figure 1). We proposed the replacement of the unconditional k-mer frequency 
by the conditional frequency, which represents the relative frequency of the n-mer conditioned to the relative frequency of m-mer plus n-mer (for more details and performance comparisson please see Andrade et al., 2022 (in press)). 

![Fig1](https://user-images.githubusercontent.com/57667417/191081859-0b0ae464-f257-4c82-9dea-8d4629605357.png)

**Fig 1.** Comparing k-mer to mn-mer relative frequency.

According to Figure 2, the k-mers are represented as (0,k) and the mn-mers as (m,n).

![Fig2](https://user-images.githubusercontent.com/57667417/191081936-1aed5ca6-9c88-4d4d-a46b-e1ccae0bcafe.png)

**Fig 2.** Numeric representation.

The output table (Figure 3) includes the fasta file accession numbers as an ID column, the relative frequency of mn-mers up to 4k columns, and class information. 

![Fig3](https://user-images.githubusercontent.com/57667417/191082016-b6835c4c-c115-498d-a2d1-c7d93ec20fe5.png)
**Fig 3.** Output example.


## Dependencies



## Installation

```
library(devtools)

devtools::install_github("labinfo-lncc/mnmer", ref="main")
```



## Quick Start: Running (m,n)-mer on example dataset

```
library("mnmer")
```

#### Summary

Function to summarize base content, including N and IUPAC bases

#### Filtering sequences

Function to remove sequences with a certain content of non informative bases and length 

#### Producing k-mers

#### Producing (m,n)-mers 

Other features: paralelism, optimize memory use, read fasta in blocks and zipped



## Citation

All data used in this project are stored at: (link) (FASTA files), (link) (Feature Matrices), and (link) (Results Metrics). 



## Work in progress

If you have any queries or find a bug, please submit an issue on GitHub or email atrv@lncc.br.
