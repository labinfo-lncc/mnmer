![image-20220908153252560](/home/amanda/.config/Typora/typora-user-images/image-20220908153252560.png)

# Conditional frequency distribution

The (m,n)-mer R package was created to summarize biological data into numerical characteristics. It reads a FASTA file and generates a table describing the conditional frequency distribution of the selected (m,n)-mer in the sequences. This output is combined with class information to generate a feature matrix for classification. 

(m,n)-mers are an alternative case of k-mers (Figure 1). We proposed the replacement of the unconditional k-mer frequency 
by the conditional frequency, which represents the relative frequency of the n-mer conditioned to the relative frequency of m-mer plus n-mer (for more details and performance comparisson please see Andrade et al., 2022 (in press)). 

<img src="/home/amanda/.config/Typora/typora-user-images/image-20220908154630665.png" alt="image-20220908154630665" style="zoom:67%;" />

**Fig 1.** Comparing k-mer to mn-mer relative frequency.

According to Figure 2, the k-mers are represented as (0,k) and the mn-mers as (m,n).

<img src="/home/amanda/.config/Typora/typora-user-images/image-20220908160359290.png" alt="image-20220908160359290" style="zoom:67%;" />

**Fig 2.** Numeric representation.

The output table (Figure 3) includes the fasta file accession numbers as an ID column, the relative frequency of mn-mers up to 4k columns, and class information. 

<img src="/home/amanda/.config/Typora/typora-user-images/image-20220908153711469.png" alt="image-20220908153711469" style="zoom:67%;" />

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
