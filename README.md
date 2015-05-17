# F4

Assessing confidence in introgression among four populations

#### Summary

F4 calculates the f4-statistic from allele frequencies of four populations and uses coalescent simulations to test whether this value could be the result of incomplete lineage sorting.

#### Introduction

The f4-statistic was introduced by Reich et al. (2009) and is a powerful measure to distinguish introgression from incomplete lineage sorting based on allele frequencies of four populations. With populations A, B, C, and D, and the assumed population topology (A,B),(C,D), the f4-statistic is calculated as the product of the difference of allele frequencies between A and B, and between C and D. Thus, if at a particular SNP, the frequency of the base "G" was 0.4 in population A, 1.0 in population B, 0.2 in population C, and 0.0 in population D, the f4-statistic for this SNP would be (0.4-1.0)x(0.2-0.0) = -0.6x0.2 = -0.12. With more than one SNP, the f4-statistic of the whole data set is simply the mean of the f4 values of all individual SNPs. What's interesting about this measure is that under incomplete lineage sorting alone, the allele frequency differences between A and B should be unrelated to those between C and D, and the f4-statistic should thus be 0. If there is introgression however between the two pairs of populations (e.g. A introgressed into C), this would lead to non-zero f4 values.