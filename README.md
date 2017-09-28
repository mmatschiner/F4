# F4
Assessing confidence in introgression among four populations

#### Summary

F4 calculates the f4-statistic from allele frequencies of four populations and uses coalescent simulations to test whether this value could be the result of incomplete lineage sorting.

#### Background

The f4-statistic was introduced by Reich et al. (2009) and is a powerful measure to distinguish introgression from incomplete lineage sorting, based on allele frequencies of four populations. With populations A, B, C, and D, and the assumed population topology (A,B),(C,D), the f4-statistic is calculated as the product of the difference of allele frequencies between A and B, and between C and D. Thus, if at a particular SNP, the frequency of the base "G" was 0.4 in population A, 1.0 in population B, 0.2 in population C, and 0.0 in population D, the f4-statistic for this SNP would be (0.4-1.0)x(0.2-0.0) = -0.6x0.2 = -0.12. With more than one SNP, the f4-statistic of the whole data set is simply the mean of the f4 values of all individual SNPs. What's interesting about this measure is that under incomplete lineage sorting alone, the allele frequency differences between A and B should be independent of those between C and D, and the f4-statistic should thus be zero. If there is introgression however between the two pairs of populations (e.g. A introgressed into C), this would lead to non-zero f4 values.

The tricky part however is not to calculate the f4-statistic, but to assess support for it being different from zero, and thus for introgression. As implemented in the program fourpop which comes as part of the Treemix package (Pickrell & Pritchard 2012), block jackknifing is one way to do this, which is commonly used. This means that the data set is chopped into blocks of a particular size, and the f4-statistic is calculated individually for each of these blocks. By then comparing the overall f4 value to the standard error of f4 values taken from jackknife blocks, a z-score is calculated and serves as a measure of support. There are a few problems to this approach. One is that SNP data sets, e.g. those obtained by RAD sequencing often contain linked groups of SNPs that can confound the standard error of jackknife block f4 values. Further, the use of z-scores as support assumes that the underlying data is normally distributed, but often, f4 values of jackknife blocks are not. This is due to a large proportion of these values being exactly zero if the block from which they are taken does not show any evidence of either introgression or  incomplete lineage sorting (or worse, if monomorphic SNPs were not removed from the data set, which they should be for this analysis). Another issue is that the chosen value of the jackknife block size can (and often does) influence the standard error and thus the z-score supporting introgression.

The approach chosen by F4 is therefore to run simulations to assess support for introgression. After all, we're not so much interested in whether the observed f4-statistic is different from zero, but in whether or not it could be produced by incomplete lineage sorting alone, without any introgression. F4 uses the coalescent software fastsimcoal2 (Excoffier et al. 2013) to produce SNP data sets that resemble the actual SNP data set, with migration rates set to zero and therefore strictly without introgression. All simulated data sets have the same number of individuals and SNPs, and are masked to include the same amount of missing data per population and SNP as the original data set. Further, a burn-in phase is used to automatically adjust simulation parameters (effective population size and relative divergence times) so that the resulting SNP variation matches the observed, both between all populations and between the two pairs of populations. Thus, the simulated data sets should be equivalent to the original data set, but have been produced completely without introgression. Finally, F4 calculates the f4-statistic for each simulated data set and reports the proportion of these that is more extreme than the observed. If the number of simulations is large enough, this proportion of more extreme values can be taken as the probability of obtaining this f4 value (or a more extreme one) without introgression. Thus if this probability is small enough (e.g. < 0.05), it supports introgression between the two pairs of populations. It does not, however, indicate the directionality of introgression.

#### Requirements for F4
F4 is written in python3, therefore this version of python must be installed. It also uses a number of python packages, most of these however are likely to be installed on your machine anyway, if you have python. Only the two packages numpy and scipy may require additional installations. Importantly, F4 uses fastsimcoal2 for coalescent simulations, thus this program needs to be installed, and needs to be executable with `fsc252`, for F4 to run.

#### F4 input
To keep things simple (and to allow comparison), F4 uses the exact same input format as Treemix. It is described in the Treemix manual that you can download from [here](https://bitbucket.org/nygcresearch/treemix/wiki/Home ), but it also is simple enough to describe: It starts with a header line that lists the four population names, separated by spaces. It is followed by one line per SNP, listing the allele frequencies per population, again separated by spaces. This is how it goes:

```
pop1 pop2 pop3 pop4  
4,2 0,8 5,3 10,0  
0,6 0,8 0,8 3,7  
...
```  

If your data does not consist of (ideally) unlinked SNPs, but of different blocks of linked SNPs, then these can be specified in the input file as above, but simply with additional empty lines to separate the individual blocks. This will not influence the calculations of the observed f4 value, but simulations will be performed accordingly, so that the number and lengths of simulated linkage blocks matches exactly that of the blocks in the input file. Thus, the reported proportion of simulated f4 values that are more extreme than the observed f4 value should still be a reliable measure of introgression, just as it is when the data set consists only of unlinked SNPs.

#### Example file

The example file `example.txt` can be found in the same directory as f4.py. The data come from a RAD sequencing data set of Sarotherodon cichlid fishes, published by Martin et al. (2015), on the Dryad digital repository (file Sarotherodon-pstacks-map3-mind95-geno5.ped in archive Sarotherodon-plink-nexus-STRUCTURE-files.zip, [here](http://datadryad.org/resource/doi:10.5061/dryad.b28p1 )). As you'll see by running this example file, there is good support for introgression in this data set, despite small samples sizes (betwen 1 and 4 individuals per population) and large amounts of missing data for some of the individuals. Some of the SNPs in this data set appear to be linked (as apparently multiple SNPs per RAD tag were used), but the pattern of introgression persists even after downweighting clusters of f4 outlier SNPs.

#### Running F4
To run F4, simply type

`python3 f4.py example.txt`

With default settings, though, F4 will only calculate the f4-statistic of the original data set. To perform block jackknifing and simulations, use

`python3 f4.py -k 10 -s 1000 example.txt`

This will use jackknife blocks of 10 SNPs each, and conduct 1000 coalescent simulations. If simulated data sets should be written, cause you would also like to use them with Treemix, use the `-o` option, followed by the name of a directory into which files should be written:

`python3 f4.py -k 10 -s 1000 -o temp example.txt`

Finally, the `-l` option, followed by a file name (and only in combination with the `-s` option) will write a log file that contains for each simulation the effective population size parameter and the time of the second divergence used for the simulation, as well as the resulting proportion of SNPs in the simulated data set that are variable in more than one population, and the proportion of SNPs that are variable on both sides of the root split. The first two parameters (effective population size and time of the second divergence) can't easily be translated into anything biologically meaningful since the time of the first divergence event (the root) is arbitrarily fixed to 1000 time units. Nevertheless, changing these two parameters also changes the resulting proportions of variable SNPs, and thus, the first two parameters can be iteratively optimized to produce a data set that matches the observed data set in terms of SNP variability. And this optimization is exactly what's happening automatically during the burnin phase of the simulations. The log file allows you to assess convergence of these parameters, as well as of the simulated f4 values. It is written in a format that can be read by [Tracer](http://tree.bio.ed.ac.uk/software/tracer/).

#### F4 output
Except when the `-o` or `-l` options are used, F4 does not write any files, but only reports to stdout.

#### How to cite F4
F4 is described and applied in our publication on introgression between tribes of the East African cichlid radiation:

[Meyer BS, Matschiner M, Salzburger W (2016) Disentangling incomplete lineage sorting and introgression to refine species-tree estimates for Lake Tanganyika cichlid fishes. Systematic Biology, 66, 531-550.](https://academic.oup.com/sysbio/article-abstract/66/4/531/2670093/Disentangling-Incomplete-Lineage-Sorting-and)

#### References
* Reich D, Thangaraj K, Patterson N, Price AL, Singh L (2009) Reconstructing Indian population history. Nature 461, 489-494.
* Pickrell JK, Pritchard JK (2012) Inference of population splits and mixtures from genome-wide allele frequency data. PLoS Genetics 8, e1002967.
* Excoffier L, Dupanloup I, Huerta-Sañchez E, Sousa VC, Foll M (2013) Robust demographic inference from genomic and SNP data. PLoS Genetics 9, e1003905.
* Martin CH, Cutler JS, Friel JP, Dening T C, Coop G, Wainwright PC (2015) Complex histories of repeated gene flow in Cameroon crater lake cichlids cast doubt on one of the clearest examples of sympatric speciation. Evolution, 69, 1406-1422.