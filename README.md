correctKin : Kinship coefficient estimation from low coverage ancient (or degraded forensic) genome data.
======================================================
This is a repository containing the tools to infer the degree of relationship from ancient/low coverage DNA. The approach is described detail in our manuscript:
Emil Nyerki, Tibor Kalmár, Oszkár Schütz, Rui M. Lima, Endre Neparáczki, Tibor Török and Zoltán Maróti: "An optimized method to infer relatedness up to the 4th degree from low coverage ancient human genomes"

doi: https://doi.org/10.1186/s13059-023-02882-4

The repository contains a minimal library to read/write binary PLINK, plain text EIGENSTRAT and binary PACKEDANCESTRYMAP data sets and command-line tools to aid importing haploid genotype data, calculate the pairwise marker overlap fraction between samples, calculate the corrected kinship coefficient matrix from the PCangsd numpy kinship coefficient estimation output file and the marker overlap fraction file, and filter the potential relatives from the corrected kinship coefficient matrix. Additionally it contains the tools that were used in our original manuscript to simulate partially genotyped data from fully typed data sets.

### Installation and Building
To install and build, ensure you have `go` and `make` installed. Clone (or `go get`) this repo into your `$GOPATH`:
```sh
go get github.com/zmaroti/correctKin
```

Enter the package directory and type
```sh
make
```
to build all command line tools in the bin directory.

Alternatively you can only build the selected tools manually:

```sh
go build -o bin/cmdName ./cmd/cmdName
```

Copy all the binary files or include the bin directory into your `$PATH` to use the CLI tools from shell.

### Usage
Each command line tool will print a help if issued without the proper options or with the '-help' flag.

```sh
cmdName -help
```

### Quickstart with practical examples
A [PDF](correctKin_quickstart.pdf) is available for detailed description and considerations for analyzing various data sources.

## Description of individual tools used in kinship analysis
# importHaploCall
This tool is to import ANGSD random pseudo haploid call genotype data into the standard binary PLINK or plain text EIGENSTRAT genotype data format.

```sh
importHaploCall <EIGENSTRAT.snp> <outfile.(bed|geno)> <list of ANGSD.haplo.gz>
```
As input it requires an EIGENSTRAT.snp file that contains the required marker informations (CHR, POS, REF, ALT) to generate the genome data. This file must be the same that were used to create the ANGSD sites file. 
An output filename ending with the ('.bed' or '.geno' suffix), and one more more ANGSD haploid call genotype files.

In case more haploid genotype files are provided the tool will read all files considering ANGSD typing was done using a scatter gathering method and the genotype data is split between several files. In case the genotype of a given individual does not match the expected REF or ALT alleles of a marker defined in the '.snp' file the chromosome, position, the expected REF/ALT and the GT of the individual will be written in a log file (indX.log).

NOTE: The tool only creates the '.bed' or '.geno' data file and the corresponding '.snp' or '.bim' files (from the provided .snp file). The corresponding '.ind' or '.fam' file has to be manually prepared based on SEX, and the actual IDs of the individuals as this data is not contained and also that ANGSD will use a general ind0, ind1... naming scheme.

# markerOverlap
```sh
markerOverlap [-threads n|default ALL] [-maf M|default 0.05] <DATA.(bed|geno)>
```

This tool calculates the pairwise genotyping overlap fractions between all samples of a genome dataset. The marker overlap fraction is defined between two samples as: (count of markers genotyped in both samples)/(count of all markers in the dataset) /considering only the markers above the maf threshold for both counts/.

As input it requires either a binary PLINK (.bed, .bim, .bam), an EIGENSTRAT or a PACKEDANCESTRYMAP (.geno, .ind, .snp) dataset. The format is guessed by the provided filetype (either the .bed or .geno).

The output is a plain text file (DATA.overlap) containing the pairwise marker overlap fraction matrix. Only markers above the minor allele frequency threshold (default 0.05) are counted. This maf threshold is the same as the default maf threshold in PCangsd -kinship calculation, where low AF markers pruned prior to PCA analysis. In case kinship estimation is done with a different threshold, use the same threshold for calculating marker overlap fractions. 

optional flags:
-threads value    the number opf threads to use |default ALL
-maf     value    minor allele frequency threshold (markers below the theshold are not counted) |default 0.05 

# pseudoHaploidize
```sh
pseudoHaploidize [-seed NUMBER|default random] [-out OUT_PREFIX] <DATA.(bed,geno)>
```
The tool will random pseudo haploidize a binary PLINK (.bed, .fam, .bim), EIGENSTRAT or PACKEDANCESTRYMAP (.geno, .snp, .ind) dataset and write the output into the appropriate OUTPREFIX.(bed, .geno) file.

NOTE: only the pseudo-haploidized data file is created without the corresponding '.snp', '.ind' files (or '.bim', '.fam' in case of PLINK input). Since family and SNP marker infromation is same in both data sets those can be copied manually from the original diploid data set.

optional flags:
-seed positive_number   the randomization will be performed with this seed, DEFAULT generate from unix time
-out  OUT_PREFIX        outprefix of the output, if omited the default is DATA_haploid

# filterRelates
```sh
filterRelates [OPTIONS] <marker_overlap_fraction.overlap>  <PCangsd kin coeff matrix.npy>
```
filterRelates performs correction of estimated kinship coefficient, writes out statistics and filters the putative kins from the pairwise corrected kinship coefficient matrix.

The input is a PCangsd v.0.99 kinship coefficient matrix (PREFIX.kinship.npy) and the pairwise marker overlap fraction of the samples generated by markerOverlap (PREFIX.overlap). Based on the two matrices it corrects the estimated kinship coefficients by the marker overlap fraction and writes out the corrected kinship coefficient file (PREFIX.corr.tsv) as a flat TAB separated file (upper triangle corrected kinship coeffients, lower triangle the original kinship coefficient estimate).

The tool calculates statistics on the whole pairwise kinship coefficient matrix where presumably the majority of connections between the individuals are unrelated. Based on this statistics and the overlap fractions it estimates the standard deviation of corrected kinship coefficient for different bins of overlap fractions and writes the statistics (PREFIX.stats.tsv) including the overlap fraction range, mean corrected kinship coeff, SD of corrected kinship coeff, N sigma threshold, 95% confidence interval.

And lastly filters close relatives based on sigma threshold and/or kinship coefficient threshold and prints in FLAT or MATRIX format to the STDOUT.

As default it prints the kinship coefficients for the group of relatives in flat format where the columns contains the following information:
    GROUPID                  number of the kin group
    ID1                      sample ID 1
    ID2                      sample ID 2
    uncorrKinCoeff           the estimated kinship coefficient between the two samples
    overlapFrac              overlap marker fraction between the two samples
    corrKinCoeff             the corrected kinship coefficient between the two samples
    N sigma Thresh           the N sigma threshold of the overlap bin (containing ovelapFrac)
    95% conf (lower)         95% confidence interval lower limit of the corrected kinship coeff
    95%conf(upper)           95% confidence interval upper limit of the corrected kinship coeff
    estimated Relatedness    estimated relatedness (with text)

Alternatively it can also print the values grouped by the relatives as a simple kinship coefficient matrix (upper triangle corrected, lower triangle uncorrected coeffs).

optional flags:
-sigma FLOAT
    A positive number representing the how many SD above the mean corrected kinship coefficient is considered a biological difference (default 6 SD)

-thresh FLOAT
    A positive number to use a corrected kinship coefficient threshold for filtering relatives. If omitted then all pairs above N sigma are considered as kins, if provided the above rule still applies however kins bellow this threshold are not filtered (can be used to limit kins to closer relatives)

-matrix bool|DEFAULT false
    The default is false, meaning we have a FLAT output. If this flag is set only the matrix of coeffs between individuals of each group of relatives (upper triangle holds the corrected, lower triangle holds the uncorrected kinship coefficients) is printed

-showself bool|DEFAULT false
    DEFAULT we don't print kiship coefficients "between" SELF. To override this use the -showself flag. In case of matrix ouput this flag affects whether to display the (self) kinship coefficients at the DIAGONAL or not.

-ids ID1,ID2:ID3,ID4,ID5
    when ID groups are explictly provided the kinship coefficient matrix will be only tested for the ID combinations in each ID groups. So ID1 vs ID2, and ID3 vs ID4; ID3 vs ID5 and ID4 vs ID5

## Description of individual tools used in simulating partially genotyped data sets in the original manuscript.
# depleteIndivs
```sh
depleteIndivs [-out outpref] [-seed randomseed] <DATA_PREFIX.(bed,geno)> <minMarker> [maxMarker]
```
CLI tool to randomly deplete genotyped markers between min and max marker counts of each individuals in the data set. NOT REQUIRED FOR KINSHIP ANALYSISű, only used for simulating partially overlapped data from fully typed data in our manuscript.

# depleteMarkers
```sh
depleteMarkers [-out outpref] [-seed randomseed] <DATA_PREFIX.(bed,geno)> <frac1,s1,s2:frac2,s3,s4,s5:...>
```
CLI tool to deplete groups of samples to a specified marker overlap fraction in a binary PLINK (.bed, .fam, .bim), EIGENSTRAT or PACKEDANCESTRYMAP (.geno, .snp, .ind) dataset and writes the output into the appropriate OUTPREFIX.(bed, .geno) file. NOT REQUIRED FOR KINSHIP ANALYSIS, only used for simulating partially overlapped data from fully typed data in our manuscript.

The tool deplete markers between different sample groups to the desired marker overlap fractions in one pass. The last argument codes the sampleID and marker overlap fraction information for each groups with the following syntax:
    groupinfo1:groupinfo2:groupinfo3...

Where groupinfoN is a comma delimited string with the syntax:
    markerOverlapFractionN,sampleId1,sampleId2,...

Each sample groups consists of at least 2 samples. Each sample Id can be included in only one sample group. Marker overlap fraction depletion is done on random markers till shared marker count is equal with the desired marker overlap fraction in all sample Ids of the given sample group. All remaining genotyped markers that are typed in more sampleIds of the sample group will be kept only in one random individual. If not started from fully typed there might be not overlapping markers between sample groups, in that case a warning with the sample IDs and the expected and maximum overlap fraction values will be printed on STDERROR.

Note, only the appropriate data file is created while the marker (.bim, .snp) and family data (.fam, .ind) files are not copied to the outprefix, these have to be copied manually.

