// CLI tool to test robustness of the methodolgy against aDNS specific genotyping errors
// NOT REQUIRED IN THE correctKin workflow, only used to perform simulations in the manuscript 
package main

import (
    "fmt"
    "os"
    "math/rand"
    "flag"
    "time"
    "path"
    "github.com/zmaroti/correctKin"
)

// matrix that holds genotype data of samples for the markers
var GENOTYPES [][]uint8

// this struct holds the marker information of EIG or PLINK data
var SNPs   []correctKin.SNP

var sampleCount int

// simulate post mortem damage (C->T or G->A conversions)
// and REF bias (exogenous DNA contamination)
// in PLINK data set
func simContam(conv, ref, cont float64, endoGENOTYPES [][]uint8) {
    // in the 2D matrix of GENOTYPES
    // we have the different markers in the rows (indexed by i) and data for different samples (indexded by j) in the columns

    // GENOTYPES holds data in EIGENSTRAT internal format
    // 9 - missing
    // 0 - HOM_ALT (0 REF ALLELE)
    // 1 - HET     (1 REF ALLELE)
    // 2 - HOM_REF (2 REF ALLELES)

    // first pass, check if we have NO HET alleles (1) so we are pseudo haplopid (9, 0, 2)
    for _, dat := range GENOTYPES {
        for _, geno := range dat {
            // EIG style HET_ALT = 1 (1 REF ALLELE)
            // only haploid data sets are considered, so exit with error
            if geno == 1 {
                fmt.Fprintln(os.Stderr, "HET allele in the dataset, this tool is intended to used for haploid data onls. EXITING!")
                os.Exit(1)
            }
        }
    }

    // in case we have reference bias error, we have to randomly change given fraction of
    // markers to HOM REF state in each samples regardless of its previous state (missing, ref or alt)
    if ref > 0.0 {
        // create a list of marker indexes to shuffle randomly
        allMarkerIdxs := make([]int, len(GENOTYPES))

        for i, _ := range allMarkerIdxs {
            allMarkerIdxs[i] = i
        }

        // count of markers to change
        count := int(float64(len(allMarkerIdxs)) * ref)

        // set to HOM REF a random set of markers for each samples
        for sampleIdx := 0; sampleIdx < sampleCount; sampleIdx++ {
            // shuffle the indexes randomly for each sample
            rand.Shuffle(len(allMarkerIdxs), func(i, j int) {
                allMarkerIdxs[i], allMarkerIdxs[j] = allMarkerIdxs[j],allMarkerIdxs[i]
            })

            errCount := rand.Intn(count)

            // set random markers GT to HOM REF
            for c := 0; c < errCount; c++ {
                 GENOTYPES[allMarkerIdxs[c]][sampleIdx] = 0
            }
        }
    }

    // in case we have endogene contamination (human vs human)
    // we flip the GT for a fraction of markers as we have random allele calling
    if cont > 0.0 {
        // check that we don'T have HET calls in contaminant samples
        for _, dat := range endoGENOTYPES {
            for _, geno := range dat {
                // EIG style HET_ALT = 1 (1 REF ALLELE)
                // only haploid data sets are considered, so exit with error
                if geno == 1 {
                    fmt.Fprintln(os.Stderr, "HET allele in the dataset, this tool is intended to used for haploid data onls. EXITING!")
                    os.Exit(1)
                }
            }
        }

        // check if marker set is the same length
        if len(GENOTYPES) != len(endoGENOTYPES) {
            fmt.Fprintln(os.Stderr, "HET allele in the contamination dataset, this tool is intended to used for haploid data onls. EXITING!")
            os.Exit(1)
            
        }

        // create a list of marker indexes to shuffle randomly
        allMarkerIdxs := make([]int, len(GENOTYPES))

        for i, _ := range allMarkerIdxs {
            allMarkerIdxs[i] = i
        }

        // count of markers to change
        count := int(float64(len(allMarkerIdxs)) * cont)

        // set to HOM REF a random set of markers for each samples
        for sampleIdx := 0; sampleIdx < sampleCount; sampleIdx++ {
            // randomize a contaminant from the contaminants
            contamSampleIdx := rand.Intn(len(endoGENOTYPES[0]))

            // shuffle the indexes randomly for each sample
            rand.Shuffle(len(allMarkerIdxs), func(i, j int) {
                allMarkerIdxs[i], allMarkerIdxs[j] = allMarkerIdxs[j],allMarkerIdxs[i]
            })

            errCount := rand.Intn(count)

            // change random markers to contamination sample state
            for c := 0; c < errCount; c++ {
                GENOTYPES[allMarkerIdxs[c]][sampleIdx] = endoGENOTYPES[allMarkerIdxs[c]][contamSampleIdx]
            }
        }
    }
 
    // in case we have conversion error rate, introduce this error
    if conv > 0.0 {
        // check the markers which can have conversion at all (ie have C or G as minor or major allele)
        convMarkerIdxs := make([]int, 0, len(SNPs))
        
        for idx, snp := range SNPs {
            if snp.REF == "C" || snp.REF == "G" || snp.ALT == "C" || snp.ALT == "G" {
                convMarkerIdxs = append(convMarkerIdxs, idx)
            }
        }

        count := int(float64(len(convMarkerIdxs)) * conv)

        // for each sample randomize the list of conversable markers
        // and apply PMD
        for sampleIdx := 0; sampleIdx < sampleCount; sampleIdx++ {
            // shuffle the indexes randomly for each sample
            rand.Shuffle(len(convMarkerIdxs), func(i, j int) {
                convMarkerIdxs[i], convMarkerIdxs[j] = convMarkerIdxs[j],convMarkerIdxs[i]
            })

            errCount := rand.Intn(count)

            // apply PMD on first count of random marker positions for the sample
            for c := 0; c < errCount; c++ {
                sampleEigGT := GENOTYPES[convMarkerIdxs[c]][sampleIdx]

                // sample is HOM REF
                if sampleEigGT == 0 {
                     // sample allele is hom C, that is converted to hom T
                    if SNPs[convMarkerIdxs[c]].REF == "C" {
                        // if marker ALT allele is T then sample flipped from C to T
                        if SNPs[convMarkerIdxs[c]].ALT == "T" {
                            GENOTYPES[convMarkerIdxs[c]][sampleIdx] = 2 // that is HOM ALT
                        } else {
                        // otherwise sample GT does not conforms with REF or ALT so it is missing
                            GENOTYPES[convMarkerIdxs[c]][sampleIdx] = 9
                        }
                     // sample allele is hom G, that is converted to hom A
                    } else if SNPs[convMarkerIdxs[c]].REF == "G" {
                        // if marker ALT allele is T then sample flipped from C to T
                        if SNPs[convMarkerIdxs[c]].ALT == "A" {
                            GENOTYPES[convMarkerIdxs[c]][sampleIdx] = 2 // that is HOM ALT
                        } else {
                        // otherwise sample GT does not conforms with REF or ALT so it is missing
                            GENOTYPES[convMarkerIdxs[c]][sampleIdx] = 9
                        }
                    }
                    // else nothing to do, our sample had no alleles that can suffer PMD at this locus
                // sample is HOM ALT
                } else if sampleEigGT == 2 {
                     // sample allele is hom C, that is converted to hom T
                    if SNPs[convMarkerIdxs[c]].ALT == "C" {
                        // if marker REF allele is T then sample flipped from C to T
                        if SNPs[convMarkerIdxs[c]].REF == "T" {
                            GENOTYPES[convMarkerIdxs[c]][sampleIdx] = 0 // that is HOM REF
                        } else {
                        // otherwise sample GT does not conforms with REF or ALT so it is missing
                            GENOTYPES[convMarkerIdxs[c]][sampleIdx] = 9
                        }
                     // sample allele is hom G, that is converted to hom A
                    } else if SNPs[convMarkerIdxs[c]].ALT == "G" {
                        // if marker ALT allele is T then sample flipped from C to T
                        if SNPs[convMarkerIdxs[c]].REF == "A" {
                            GENOTYPES[convMarkerIdxs[c]][sampleIdx] = 0 // that is HOM REF
                        } else {
                        // otherwise sample GT does not conforms with REF or ALT so it is missing
                            GENOTYPES[convMarkerIdxs[c]][sampleIdx] = 9
                        }
                    }
                    // else nothing to do
                }
                // even though this site could be convesed theoretically if sample is untyped at
                // this site, we have nothing to do
            }
        }
    }
}


func printHelp() {
    fmt.Fprintln(os.Stderr,
`USAGE
simContam [-seed NUMBER|default random] [-out OUT_PREFIX|Default DATAPREFIX_err] [-conv conversion error fraction|Default 0.05] [-ref reference error fraction|Default 0.05] <DATAPREFIX.(bed,geno)> 

The tool will introduce random error in a HAPLOID binary PLINK (.bed, .fam, .bim), EIGENSTRAT or PACKEDANCESTRYMAP (.geno, .snp, .ind) dataset and write the output into the appropriate OUTPREFIX.(bed, .geno) file.
NOTE: this tool is written for applying on haploid data set (see bellow the limitations of PLINK/EIG data sets). Simulating diploid data for such error would require additional information, what was the average coverage what is the overal distribution of coverage and whether at a given coverage changing one read would change the infered genotype (this is also genotype caller dependent). However when we have random alele call pseudo haploid gnotypes from low coverage aDNA we can assume the worst scenario, that is we have only 1 read at the given marker position and any error (PMD or ref bias) would change the genotype (see bellow for constraints).

IMPORTANT CONSTRAINTS OF the DATA FORMATs:
Since PLINK and EIGENSTRAT data formats are bialellic, and cannot store iformation on a second ALT allele like VCF,the samples can have only the following states: HOM MAJOR, HET (MINOR/MAJOR) and HOM MINOR, or MISSING information at each marker positions. Thus any other genotype (be it gnotype error in a given sample or a valid different ALT2 allele) cannot be imported into the data set. Accordingly, when at one position there are more than 2 alleles present importing data into these format has two options: #1) omit the marker alltogether, #2) samples that has different than MINOR/MAJOR alleles are set to MISSING information for the given marker. Thus simulating genotype errors in the PLINK and EIGENSTRAT data can only be performed by either removing markers which have any conflicting genotype in any samples (this is unfeasible for large sample cophort as we could remove all markers even at low genotyping error rate when many samples exists in the dataset) or setting the given sample at the given marker to missing information. Usually this second option is used as random missing informatiom in samples does not alter the result of projected PCA or other population genetic analyses that accounts for missing data points in individuals.

In typical NGS sequencing the technical error rate is very low (typically bmajority of bases have phred score Q30 or higher meaning the chance of having a genotype error is less than 1:1000).
However in case of aDNA post mortem damage (PMD) specific nucleotids have higher chance of conversions (C->T or a G->A). And also bad QC of aligned reads could lead to so called reference bias, when exogenous non human DNA is forced on the human reference genome (when the read has sufficient length of matching DNA) and since non human (evolutionary older) organism have usually the ancestral (usually major allele) state, this could lead to having genotypes with HOM REF state that are the result of aligned exogenous reads.

optional flags:
-seed positive_number   the randomization will be performed with this seed, DEFAULT (0) generate from unix time
-out  OUT_PREFIX        outprefix of the output, if omited the default is DATAPREFIX_err
-conv FLOAT             the percentage of conversion error
-ref  FLOAT             the percentage of reference error
-cont FLOAT             the percentage of flip GT (contamination) error
`)

    os.Exit(0)
}

func main() {
    var help bool
    var seed int
    var outPref, endoFn string
    var conv, ref, cont float64

    flag.BoolVar(&help,      "help", false, "print help")
    flag.IntVar(&seed,       "seed", 0, "The random seed to use, DEFAULT 0 -> generate from unix time")
    flag.StringVar(&outPref, "out", "", "output file prefix to use for the randomly pseudo haploized dataset DEFAULT: DATA_PREFIX_haploid")
    flag.StringVar(&endoFn, "efile", "", "BED file for outlier contaminator") 
    flag.Float64Var(&conv,       "conv", 0.05, "The percentage of conversion error introduced|DEDFAUL 0.05")
    flag.Float64Var(&ref,        "ref",  0.05, "The percentage of reference error introduced|DEFAULT 0.05")
    flag.Float64Var(&cont,       "cont",  0.05, "The percentage of contamination (flip GT) error introduced|DEFAULT 0.05")

    flag.Parse()

    args := flag.Args()

    // missing/too much arguments
    if help || (len(args) != 1) {
        printHelp()
    }

    fext := path.Ext(args[0])

    // not the data genom files provided?
    if fext != ".geno" && fext != ".bed" {
        printHelp()
    }

    // seed random generator
    if seed == 0 {
        rand.Seed(time.Now().UTC().UnixNano())
    } else {
        rand.Seed(int64(seed))
    }

    // figure out prefix from provided filename
    prefix := args[0][0:len(args[0]) - len(fext)]

    // have default out prefix if omitted
    if outPref == "" {
        outPref = prefix + "_haploid"
    // have a warning that outpref is same as prefix, so file will be overwritten
    } else if outPref == prefix {
        fmt.Fprintf(os.Stderr, "Output file will overwrite the input file %s\n", args[0])
    }

    var endoGENOTYPES [][]uint8

    if cont > 0.0 {
        if endoFn != "" {
            endoFext := path.Ext(endoFn)

            // not the data genom files provided?
            if endoFext != ".geno" && fext != ".bed" {
                printHelp()
            }

             // figure out prefix from provided filename
            endoPrefix := endoFn[0:len(endoFn) - len(endoFext)]

            if fext == ".geno" {
                endoSamples  := correctKin.ReadIND(endoPrefix + ".ind")
                endoSNPs := correctKin.ReadSNP(endoPrefix + ".snp")

                // used to validate binary dataset size
                endoSampleCount := len(endoSamples)
                endoMarkerCount := len(endoSNPs)

                _, endoGENOTYPES  = correctKin.ReadEIG(endoFn, endoSampleCount, endoMarkerCount)
            } else {
                endoSamples     := correctKin.ReadFAM(endoPrefix + ".fam")
                endoSNPs := correctKin.ReadBIM(endoPrefix + ".bim")

                // used to validate binary dataset size
                endoSampleCount := len(endoSamples)
                endoMarkerCount := len(endoSNPs)

                _, endoGENOTYPES  = correctKin.ReadBED(endoFn, endoSampleCount, endoMarkerCount)
            }
        } else {
            fmt.Fprintln(os.Stderr, "Endogene contamination file not provieded, exiting")
            os.Exit(1)
        }
    }

    var header []byte
    var samples []string
    var markerCount int

    if fext == ".geno" {
        samples     = correctKin.ReadIND(prefix + ".ind")
        SNPs        = correctKin.ReadSNP(prefix + ".snp")

        // used to validate binary dataset size or line count in flat format
        sampleCount = len(samples)
        markerCount = len(SNPs)

        header, GENOTYPES  = correctKin.ReadEIG(args[0], sampleCount, markerCount)
    } else {
        samples     = correctKin.ReadFAM(prefix + ".fam")
        SNPs = correctKin.ReadBIM(prefix + ".bim")

        // used to validate binary dataset size
        sampleCount = len(samples)
        markerCount = len(SNPs)

        header, GENOTYPES  = correctKin.ReadBED(args[0], sampleCount, markerCount)
    }
  
    simContam(conv, ref, cont, endoGENOTYPES)

    if fext == ".geno" {
        if header == nil {
            correctKin.WriteEIG(outPref + ".geno", GENOTYPES)
        } else {
            correctKin.WritePackedEIG(outPref + ".geno", header, GENOTYPES)
        }
    } else {
        correctKin.WriteBED(outPref + ".bed", header, GENOTYPES)
    }
}
