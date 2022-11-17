// CLI tool to pseudo-haploidize (randomly change HET alleles to either HOM REF
// or HOM ALT) diploid genome data from binary PLINK, PACKEDANCESTRYMAP or plain
// text EIGENSTRAT data.
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

var sampleCount int

func pseudoHaploidize() {
    for i, dat := range GENOTYPES {
        for j, geno := range dat {
            // EIG style HET_ALT = 1 (1 REF ALLELE)
            // we need to randomize this to EIG style HOM_ALT (0) or HOM_REF (2)
            if geno == 1 {
                if rand.Intn(2) == 1 {
                    GENOTYPES[i][j] = 0
                } else {
                    GENOTYPES[i][j] = 2
                }
            }
        }
    }
}


func printHelp() {
    fmt.Fprintln(os.Stderr,
`USAGE
pseudoHaploidize [-seed NUMBER|default random] [-out OUT_PREFIX] <DATA.(bed,geno)>

The tool will random pseudo haploidize a binary PLINK (.bed, .fam, .bim), EIGENSTRAT or PACKEDANCESTRYMAP (.geno, .snp, .ind) dataset and write the output into the appropriate OUTPREFIX.(bed, .geno) file.

NOTE: only the pseudo-haploidized data file is created without the corresponding '.snp', '.ind' or '.bim', '.fam' files, those files has to be copied manually.

optional flags:
-seed positive_number   the randomization will be performed with this seed, DEFAULT generate from unix time
-out  OUT_PREFIX        outprefix of the output, if omited the default is DATA_haploid
`)

    os.Exit(0)
}

func main() {
    var help bool
    var seed int
    var outPref string

    flag.BoolVar(&help,      "help", false, "print help")
    flag.IntVar(&seed,       "seed", 0, "The random seed to use, DEFAULT: generate from unix time")
    flag.StringVar(&outPref, "out", "", "output file prefix to use for the randomly pseudo haploized dataset DEFAULT: DATA_PREFIX_haploid")

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

    var header []byte
    var samples []string
    var markerCount int

    if fext == ".geno" {
        samples     = correctKin.ReadIND(prefix + ".ind")
        markerCount = correctKin.LineCount(prefix + ".snp")

        sampleCount = len(samples)

        header, GENOTYPES  = correctKin.ReadEIG(args[0], sampleCount, markerCount)
    } else {
        samples     = correctKin.ReadFAM(prefix + ".fam")
        markerCount = correctKin.LineCount(prefix + ".bim")

        sampleCount = len(samples)

        header, GENOTYPES  = correctKin.ReadBED(args[0], sampleCount, markerCount)
    }
  
    pseudoHaploidize()

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
