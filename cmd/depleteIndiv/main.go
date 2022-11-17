// CLI tool to randomly deplete genotyped markers between min and max marker counts
// of each individuals in the data set.
// Only used for simulating partially overlapped data from fully typed data

package main

import (
    "fmt"
    "os"
    "strconv"
    "math/rand"
    "flag"
    "time"
    "path"
    "github.com/zmaroti/correctKin"
)

// matrix that holds genotype data of samples for the markers
var GENOTYPES [][]uint8

var sampleCount int
    
func depleteIndivs(minMarkerCount, maxMarkerCount int, samples []string) {
    // get a list of marker indexes where all deplete IDs are typed (thats we can go with random depleting)
    markerCount := len(GENOTYPES)

    // check wich markers are fully typed in the selected individuals
    for sampleIdx, sampleName := range samples {
        typedMarkerIdxs := make([]int, 0, markerCount)

        for markerIdx, dat := range GENOTYPES {
            if dat[sampleIdx] != 3 {
                typedMarkerIdxs = append(typedMarkerIdxs, markerIdx)
            }
        }

        // random shuffle genotyped SNP indexes
        rand.Shuffle(len(typedMarkerIdxs), func(i, j int) { typedMarkerIdxs[i], typedMarkerIdxs[j] = typedMarkerIdxs[j],typedMarkerIdxs[i] })

        typedCount := len(typedMarkerIdxs)

        if minMarkerCount > typedCount {
            fmt.Fprintf(os.Stderr, "Less typed markers (%d) than minMarkerCount (%d) in %s\n", typedCount, minMarkerCount,  sampleName)
        } else {
            if maxMarkerCount == 0 {
                maxMarkerCount = typedCount
            } else if maxMarkerCount > typedCount {
                fmt.Fprintf(os.Stderr, "Less typed markers (%d) than maxMarkerCount (%d) in %s\n", typedCount, maxMarkerCount, sampleName)
                maxMarkerCount = typedCount
            }
            keep := minMarkerCount + rand.Intn(maxMarkerCount - minMarkerCount)

            if keep < typedCount {
                for i := keep; i < typedCount; i++ {
                    markerIdx := typedMarkerIdxs[i]

                    GENOTYPES[markerIdx][sampleIdx] = 3
                }
            }
        }

    }
}

func printHelp() {
    fmt.Fprintln(os.Stderr,
`USAGE
depleteIndivs <DATA.(bed,geno)> <minMarker> [maxMarker]

The tool will deplete ALL individual sample's SNPs to a random marker count between minMarker and the maxMarker SNP counts in the PLINK (.bed, .fam, .bim), EIGENSTRAT or  PACKEDANCESTRYMAP (.geno, .snp, .ind) dataset and writes the output into the appropriate OUTPREFIX.(bed, .geno) file. If maxMarker count is not provided then the maximum is the SNP count of the given individual.

Note, only the appropriate data file is created while the marker (.bim, .snp) and family data (.fam, .ind) files are not copied to the outprefix, these have to be copied manually.

optional flags:
 -help         print this help
 -out  value   output file prefix to use for the output (depleted) dataset DEFAULT: DATA_dep
 -seed value   The random seed to use, DEFAULT: generate from unix time
`)

    os.Exit(0)
}

func main() {
    var help bool
    var seed int
    var outPref string
    var samples []string
    var minMarker, maxMarker int
    var err error

    flag.BoolVar(&help,      "help", false, "print help")
    flag.IntVar(&seed,       "seed", 0, "The random seed to use, DEFAULT: generate from unix time")
    flag.StringVar(&outPref, "out", "", "output file prefix to use for the randomly pseudo haploized dataset DEFAULT: DATA_PREFIX_dep")

    flag.Parse()

    args := flag.Args()

    // missing/too much arguments
    if help || (len(args) < 2 || len(args) > 3) {
        printHelp()
    }

    // parse min/max marker counts
    minMarker, err = strconv.Atoi(args[1])
    if err != nil {
        panic(err)
    }

    if len(args) == 3 {
        maxMarker, err = strconv.Atoi(args[2])
        if err != nil {
            panic(err)
        }

        if maxMarker < minMarker {
            fmt.Fprintf(os.Stderr, "minMarker (%d) must less than equal to maxMarker (%d) count.\n", minMarker, maxMarker)
            os.Exit(1)
        }
    }

    // seed random generator
    if seed == 0 {
        rand.Seed(time.Now().UTC().UnixNano())
    } else {
        rand.Seed(int64(seed))
    }

    fext := path.Ext(args[0])


// figure out prefix from provided filename
    prefix := args[0][0:len(args[0]) - len(fext)]

    // read sample IDs from data sets
    if fext == ".geno" {
        samples     = correctKin.ReadIND(prefix + ".ind")
    } else if fext == ".bed" {
        samples     = correctKin.ReadFAM(prefix + ".fam")
    } else {
    // not the data files provided?
        fmt.Fprintf(os.Stderr, "Either a '.bed' or a '.geno' file has to be provided as input (%s)\n", args[0])
        os.Exit(1)
    }

    // have default out prefix if omitted
    if outPref == "" {
        outPref = prefix + "_dep"
    // have a warning that outpref is same as prefix, so file will be overwritten
    } else if outPref == prefix {
        fmt.Fprintf(os.Stderr, "WARNING: output file will overwrite the input file %s\n", args[0])
    }

    var header []byte
    var markerCount int

    // read actual big data files
    // we have either .geno or .bed file extensions (previously checked)
    if fext == ".geno" {
        markerCount = correctKin.LineCount(prefix + ".snp")

        sampleCount = len(samples)

        header, GENOTYPES  = correctKin.ReadEIG(args[0], sampleCount, markerCount)
    } else {
        markerCount = correctKin.LineCount(prefix + ".bim")

        sampleCount = len(samples)

        header, GENOTYPES  = correctKin.ReadBED(args[0], sampleCount, markerCount)
    }

    depleteIndivs(minMarker, maxMarker, samples)

    // write out datasets
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
