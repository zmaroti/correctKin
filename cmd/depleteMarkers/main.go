// CLI tool to randomly deplete genotyped markers between group of individuals so
// they share only a desired marker fraction of the whole data set.
// Only used for simulating partially overlapped data from fully typed data

package main

import (
    "fmt"
    "os"
    "strings"
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

type depInfo struct {
    frac float64
    ids  []string
}
    
func depleteMarkers(dinfo depInfo, samples []string) {
    depIdxs := make([]int, len(dinfo.ids), len(dinfo.ids))

    // find the indexes of samples to deplete
    // previously we checked that all depIds are contained in samples
    for i, depId := range dinfo.ids {
        for j, sampleId := range samples {
            if depId == sampleId {
                depIdxs[i] = j
                break
            }
        }
    }

    // get a list of marker indexes where all deplete IDs are typed (thats we can go with random depleting)
    markerCount := len(GENOTYPES)

    // indexes that are fully typed in selected samples and can be potentially depleted
    fullIdxs := make([]int, 0, markerCount)

    // indexes that are not fully typed in the selected samples
    var nonfullIdxs  []int

    // check wich markers are fully typed in the selected individuals
    for idx, dat := range GENOTYPES {
        var missing bool = false

        for _, depIdx := range depIdxs {
            if dat[depIdx] == 3 {
                missing = true
                break
            }
        }

        if missing == false {
            fullIdxs = append(fullIdxs, idx)
        } else {
            nonfullIdxs = append(nonfullIdxs, idx)
        }
    }
    
    // random shuffle indexes

    rand.Shuffle(len(fullIdxs), func(i, j int) { fullIdxs[i], fullIdxs[j] = fullIdxs[j], fullIdxs[i] })
    
    keepCount  := int(float64(markerCount) * dinfo.frac)
    typedCount := len(fullIdxs)

    // check if we have less markers than needed
    if keepCount > typedCount {
        maxFrac := float64(typedCount) / float64(markerCount)

        fmt.Fprintf(os.Stderr, "Not enough overlapping markers in the selected individuals (%s) to have %f overlap fraction max fraction is %f (%d/%d)\n",
            strings.Join(dinfo.ids, ", "), dinfo.frac, maxFrac, typedCount, markerCount)

        os.Exit(1)
    } else {
        // deplete excess fully typed markers -> keep typing in one random sample and set it to missing in other selected
        for i := keepCount; i < typedCount; i++ {
            keepidx := rand.Intn(len(depIdxs))

            for j, sampleIdx := range depIdxs {
                if j != keepidx {
                    GENOTYPES[fullIdxs[i]][sampleIdx] = 3 // set to missing genotype
                }
            }
        }
    }

    // random deplete non fully typed markers as well for selected individuals
    for _, markerIdx := range nonfullIdxs {
        // get the list of samples where we have genotype
        var typedSampleIdxs []int

        for _, sampleIdx := range depIdxs {
            if GENOTYPES[markerIdx][sampleIdx] != 3 {
                typedSampleIdxs = append(typedSampleIdxs, sampleIdx)
            }
        }
        // if more than one typed, randomize one to keep, delete the rest
        if len(typedSampleIdxs) > 1 {
            keepidx := rand.Intn(len(typedSampleIdxs))

            for j, sampleIdx := range typedSampleIdxs {
                if j != keepidx {
                    GENOTYPES[markerIdx][sampleIdx] = 3 // set to missing genotype
                }
            }
        }
    }
}

// parses the input format and validates that it contains ok float fractions
// and sample IDs that exist in the ind/fam files
// the format is 
//        group1;group2;group3
// where groupN is a comma separated list of
//        frac1,sampleid1,sampleid2,...
// denoting the IDs and the marker overlap fractions to deplete to
func parseDepString(depstring string, samples []string) []depInfo {
    var depInfos []depInfo

    // hash to track all sample IDs so we don't have duplicate sample IDs
    selectedSampleIds := make(map[string]bool, len(samples))

    // hash for easy check on whether a sample is in the dataset's sampleIDs
    allSampleIds := make(map[string]bool, len(samples))

    for _, sampleId := range samples {
        allSampleIds[sampleId] = true
    }

    depGroups := strings.Split(depstring, ":")

    for _, groupString := range depGroups {
        depData := strings.Split(groupString, ",")

        // we expect a fraction (float64) and minimum 2 sample Ids
        if len(depData) < 3 {
            fmt.Fprintf(os.Stderr, "Invalid syntax (%s) we expect a comma separated list of a fraction (float64) and minimum two sample Ids\n",
                groupString)

            os.Exit(1)
        }

        // parse fraction
        depfrac, err := strconv.ParseFloat(depData[0], 64)

        if err != nil {
            fmt.Fprintln(os.Stderr, err)

            os.Exit(1)
        }

        if depfrac < 0.0 || depfrac > 1.0 {
            fmt.Fprintf(os.Stderr, "Invalid fraction at (%s), it must be (0<=frac<=1)\n",
                groupString)

            os.Exit(1)
        }

        selectedIds := make([]string, 0, len(depData) - 1)

        // check if all sample IDs are unique and in the flobal samplelist
        for i := 1; i < len(depData); i++ {
            selectedId := depData[i]

            if _, ok := allSampleIds[selectedId]; !ok {
                fmt.Fprintf(os.Stderr, "Sample ID at deplete samples (%s) is not found in the dataset\n",
                    selectedId)

                os.Exit(1)
            }

            if _, ok := selectedSampleIds[selectedId]; ok {
                fmt.Fprintf(os.Stderr, "Sample ID at deplete samples (%s) is duplicated\n",
                    selectedId)

                os.Exit(1)
            }
            // update tracing hash
            selectedSampleIds[selectedId] = true

            // make an id list for the struct
            selectedIds = append(selectedIds, selectedId)
        }

        // append data on the group
        depInfos = append(depInfos, depInfo{frac: depfrac, ids: selectedIds} )
    }

    return depInfos
}

func printHelp() {
    fmt.Fprintln(os.Stderr,
`USAGE
depleteMarkers <DATA.(bed,geno)> <frac1,s1,s2:frac2,s3,s4,s5:...>

The tool will deplete groups of samples to a specified marker overlap fraction in a binary PLINK (.bed, .fam, .bim), EIGENSTRAT or PACKEDANCESTRYMAP (.geno, .snp, .ind) dataset and writes the output into the appropriate OUTPREFIX.(bed, .geno) file.

The command can deplete markers between different sample groups to the desired marker overlap fractions in one pass. The last argument codes the sampleID and marker overlap fraction information for each groups with the following syntax:
    groupinfo1:groupinfo2:groupinfo3...

Where groupinfoN is a comma delimited string with the syntax:
    markerOverlapFractionN,sampleId1,sampleId2,...

Each sample groups consists of at least 2 samples. Each sample Id can be included in only one sample group. Marker overlap fraction depletion is done on random markers till shared marker count is equal with the desired marker overlap fraction in all sample Ids of the given sample group. All remaining genotyped markers that are typed in more sampleIds of the sample group will be kept only in one random individual. If not started from fully typed there might be not overlapping markers between sample groups, in that case a warning with the sample IDs and the expected and maximum overlap fraction values will be printed on STDERROR.

Note, only the appropriate data file is created while the marker (.bim, .snp) and family data (.fam, .ind) files are not copied to the outprefix, these have to be copied manually.

optional flags:
 -help         print this help
 -out  value   output file prefix to use for the output (depleted) dataset DEFAULT: DATA_dep
 -seed value   The random seed to use, DEFAULT: generate from unix time
`)

    os.Exit(0)
}

func main() {
    var help    bool
    var seed    int
    var outPref string

    flag.BoolVar(&help,      "help", false, "print help")
    flag.IntVar(&seed,       "seed", 0, "The random seed to use, DEFAULT: generate from unix time")
    flag.StringVar(&outPref, "out", "", "output file prefix to use for the randomly pseudo haploized dataset DEFAULT: DATA_PREFIX_dep")

    flag.Parse()

    args := flag.Args()

    // missing/too much arguments
    if help || (len(args) != 2) {
        printHelp()
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

    var samples []string

    // read samples (small files) so we can test whether deplete sample IDs are OK
    if fext == ".geno" {
        samples     = correctKin.ReadIND(prefix + ".ind")
    } else if fext == ".bed" {
        samples     = correctKin.ReadFAM(prefix + ".fam")
    } else {
    // not the data genom files provided?
        printHelp()
    }

    // have default out prefix if omitted
    if outPref == "" {
        outPref = prefix + "_dep"
    // have a warning that outpref is same as prefix, so file will be overwritten
    } else if outPref == prefix {
        fmt.Fprintf(os.Stderr, "WARNING: output file will overwrite the input file %s\n", args[0])
    }

    // parse depstring and die if format is invalid or missing sample, bad fractions provided
    depleteInfo := parseDepString(args[1], samples)

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
    // deplete sample groups to expected values
    for _, depData := range depleteInfo {
        depleteMarkers(depData, samples)
    }

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
