// CLI tool to calculate the pairwise marker overlap fraction between samples
// from a binary PLINK, PACKEDANCESTRYMAP or a plain text EIGENSTRAT genome data.
// We defined marker ovelap fraction between two samples as the markers that are
// genotyped in both samples divided by the total number of markers in the data set.
// NOTE, only variants > 0.05 MAF (by default or the provided -maf TRESHOLD) markers
// are considered in the calculation. This is the default MAF threshold applied in
// "PCangsd -kinship" calculation.
package main

import (
    "fmt"
    "bufio"
    "os"
    "strings"
    "sync"
    "runtime"
    "flag"
    "path"
    "github.com/zmaroti/correctKin"
)

// matrix that holds typing data of samples for the markers in EIGENSTRAT format
// 0 -> HOM_ALT (0 REF ALLELE)
// 1 -> HET     (1 REF ALLELE)
// 2 -> HOM_REF (2 REF ALLELES)
// 3 -> MISSING DATA
var GENOTYPES [][]uint8

var sampleCount int

// non overlapping ranges of marker chunks (for paralelization)
type typeRange struct {
    start   int    // start of range
    end     int    // end of range
    matrix [][]int // result matrix
}

// calculate the marker overlap of samples in a subset of markers
// store the result in a matrix of counts
func calcOverlap(r <-chan typeRange, done chan<- typeRange) {
    for rx := range r {
        var data []uint8

        for n := rx.start; n <= rx.end; n++ {
            // int slice of sample typing for given marker
            data = GENOTYPES[n]

            // check that a sample pair is both typed (in the upper triangle)
            for i, c1 := range data {
                // if first sample is not typed, skip as both markers never will be typed
                if c1 == 3 {
                    continue
                }

                // otherwise test in the upper triangle whether other sample is typed
                for j := i; j < sampleCount; j++ {
                    if data[j] != 3 {
                        rx.matrix[i][j]++
                    }
                }
            }
        }

        done<-rx
    }
}

func excludeLowMAF(mafTresh float64, sampleCount, markerCount int) int {
    var snpMaf float64 // MAF of the actual SNP

    for i, rowData := range GENOTYPES {
        var typedAlleles, missingSamples int

        // our data is in PACKED ANCESTRYMAP coding
        // 3 - missing
        // 0 - HOM_ALT (0 REF ALLELE)
        // 1 - HET     (1 REF ALLELE)
        // 2 - HOM_REF (2 REF ALLELES)
        for _, gt := range rowData {
            if gt == 0 {
                typedAlleles += 2
            } else if gt == 1 {
                typedAlleles += 1
            } else if gt == 3 {
                missingSamples++
            }
        }

        typedSamples := sampleCount - missingSamples
        alleleCount  := typedSamples * 2
        
        // the minor allele is the REF allele
        if typedAlleles > typedSamples {
            snpMaf = float64(alleleCount - typedAlleles) / float64(alleleCount)
        } else {
            snpMaf = float64(typedAlleles) / float64(alleleCount)
        }

        // set all data in the row for this SNP to unknown, so it will not add to marker overlap
        // decrease total markerCount accordingly
        if snpMaf < mafTresh {
            for j, _ := range rowData {
                GENOTYPES[i][j] = 3
            }

            markerCount--
        }
    }
    return markerCount
}

// generate worker number of ranges to cover all markers
func generateRanges(workers int) <-chan typeRange {
    r := make(chan typeRange)

    max  := len(GENOTYPES)
    step := float64(max) / float64(workers)
    var start, end int
    
    tasks := make([]typeRange, workers)

    for i := 0; i < workers; i++ {
        end = int((float64(i) + 1) * step) - 1

        tasks[i].start = start
        tasks[i].end   = end

        tasks[i].matrix = make([][]int, sampleCount)

        for j := range tasks[i].matrix {
            tasks[i].matrix[j] = make([]int, sampleCount)
        }

        start = end + 1
    }

    tasks[workers-1].end = max - 1

    go func() {
        defer close(r)

        for i := 0; i < workers; i++ {
            r<-tasks[i]
        }
    }()

    return r
}

// parallelize overlap counting
func doWork(workers int) [][]int {
    input := generateRanges(workers)

    var wg sync.WaitGroup

    wg.Add(workers)

    result := make(chan typeRange)

    for i := 0; i < workers; i++ {
        go func() {
            calcOverlap(input, result)

            wg.Done()
        }()
    }

    go func() {
        wg.Wait()
        close(result)
    }()

    // initialize the final overlap count matrix
    final := make([][]int, sampleCount)
    for i := range final {
        final[i] = make([]int, sampleCount)
    }

    // add the counts of overlap chunks to the final matrix
    for res := range result {
        for i:= 0; i < sampleCount; i++ {
            for j:= i; j < sampleCount; j++ {
                final[i][j] += res.matrix[i][j]
            }
        }
    }

    return final
}

// calculate the fraction of overlap for one row from the final count matrix
// and return it as a string slice for printing
func calcFracs(tot int, vals []int) []string {
    fracs := make([]string, len(vals))

    totalCount := float64(tot)

    for i, val := range vals {
        fracs[i] = fmt.Sprintf("%f", float64(val) / totalCount)
    }

    return fracs
}

func printHelp() {
    fmt.Fprintln(os.Stderr,
`USAGE
markerOverlap [-threads n|default ALL] <DATA.(bed|geno)>

The tool calculates the pairwise genotyping overlap fractions between all samples of a genome dataset. As input it requires either a binary PLINK (.bed, .bim, .bam), an EIGENSTRAT or a PACKEDANCESTRYMAP (.geno, .ind, .snp) dataset. The format is guessed by the provided filetype (either the .bed or .geno).

The output is a plain text file (DATA.overlap) containing the pairwise marker overlap fraction matrix. Only markers above the minor allele frequency threshold (default 0.05) are counted. This maf threshold is the same as the default maf threshold in PCangsd -kinship calculation, where low AF markers pruned prior to PCA analysis. In case kinship estimation is done with a different threshold, use the same threshold for calculating marker overlap fractions. The marker overlap fraction is defined between two samples as: (count of markers genotyped in both samples)/(count of all markers in the dataset) /considering only the markers above the maf threshold for both counts/.

optional flags:
-threads value    the number opf threads to use |default ALL
-maf     value    minor allele frequency threshold (markers below the theshold are not counted) |default 0.05 
`)

    os.Exit(0)
}

func main() {
    var help    bool
    var maf     float64
    var threads int

    flag.BoolVar(&help,   "help", false, "print help")
    flag.Float64Var(&maf, "maf", 0.05, "MAF threshold, markers below MAF are not counted in the overlap, DEFULT: 0.05 (default in PCangsd)")
    flag.IntVar(&threads, "threads", 0, "Number of threads to use, DEFAULT=ALL available")

    flag.Parse()

    args := flag.Args()

    if help || (len(args) != 1) {
        printHelp()
    }
    
    fext := path.Ext(args[0])

    if fext != ".geno" && fext != ".bed" {
        printHelp()
    }

    prefix := args[0][0:len(args[0]) - len(fext)]

    var samples []string
    var markerCount int

    if fext == ".geno" {
        samples     = correctKin.ReadIND(prefix + ".ind")
        markerCount = correctKin.LineCount(prefix + ".snp")

        sampleCount = len(samples)

        _, GENOTYPES  = correctKin.ReadEIG(args[0], sampleCount, markerCount)
    } else if fext == ".bed" {
        samples     = correctKin.ReadFAM(prefix + ".fam")
        markerCount = correctKin.LineCount(prefix + ".bim")

        sampleCount = len(samples)

        _, GENOTYPES  = correctKin.ReadBED(args[0], sampleCount, markerCount)
    } else {
        printHelp()
    }

    // if we have to exclude markers by MAF, recalculate the markerCount that are above this MAF threshold
    if maf > 0 {
        markerCount = excludeLowMAF(maf, sampleCount, markerCount)
    }

    // concurrently calculate the overlap matrix between samples
    if threads == 0 {
        threads = runtime.NumCPU()
    }

    final := doWork(threads)

    // create overlap file
    overlapFn := prefix + ".overlap"

    overlapFile, err := os.Create(overlapFn)
    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        os.Exit(2)
    }
    defer overlapFile.Close()

    writer := bufio.NewWriter(overlapFile)

    writer.WriteString(fmt.Sprintf("ID\t%s\n", strings.Join(samples, "\t")))

    for i, s := range(final) {
        fracs := calcFracs(markerCount, s)

        writer.WriteString(fmt.Sprintf("%s\t%s\n", samples[i], strings.Join(fracs, "\t")))
    }

    writer.Flush()
}
