// CLI tool to
//     create the corrected kinship coefficient matrix from marker overlap and kinship coeff estimation matrices
//     statistical analysis of corrected kinship coeff matrix to infer confidence intervals
//     filter potential relatives 
 
package main

import (
    "fmt"
    "os"
    "path"
    "bufio"
    "strings"
    "strconv"
    "flag"
    "sort"
    "math"
    "github.com/kshedden/gonpy"
)


type overlapKinCoeff struct {
    overlapFrac     float64
    corrKinCoeff    float64
}

// sort by overlap
type overlapSort []overlapKinCoeff
func (slice overlapSort) Len() int {
    return len(slice)
}
func (slice overlapSort) Less(i, j int) bool {
    return slice[i].overlapFrac < slice[j].overlapFrac
}
func (slice overlapSort) Swap(i, j int) {
    slice[i], slice[j] = slice[j], slice[i]
}

// sort by kinship
type kinSort []overlapKinCoeff
func (slice kinSort) Len() int {
    return len(slice)
}
func (slice kinSort) Less(i, j int) bool {
    return slice[i].corrKinCoeff < slice[j].corrKinCoeff
}
func (slice kinSort) Swap(i, j int) {
    slice[i], slice[j] = slice[j], slice[i]
}

type overlapBin struct {
    minOverlap          float64
    medianOverlap       float64
    maxOverlap          float64
    corrKinSD           float64
    corrKinAvg          float64
    sortedOvKin         kinSort
}

// converts an array of strings (representing floats) to float type slice
func strArr2Float(strArr []string) []float64 {
    floatArr := make([]float64, len(strArr))

    var err error

    for i, str := range strArr {
        floatArr[i], err = strconv.ParseFloat(str, 64)

        if err != nil {
            panic(err)
        }
    }

    return floatArr
}


// check that CLI provided sample IDs are present in the data set
func getSampleIdxs(idsarr []string, sampleIds []string) []int {
    idxs := make([]int, len(idsarr))

    for i, id := range idsarr {
        var found bool

        for j, sampleId := range sampleIds {
            if sampleId == id {
                idxs[i] = j
                found = true

                break
            }
        }

        if !found {
            fmt.Fprintf(os.Stderr, "Sample ID %s is not in the kinship matrix file\n", id)
            os.Exit(1)
        }
    }

    // sort indexes by sample order
    sort.Ints(idxs)

    return idxs
}

// sample IDs are not provided, only the float64 numbers of the estimated kin coeffs
// the dimensions should be the same as in the overlap matrix
func readNumpyKinMat(fn string) [][]float64 {
    r, err := gonpy.NewFileReader(fn)

    if err != nil {
        panic(err)
    }
    
    if len(r.Shape) != 2 || (r.Shape[0] != r.Shape[1]) {
        panic("Invalid kinship coeff file (it should be a 2 dimensional symmetrc matrix)")
    }
   
    data, err := r.GetFloat64()

    if err != nil {
        panic(err)
    }

    kinCoeffMat := make([][]float64, r.Shape[0])

    // gonpy handles col/row major but it is symmetric anyway
    // so the flattened data is always OK
    // copy flat numpy to 2d coeff matrix
    for i := 0; i < r.Shape[0]; i++ {
        kinCoeffMat[i] = make([]float64, r.Shape[0])

        for j := 0; j < r.Shape[1]; j++ {
            kinCoeffMat[i][j] = data[i * int(r.Shape[0]) + j]
        }
    }

    return kinCoeffMat
}

// read the marker overlap pairwise matrix, this data also holds the sample IDs in
// first row and column that is also returned
func readOverlap(fn string) ([][]float64, []string) {
    inFile, err := os.Open(fn)
    defer inFile.Close()

    if err != nil {
        // handle error
        fmt.Println(err)
        os.Exit(2)
    }

    scanner := bufio.NewScanner(inFile)

    buf := make([]byte, 1024*1024)

    // we may have large matrices
    scanner.Buffer(buf, bufio.MaxScanTokenSize)

    scanner.Split(bufio.ScanLines)

    // read first line
    scanner.Scan()

    // first entry is "ID"
    sampleIds := strings.Split(scanner.Text(), "\t")

    sampleCount := len(sampleIds) - 1

    matrix := make([][]float64, sampleCount)

    var lineCount int

    for scanner.Scan() {
        arr := strings.Split(scanner.Text(), "\t")

        if len(arr) -1  != sampleCount {
            fmt.Fprintf(os.Stderr, "Invalid number of data in %s (expected %d, got %d)\n", fn, sampleCount, len(arr) - 1)
            os.Exit(1)
        }

        if lineCount > sampleCount {
            fmt.Fprintf(os.Stderr, "Extra number of lines in %s (expected %d)\n", fn, sampleCount + 1)
            os.Exit(1)
        }

        matrix[lineCount] = strArr2Float(arr[1:])

        lineCount++
    }

    if err := scanner.Err(); err != nil {
        fmt.Fprintln(os.Stderr, "reading marker overlap file:", err)

        os.Exit(1)
    }

    return matrix, sampleIds[1:]
}

// calculates the corrected kinship coefficient matrix from the marker overlap
// and estimated kinship coeff matrices in the upper triangle, lower triangle holds
// the original uncorrected kinship coeff estimates
func calcCorrKinMat(overlapMat, kinCoeffMat [][]float64) [][]float64 {
    corrKinMat := make([][]float64, len(kinCoeffMat))


    // we fill upper triangle with corrected, and leave the lower triangle with
    // the uncorrected kinship coeff estimates
    for i, estCoeffs := range kinCoeffMat {
        corrKinMat[i] = make([]float64, len(estCoeffs))

        for j, estCoeff := range estCoeffs {
            if j > i {
                corrKinMat[i][j] = estCoeff / overlapMat[i][j]
            } else {
                corrKinMat[i][j] = estCoeff
            }
        }
    }

    return corrKinMat
}

// write the matrix to a file with sample Ids
func writeCorrKinMat(fn string, sampleIds []string, corrKinMat [][]float64) {
    outFile, err := os.Create(fn)
    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        os.Exit(2)
    }
    defer outFile.Close()

    writer := bufio.NewWriter(outFile)

    // write header
    writer.WriteString(fmt.Sprintf("ID\t%s\n", strings.Join(sampleIds, "\t")))

    strDat := make([]string, len(sampleIds))

    for i, corrCoeffs := range corrKinMat {
        for j, corrCoeff := range corrCoeffs {
            strDat[j] = fmt.Sprintf("%f", corrCoeff)
        }

        writer.WriteString(fmt.Sprintf("%s\t%s\n", sampleIds[i], strings.Join(strDat, "\t")))
    }

    writer.Flush()
}

// write statistics from the marker overlap bins
func writeStats(fn string, overlapBins []overlapBin, sigma float64) {
    outFile, err := os.Create(fn)
    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        os.Exit(2)
    }
    defer outFile.Close()

    writer := bufio.NewWriter(outFile)

    // write header
    writer.WriteString(fmt.Sprintf("bin ID\tmin overlap frac\tmax overlap frac\tmedian overlap frac\tmean (corr. coeff)\tSD (corr. coeff)\t%.2f sigma (corr. coeff)\t95%% confidence (corr. coeff)\tbin size\n", sigma))

    for i, ob := range overlapBins {
        Nsigma := ob.corrKinSD * sigma
        conf95 := ob.corrKinSD * 1.96

        writer.WriteString(fmt.Sprintf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",
            i, ob.minOverlap, ob.maxOverlap, ob.medianOverlap, 
            ob.corrKinAvg, ob.corrKinSD, Nsigma, conf95, len(ob.sortedOvKin)))
    }

    writer.Flush()
}

// returns the average and SD of the kinship coefficients of
// a kinship bin data
func calcKinBin(data kinSort) (float64, float64) {
    var corrAvg, corrSum, corrSqDiffSum float64

    fCount := float64(len(data))

    for _, okc := range data {
        corrSum += okc.corrKinCoeff
    }

    corrAvg = corrSum / fCount

    for _, okc := range data {
        corrDiff := okc.corrKinCoeff - corrAvg

        corrSqDiffSum += corrDiff * corrDiff
    }

    corrSD := math.Sqrt(corrSqDiffSum / fCount)

    return corrAvg, corrSD
}

func copyOverlapBin(minIdx, maxIdx int, allOverlapKinCoeff overlapSort) overlapBin {
    var ob overlapBin

    ob.sortedOvKin = make(kinSort, 0, maxIdx - minIdx)

    for idx := minIdx; idx < maxIdx; idx++ {
        ob.sortedOvKin = append(ob.sortedOvKin, allOverlapKinCoeff[idx])
    }

    // we have the data sorted by overlap fraction so we can get the median overlap fraction
    // since we have 1000+ data, skip accounting for odd/even bins, just get close to median is ok
    // this overlap median is only for plotting anyway
    ob.medianOverlap = ob.sortedOvKin[len(ob.sortedOvKin)/2].overlapFrac
    ob.minOverlap = ob.sortedOvKin[0].overlapFrac
    ob.maxOverlap = ob.sortedOvKin[len(ob.sortedOvKin)-1].overlapFrac

    // re sort it by kinship coeff
    sort.Sort(ob.sortedOvKin)

    return ob
}

func calcOverlapBins(overlapMat, corrKinMat [][]float64) []overlapBin {
    // in the matrix we have N*(N-1)/2 coefficients
    numCoeff := len(corrKinMat) * (len(corrKinMat)-1) / 2

    allOverlapKinCoeff := make(overlapSort, 0, numCoeff)

    // get upper triangle corrected kinship coefficients
    for i, corrCoeffs := range corrKinMat {
        for j := i + 1; j < len(corrKinMat); j++ {
            allOverlapKinCoeff = append(allOverlapKinCoeff,
                overlapKinCoeff{overlapFrac: overlapMat[i][j], corrKinCoeff: corrCoeffs[j]})
        }
    }

    sort.Sort(allOverlapKinCoeff)

    var ovStep float64 = 0.01

    nextOvFrac  := allOverlapKinCoeff[0].overlapFrac + ovStep

    var minIdx int
    maxIdx := len(allOverlapKinCoeff)

    ovB := make([]overlapBin, 0, 100)

    for idx, ovk := range allOverlapKinCoeff {
        if ovk.overlapFrac >= nextOvFrac {
            if idx - minIdx > 499 {
                ovB = append(ovB, copyOverlapBin(minIdx, idx, allOverlapKinCoeff))

                minIdx = idx
                nextOvFrac= allOverlapKinCoeff[idx].overlapFrac + ovStep
            }
        } else if maxIdx - idx < 1000 {
            ovB = append(ovB, copyOverlapBin(minIdx, maxIdx, allOverlapKinCoeff))
            break
        }
    }

    // calculate overlap frac median, kin coeff SD, etc for bins
    // our previous experiments showed that the technical error from the whole workflow has
    // symetric distribution
    // so for SD calculations we check the most negative kinship coeff in the bin and exclude positive
    // kinship coeffs above the absolute value of that (that might be real kins, or bad reference pop
    // over estimated kins) so we can estimate the SD for technical errors
    for i, ob := range ovB {
        // we have kinship coeff sorted so first is the most negative
        minKinCoeff := ob.sortedOvKin[0].corrKinCoeff

        var maxIdx int = len(ob.sortedOvKin) - 1

        // we should have +/- values but be on safe
        if minKinCoeff < 0 {
            for maxIdx = len(ob.sortedOvKin) - 1; maxIdx > 0; maxIdx-- {
                if ob.sortedOvKin[maxIdx].corrKinCoeff <= -minKinCoeff {
                    break
                }
            }
        }

        ovB[i].corrKinAvg, ovB[i].corrKinSD = calcKinBin(ob.sortedOvKin[0:maxIdx])
    }

    return ovB
}

// add indexes of a pair of relatives in a slice of grouped indexes
//   (one idx group contains all the indexes of individuals who are related)
// since we traverse only upper triangle all idx pairs (idx1, idx2) are unique and not repeated in reverse order
// there are only three cases
// #1 idx1 and idx2 are not in any groups -> add a new group containing idx1, idx2
// #2 either idx1 or idx2 is in one group but the other is not, -> add other idx to this group
// #3 idx1 and idx2 are in different groups -> join the two groups so now it contains everyone including idx1 and idx2
func addIdxs(idx1, idx2 int, idxGroups *[][]int) {
    var groupIdx1 int = -1
    var groupIdx2 int = -1

    for grIdx, idxGroup := range *idxGroups {
        for _, idx := range idxGroup {
            if idx == idx1 {
                groupIdx1 = grIdx
            } else if idx == idx2 {
                groupIdx2 = grIdx
            }
        }
    }

    // idx1 and idx2 are not in any groups (yet)
    if groupIdx1 == -1 && groupIdx2 == -1 {
        *idxGroups = append(*idxGroups, []int{idx1, idx2})

    // idx1 is an elem of the groupIdx1 th group, so add idx2 to this group
    } else if groupIdx1 != -1 && groupIdx2 == -1 {
        (*idxGroups)[groupIdx1] = append((*idxGroups)[groupIdx1], idx2)

    // idx1 is an elem of the groupIdx1 th group, so add idx2 to this group
    } else if groupIdx1 == -1 && groupIdx2 != -1 {
        (*idxGroups)[groupIdx2] = append((*idxGroups)[groupIdx2], idx1)
    // if the two groups ate not the same, join two groups + delete the one joined in
    } else if groupIdx1 != groupIdx2 {
        // append elements to group1 from group2
        (*idxGroups)[groupIdx1] = append((*idxGroups)[groupIdx1], (*idxGroups)[groupIdx2]...)

        // remove group2
        (*idxGroups)[groupIdx2] = (*idxGroups)[len(*idxGroups)-1]
        *idxGroups = (*idxGroups)[:len(*idxGroups)-1]
    }
    // otherwise nothing to do
}

// get the overlap bin for error model
func findOverlapBin(overlapBins []overlapBin, overlapFrac float64) overlapBin {
    var ob overlapBin
    for _, ob = range overlapBins {
        if overlapFrac >= ob.minOverlap && overlapFrac <= ob.maxOverlap {
            return ob
        }
    }

    // shouldn't happen, exit if we get here
    fmt.Fprintf(os.Stderr, "Overlap fraction not in overlap bins: %f\n", overlapFrac)
    os.Exit(1)

    return ob
}

// returns a slice of []int containing gropup of indexes of the related individuals
// based on N sigma and kinship coeff threshold
func filterGroups(corrKinMat, overlapMat [][]float64, overlapBins []overlapBin,
    sigma, coeffThresh float64) [][]int {

    var idxGroups [][]int

    sampleCount := len(corrKinMat)
    // traverese upper triangle (corrected kinship values)
    for idx1 := 0; idx1 < sampleCount; idx1++ {
        for idx2 := idx1 + 1; idx2 < sampleCount; idx2++ {

            // filter ids above the kinship coeff threshold
            if corrKinMat[idx1][idx2] > coeffThresh {
                overlapFrac := overlapMat[idx1][idx2]

                ob := findOverlapBin(overlapBins, overlapFrac)

                corrSDthresh := ob.corrKinAvg + (sigma * ob.corrKinSD)

                // filter ids above the N sigma threshold of the overlapbin
                if corrKinMat[idx1][idx2] > corrSDthresh {
                    addIdxs(idx1, idx2, &idxGroups)
                }
            }
        }
    }

    // sort the indexes in groups (sample order)
    for i, _ := range idxGroups {
        sort.Ints(idxGroups[i])
    }

    return idxGroups
}

// filterIdxs(corrKinMat, overlapMat, overlapBins, sigma, coeffThresh, testIdxGroups)
// returns a slice of []int containing gropup of indexes of the filtered individuals
func filterIdxs(corrKinMat, overlapMat [][]float64, overlapBins []overlapBin,
    sigma, coeffThresh float64, testIdxGroups [][]int) [][]int {

    var idxGroups [][]int

    for _, testIdxs := range testIdxGroups {
        // indexes are sorted by sample order so we can traverse the upper triangle
        for i, idx1 := range testIdxs {
            for j := i + 1; j < len(testIdxs); j++ {
                idx2 := testIdxs[j]

                // we have IDs and coeffThreshold provided, just add id pairs that are above the thresh
                // and also above the n * sigma threshold of variance of the given overlap bin
                if corrKinMat[idx1][idx2] > coeffThresh {
                    overlapFrac := overlapMat[idx1][idx2]

                    ob := findOverlapBin(overlapBins, overlapFrac)

                    corrSDthresh := ob.corrKinAvg + (sigma * ob.corrKinSD)

                    if corrKinMat[idx1][idx2] > corrSDthresh {
                        addIdxs(idx1, idx2, &idxGroups)
                    }
                }
            }
        }
    }

    return idxGroups
}

// classifies the relation to a text format from kinship coefficients
// note this is an approoximation, there can be biolgocial and technical variances
// at distant kins the difference between the expected kinship estimates will be smaller
// than the confidence interval of the corrected kinship coefficient
func classifRelation(kinshipCoeff float64, kinCoeffTreshes[]float64, kinClassif []string) string {
    for i, kinCoeffTresh := range kinCoeffTreshes {
        if kinshipCoeff >= kinCoeffTresh {
            return kinClassif[i]
        }
    }

    return "uncertain"
}

type avgGroups struct {
    Avgs []float64
    Groups  [][]int
}

func (b avgGroups) Len() int           { return len(b.Avgs) }
func (b avgGroups) Less(i, j int) bool { return b.Avgs[i] > b.Avgs[j] }
func (b avgGroups) Swap(i, j int) {
    b.Groups[i], b.Groups[j] = b.Groups[j], b.Groups[i]
    b.Avgs[i], b.Avgs[j] = b.Avgs[j], b.Avgs[i]
}

// sorting groups by average corrected kin coeff
// so we get JOINT, SAMPLE Dups, close relatives first
func sortGroups(groupIdxs [][]int, corrKinMat [][]float64) [][]int {
    avgs := make([]float64, len(groupIdxs))

    // indexes in groups are sample order
    for i, group := range groupIdxs {
        var groupCorrKinSum float64
        var groupCount int

        for j, idx1 := range group {
            for k := j + 1; k < len(group); k++ {
                idx2 := group[k]

                groupCorrKinSum += corrKinMat[idx1][idx2]
                groupCount++
            }
        }

        avgs[i] = groupCorrKinSum / float64(groupCount)
    }

    sort.Sort(avgGroups{Avgs: avgs, Groups: groupIdxs})

    return groupIdxs
}

// print kin groups to STDOUT in flat format
func printFlat(sampleIds []string, groupIdxs [][]int, corrKinMat, overlapMat [][]float64,
    overlapBins []overlapBin, showSelf bool, coeffThresh, sigma float64) {

    var skip int

    if !showSelf {
        skip = 1
    }

    coeff := 0.5

    kinCoeffThreshes := make([]float64, 6)

    var kinClassif []string = []string{"DUP/MZT", "1st", "2nd", "3rd", "4th", "5th"}

    for i := 0; i <= 5; i++  {
        kinCoeffThreshes[i]   = (coeff + (coeff / 2)) / 2

        coeff = coeff / 2
    }

    // print header
    fmt.Printf("kin group\tID1\tID2\tuncorr. kin. coeff.\toverlap frac.\tcorr. kin. coeff.\t%.2f sigma threshold\t95 conf (lower)\t95 conf (upper)\test. relatedness\n", sigma)

    for i, group := range groupIdxs {
        // indexes are sorted by sample order so we can have the upper triangle data

        for j, idx1 := range group {
            for k := j + skip; k < len(group); k++ {
                idx2 := group[k]

                var classif string

                overlapFrac := overlapMat[idx1][idx2]

                ob := findOverlapBin(overlapBins, overlapFrac)

                // corrected kinship coeff SD, 95% conf interval, 95% conf int. upper/lower limit
                kinSDthresh := ob.corrKinAvg + (sigma * ob.corrKinSD)
                kinConf95   := ob.corrKinAvg + (1.96 * ob.corrKinSD)
                kinConf95Upper   := corrKinMat[idx1][idx2] + kinConf95
                kinConf95Lower   := corrKinMat[idx1][idx2] - kinConf95

                if idx1 == idx2 {
                    classif = "SELF"
                    
                // when same pseudo-haploid data is included in two samples we can have >>0.5 kinship coeff estimations
                // as ref pops and everyone else have true random haploidized calls, while
                // these samples share significant portion of same state pseudo-haplopid calls
                } else if corrKinMat[idx1][idx2] > 0.6 {
                    classif = "JOINT DATA"
                } else if  corrKinMat[idx1][idx2] > kinSDthresh {
                    classif = classifRelation(corrKinMat[idx1][idx2], kinCoeffThreshes, kinClassif)
                } else {
                    classif = "uncertain"
                }

                fmt.Printf("Group%d\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n", i + 1, sampleIds[idx1], sampleIds[idx2],
                    corrKinMat[idx2][idx1], overlapMat[idx1][idx2], corrKinMat[idx1][idx2], kinSDthresh,
                        kinConf95Lower, kinConf95Upper, classif)
            }
        }
    }
}

// print kin groups to STDOUT in matrix format
func printMatrix(sampleIds []string, groups [][]int, matrix [][]float64, showSelf bool) {
    for i, group := range groups {
        // indexes are sorted in sample oder

        strmat := make([][]string, len(group) + 1)

        for j, _ := range strmat {
            strmat[j] = make([]string, len(group) + 1)
        }
        
        for j, idx1 := range group {
            strmat[0][j+1] = sampleIds[idx1]
            strmat[j+1][0] = sampleIds[idx1]

            for k, idx2 := range group {
                if (j == k) && !showSelf {
                    strmat[j+1][k+1] = ""
                } else {
                    strmat[j+1][k+1] = fmt.Sprintf("%f", matrix[idx1][idx2])
                }
            }
        }

        fmt.Printf("Group%d\n", i + 1)

        for _, row := range strmat {
            fmt.Printf("%s\n", strings.Join(row, "\t"))
        }

        fmt.Println("")
    }
}

func printHelp() {
    fmt.Fprintln(os.Stderr,
`USAGE
filterRelates [OPTIONS] <PREFIX.overlap> <PREFIX.kinship.npy>

filterRelates will read the the pairwise marker overlap fraction of the samples generated by 'markerOverlap' and the kinship coefficient matrix (output of PCangsd v.0.99)

It then corrects the estimated kinship coefficients by the marker overlap fraction and writes out the corrected kinship coefficient file (PREFIX.corr.tsv) as a flat tab separated file (upper triangle corrected kinship coeff, lower triangle the original kinship coefficient estimate).

Calculates statistics on the whole pairwise kinship coefficient matrix where presumably the majority of relations between the individuals are unrelated. Based on this statistics and the overlap fractions it estimates the standard deviation of corrected kinship coefficient for different bins of overlap fractions and writes the statistics (PREFIX.stats.tsv) including the overlap fraction range, mean corrected kinship coeff, SD of corrected kinship coeff, N sigma threshold, 95% confidence interval.

And lastly filters close relatives based on sigma threshold and/or kinship coefficient threshold and prints in FLAT or MATRIX format to the STDOUT.

As default it prints the kinship coefficients for the group of relatives in flat format where the columns contains the following information
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
    The default is false, meaning we have a FLAT output. If this flag is set then the matrix of kinship coefficients between individuals of each group of relatives is printed (upper triangle holds the corrected, lower triangle holds the uncorrected kinship coefficients)

-showself bool|DEFAULT false
    DEFAULT we don't print kiship coefficients "between" SELF. To override this use the -showself flag. In case of matrix ouput this flag affects whether to display the (self) kinship coefficients at the DIAGONAL or not.

-ids ID1,ID2:ID3,ID4,ID5
    when ID groups are explictly provided the kinship coefficient matrix will be only tested for the ID combinations in each ID groups. In this example it would ouput kinship coefficients between ID1 vs ID2; ID3 vs ID4; ID3 vs ID5 and ID4 vs ID5
`)

    os.Exit(0)
}

func main() {
    var coeffThresh, sigma  float64
    var matOut, showSelf    bool
    var idsStr, outPref     string

    flag.StringVar(&outPref,  "out", "", "output file prefix to use. Default: fileprefix of overlap file")
    flag.BoolVar(&matOut,     "matrix",    false, "The relates are printed as the kinship coeff matrix of the related individuals")
    flag.BoolVar(&showSelf,   "showself",   false,  "Print kinship relations with self")
    flag.StringVar(&idsStr,   "ids",       "",    "Comma separated list of ID pairs to lookup in a form of ID1,ID2:ID2,ID3,ID4,...")
    flag.Float64Var(&sigma,   "sigma",    6.0,     "Statisticaly significant kinship coeff treshold (sigma in SD|DEFULAT 6)")
    flag.Float64Var(&coeffThresh, "thresh",   0.0,     "Kinship coefficient threshold to consider relatives")
 
    flag.Parse()

    args := flag.Args()

    if len(args) != 2 {
        printHelp()
    }

    // read overlap file
    overlapMat, sampleIds := readOverlap(args[0])

    // read estimated kinship coeff matrix from npy
    kinCoeffMat := readNumpyKinMat(args[1])

    // calculate the corrected kinship coeff from the two matrices
    corrKinMat := calcCorrKinMat(overlapMat, kinCoeffMat)

    fext := path.Ext(args[0])

    // get prefix from overlap file prefix
    if outPref == "" {
        outPref = args[0][0:len(args[0]) - len(fext)]
    }

    // write the corrected kinship coeff matrix into prefix + ".corr.tsv"
    writeCorrKinMat(outPref + ".corr.tsv", sampleIds, corrKinMat)

    overlapBins := calcOverlapBins(overlapMat, corrKinMat)

    // write out statistics about the corrected kin coeffs in the overlap bins
    writeStats(outPref + ".stats.tsv", overlapBins, sigma)

    // group of indexes of sample IDs that are kins
    var groupIdxs [][]int

    // explicitly given IDs to consider
    if idsStr != "" {
        groups := strings.Split(idsStr, ":")

        testIdxs := make([][]int, len(groups))

        for i, group := range groups {
            idsArr := strings.Split(group, ",")

            testIdxs[i] = getSampleIdxs(idsArr, sampleIds)
        }

        // we have threshold criteria and Ids provided, so only consider Ids from the list where
        // the corrected kinship coeff is above the threshold
        if coeffThresh > 0  {
            groupIdxs = filterIdxs(corrKinMat, overlapMat, overlapBins, sigma, coeffThresh, testIdxs)
        } else {
        // include all test idx combos regardless of kinship coeffs
            groupIdxs = testIdxs
        }
    // no explicitly given IDs, filter all kin groups based on thresh and/or bins/sigma
    } else {
        groupIdxs = filterGroups(corrKinMat, overlapMat, overlapBins, sigma, coeffThresh)
    }

    groupIdxs = sortGroups(groupIdxs, corrKinMat)

    // print the kin groups
    // matrix output
    if matOut {
        printMatrix(sampleIds, groupIdxs, corrKinMat, showSelf)
    } else {
    // flat output with statistics as well
        printFlat(sampleIds, groupIdxs, corrKinMat, overlapMat, overlapBins, showSelf, coeffThresh, sigma)
    }
}
