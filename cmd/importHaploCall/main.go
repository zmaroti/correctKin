// CLI tool to import ANGSD random allele call output file into binary PLINK
// or plain text EIEGNSTRAT database format.
package main

import (
    "bufio"
    "fmt"
    "os"
    "strings"
    "strconv"
    "path"
    "compress/gzip"
    "github.com/zmaroti/correctKin"
    "flag"
)

// matrix that holds genotype data of samples for the markers
var GENOTYPES [][]uint8

// this struct holds the expected REF and ALT alleles of the provided EIGENSTRAT.snp file
// required to code the data consistently
var SNPs   []correctKin.SNP

// ANGSD haploid call gzip file has a format like this
// chr     pos     major   ind0    ind1    ind2    ind3...
// 4       95344   G       G       N       G       N   ...
var header []string

type SNPfail struct {
    CHR   string
    POS   int
    REF   string
    ALT   string
    GT    string
}

var fail   [][]SNPfail

// default false: use the broken EIGENSTRAT .bim notation (flipped minor/major with flipped binary GTs)
// if set to true, then we don't flip minor/major and conform with the proper PLINK notation, so the plink data
// can be exported to VCF (--keep-allele-order --real-ref-alleles) and having major exported as REF, minor exported as ALT, and GTs coded accordingly 
var noflip bool

func parseGTs(gts []string, idx int) {
    SNP := SNPs[idx]

    // marker was already processed, no duplicate marker can
    // be in the same genome data set, so it may happen only because
    // more haplo files were provided with overlapping markers
    // panic
    if GENOTYPES[idx] != nil {
        fmt.Fprintf(os.Stderr, "Fatal error: SNP %s:%d %s/%s is duplicated in the haplo files.\n",
            SNP.CHR, SNP.POS, SNP.REF, SNP.ALT)

        os.Exit(1)
    }

    GENOTYPES[idx] = make([]uint8, len(gts))

    // the GTs are coded as one lette nucleotides, N standing for no data

    // internally we code data in EIG
    // HOM ALT - 0 (= 0 REF allele)
    // HOM REF - 2 (= 2 REF allele)
    // HET     - 1 (= 1 REF allele)
    // MISSING - 3
    for i, gt := range gts {
        if gt == "N" {
            GENOTYPES[idx][i] = 3
        } else if gt == SNP.REF {
            GENOTYPES[idx][i] = 2
        } else if gt == SNP.ALT {
            GENOTYPES[idx][i] = 0
        } else {
        // in ANGSD haploid random allele calling we don't have HET
        // so it is either a PMD or a triallelic call, either way we make this as missing for the given sample
            GENOTYPES[idx][i] = 3

        // save the position for a detailed log on failed positions per individual
            fail[i] = append(fail[i], SNPfail{CHR: SNP.CHR, POS: SNP.POS, REF: SNP.REF, ALT: SNP.ALT, GT: gt})
        }
    }
}

func readHaplo(fn string, SNPMap map[string]map[int]int) {
    inFile, err := os.Open(fn)
    defer inFile.Close()

    if err != nil {
        // handle error
        fmt.Println(err)
        os.Exit(1)
    }

    var haploCall *bufio.Scanner

    if len (fn) > 3 && fn[len(fn)-3:] == ".gz" {
        gz, err := gzip.NewReader(inFile)

        if err != nil {
            fmt.Println(err)
            os.Exit(1)
        }
        defer gz.Close()
        haploCall = bufio.NewScanner(gz)
    } else {
        haploCall = bufio.NewScanner(inFile)
    }

    // for large sample collections
    buf := make([]byte, 1024*1024)
    haploCall.Buffer(buf, bufio.MaxScanTokenSize)

    haploCall.Split(bufio.ScanLines) 

    // read first line with header
    haploCall.Scan()

    arr := strings.Split(haploCall.Text(), "\t")
    fieldCount := len(arr)

    if header != nil {
        if fieldCount != len(header) {
            fmt.Fprintln(os.Stderr, "FATAL: mixed number of individuals in provided ANGSD haploid call files.")
            os.Exit(1)
        }
    } else {
        header = arr

        // init fail data for samples
        fail = make([][]SNPfail, fieldCount - 3)
    }

    var count, misscount int

    for haploCall.Scan() {
        // NOTE ANGSD haplocall puts an extra TAB at the end of GT lines
        arr = strings.SplitN(haploCall.Text(), "\t", fieldCount)

        if fieldCount != len(arr) {
            fmt.Fprintf(os.Stderr, "FATAL: corrupted file: %s, at line %s\n", count + misscount + 2)
            os.Exit(1)
        }

        // fix last element with extra TAB
        arr[fieldCount-1] = strings.TrimRight(arr[fieldCount-1], "\t")

        pos, err := strconv.Atoi(arr[1])

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(1)
        }

        // if we have this in the position list
        if idx, ok := SNPMap[arr[0]][pos]; ok {
            parseGTs(arr[3:], idx)
            count++
        } else {
            fmt.Fprintf(os.Stderr, "WARNING: position %s:%d was not in the provided eigenstrat.snp. Did you provide the right -sites file to ANGSD?\n",
                arr[0], pos)

            misscount++
        }
    }

    if err = haploCall.Err(); err != nil {
        fmt.Fprintln(os.Stderr, "Reading file:", err)
        os.Exit(1)
    }

    fmt.Printf("%s processed %d markers found from %d entries.\n", fn, count, count + misscount)
}

func writeLog(outFn string, failData []SNPfail) {
    outFile, err := os.Create(outFn)
    if err != nil {
        // handle error
        fmt.Fprintln(os.Stderr, err)
        os.Exit(1)
    }
    defer outFile.Close()

    writer := bufio.NewWriter(outFile)

    // print header
    _, err = writer.WriteString("CHR\tPOS\tREF\tALT\tGT\n")

    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        os.Exit(1)
    }


    for _, snp := range failData {
        _, err = writer.WriteString(fmt.Sprintf("%s\t%d\t%s\t%s\t%s\n", snp.CHR, snp.POS, snp.REF, snp.ALT, snp.GT))

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(1)
        }
    }

    writer.Flush()
}

func writeLogs() {
    for i, failData := range fail {
        if failData == nil {
            continue
        }
        writeLog(header[i + 3] + ".log", failData)
    }
}

func writeBIM(outFn string) {
    outFile, err := os.Create(outFn)
    if err != nil {
        // handle error
        fmt.Fprintln(os.Stderr, err)
        os.Exit(1)
    }
    defer outFile.Close()

    writer := bufio.NewWriter(outFile)

    // .bim file is white space separated
    for _, snp := range SNPs {
        if noflip {
            _, err = writer.WriteString(fmt.Sprintf("%s %s %s %d %s %s\n", 
                snp.CHR, snp.ID, snp.MAP, snp.POS, snp.ALT, snp.REF))
        } else {
            _, err = writer.WriteString(fmt.Sprintf("%s %s %s %d %s %s\n", 
                snp.CHR, snp.ID, snp.MAP, snp.POS, snp.REF, snp.ALT))
        }

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(1)
        }
    }

    writer.Flush()
}

func writeSNP(outFn string) {
    outFile, err := os.Create(outFn)
    if err != nil {
        // handle error
        fmt.Fprintln(os.Stderr, err)
        os.Exit(1)
    }
    defer outFile.Close()

    writer := bufio.NewWriter(outFile)

    // .snp file is space aligned, the order of fields is
    //     SNPId              CHR       MAP             POS     REF ALT
    //     snp_22_51161070    22        0.740455        51161070 G A
    for _, snp := range SNPs {
        _, err = writer.WriteString(fmt.Sprintf("%20s %6s %16s %16d %s %s\n", 
            snp.ID, snp.CHR, snp.MAP, snp.POS, snp.REF, snp.ALT))

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(1)
        }
    }

    writer.Flush()
}

func printHelp() {
    fmt.Fprintln(os.Stderr,
`USAGE
importHaploCall [-noflip default|false] <EIGENSTRAT.snp> <outfile.(bed|geno)> <list of ANGSD.haplo.gz>

The tool will read ANGSD random pseudo haploid call genotype data and writes the corresponding PACKEDANCESTRYMAP (EIGENSTRAT '.geno' and '.snp') or binary PLINK ('.bed' and '.bim') files. As input it requires an EIGENSTRAT.snp file that contains the required marker informations (CHR, POS, REF, ALT) to generate the genotype data. In case more haploid call files are provided the software will take it that it is a scatter gather genotype of the same individuals for different genome regions. In case the genotype of a given individual does not match the expected REF or ALT alleles the chromosome, position expected REF/ALT and the GT of the individual will be written in a log file (indX.log). NOTE, the '.ind' or '.fam' files has to be created manually as ANGSD names individuals in your BAM list not from the BAM header but always as 'ind0', 'ind1', ...

-noflip
     By default importHaplocall mimicks the convertf behaviour (flipped GT/minor/major) and creates a PLINK data set that can be imported to EIGENSTRAT by convertf to result in a proper GRCh37 and AADR data set conformant '.snp' file. However, this also means that conforming with the "broken" behaviour of convertf, our PLINK data set cannot be exported to a proper VCF file with the true REF/ALT alleles. To allow this option, we added the '-noflip' option to the importHaploCall tool. Using this option the imported PLINK data set will conform with the PLINK proper major/minor order. Hence data imported with the '-nopflip' option can be exported to VCF with the proper GRCh37 REF/ALT alleles and can be compared to VCF data, or merged with VCF data.
`)

    os.Exit(0)
}


func main() {
    flag.BoolVar(&noflip,  "noflip",  false, "Flip GT and minor/major in .bim file. Default false (we conform with proper PLINK notation). Use this option to import in the 'flipped' convertf/EIGENSTRAT style.")
    flag.Parse()

    Args := flag.Args()

    if len(Args) < 3 {
        printHelp()
    }

    outFn := Args[1]
    fext := path.Ext(outFn)

    // not the data genom files provided?
    if fext != ".geno" && fext != ".bed" {
        printHelp()
    }

    if fext == ".geno" && noflip == true {
        fmt.Println("WARNING: -noflip option only affects the PLINK data set format.")
    }

    // figure out prefix from provided filename
    prefix := outFn[0:len(outFn) - len(fext)]

    fmt.Println("Parsing SNP coordinates")

    SNPs  = correctKin.ReadSNP(Args[0])

    GENOTYPES = make([][]uint8, len(SNPs))

    // create a map, where keys are chromsome, genome position and the value is the index of SNP from the
    // eigenstrat.snp file
    SNPMap := make(map[string]map[int]int)

    for i, SNP := range SNPs {
        // create submap at first occurance of chromosome
        if _, ok := SNPMap[SNP.CHR]; !ok {
            SNPMap[SNP.CHR] = make(map[int]int)
        }

        SNPMap[SNP.CHR][SNP.POS] = i
    }
    
    fmt.Println("Reading ANGSD haploid call files")

    for i := 2; i < len(Args); i++ {
        readHaplo(Args[i], SNPMap)
    }

    // fix missing SNPs to missing data
    // markers where the GT slice is nil -> slice with GTs missing (3)
    for i, GTARR := range GENOTYPES {
        if GTARR == nil {
            GTARR = make([]uint8, len(header)-3)

            for j, _ := range GTARR {
                GTARR[j] = 3
            }

            GENOTYPES[i] = GTARR
        // if noflip is true and output is bed, flip the HOMREF/HOMALT values as we also flip major/minor in the .bim file in this case
        // this way PINK data set will be coded as in PLINK specification, minor in 5th column, major in 6th column, and GTs accordingly
        } else if noflip == true && fext == ".bed" {
            for j, gt := range GTARR {
                if gt == 2 {
                    GENOTYPES[i][j] = 0
                } else if gt == 0 {
                    GENOTYPES[i][j] = 2
                }
            }
        }
    }

    if fext == ".geno" {
        correctKin.WriteEIG(outFn, GENOTYPES)

        writeSNP(prefix + ".snp")
    } else {
        // The first three bytes of PLINK 1.9 should be 0x6c, 0x1b, and 0x01 (this last is SNP-major)
        correctKin.WriteBED(outFn, []byte{0x6c, 0x1b, 0x01}, GENOTYPES)

        writeBIM(prefix + ".bim")
    }

    writeLogs()
}
