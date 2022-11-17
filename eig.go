package correctKin

import (
    "os"
    "io"
    "bufio"
    "strings"
    "strconv"
    "fmt"
)

// read an EIGENSTRAT/PACKEDANCESTRYMAP .ind file and return the sampleIDs in a string slice
func ReadIND(fn string) []string {
    inFile, err := os.Open(fn)
    defer inFile.Close()

    if err != nil {
        // handle error
        fmt.Println(err)
        os.Exit(1)
    }

    scanner := bufio.NewScanner(inFile)
    scanner.Split(bufio.ScanLines) 

    var line string
    var lineno int
    var samples []string

    f := func(c rune) bool {
        return c == ' '
    }

    for scanner.Scan() {
        line = scanner.Text()
        lineno++

        fields := strings.FieldsFunc(line, f)

        if len(fields) != 3 {
            fmt.Fprintf(os.Stderr, "Invalid IND line %s at %d\n", line, lineno)

            os.Exit(1)
        }

        // .ind fields:
        // SAMPLEID SEX POPID
        samples = append(samples, fields[0])
    }

    if err := scanner.Err(); err != nil {
        fmt.Fprintln(os.Stderr, "Reading Stdin: ", err)

        os.Exit(1)
    }

    return samples
}

// read an EIEGNSTRAT/PACKEDANCESTRYMAP .snp file and return the marker data in a slice
// that can be used to create a corresponding .bim file (PLINK)
func ReadSNP(fn string, flipMinorMajor bool) []SNP {
    inFile, err := os.Open(fn)
    defer inFile.Close()

    if err != nil {
        // handle error
        fmt.Println(err)
        os.Exit(1)
    }

    scanner := bufio.NewScanner(inFile)
    scanner.Split(bufio.ScanLines) 

    var line string
    var lineno, pos int

    var SNPs []SNP

    f := func(c rune) bool {
        return c == ' '
    }

    for scanner.Scan() {
        line = scanner.Text()
        lineno++

        fields := strings.FieldsFunc(line, f)

        if len(fields) != 6 {
            fmt.Fprintf(os.Stderr, "Invalid IND line %s at %d\n", line, lineno)

            os.Exit(1)
        }

        pos, err = strconv.Atoi(fields[3])

        if err != nil {
            fmt.Fprintln(os.Stderr, err)

            os.Exit(1)
        }

        // EIGENSTRAT .snp file format
        //            ID            CHR      MAP (cm)          POS    REF ALT
        //            rs3094315     1        0.020130          752566 G A

        // note convertf will convert PLINK .bim or .pedsnp files without flipping column 5-6
        // while in EIGENSTRAT documentation column 5 stands for REF, and 6 stands for the minor ALT allele
        // in the PLINK format column 5 stands for allale 1 (usually minor allele) and column 6 stands for allele 2 (usually major allele)
        // This is not an issue still we deal with one data set and not various datasets merged that has differently flipped minor and major alleles
        // however this option is here so you can fix this issue in case you have standard PLINK format or
        // PLINK format converted by conmvertf from EIGENSTRAT
        if flipMinorMajor {
            SNPs = append(SNPs, SNP{ID: fields[0], CHR: fields[1], MAP: fields[2], POS: pos, REF: fields[5], ALT: fields[4]})
        } else {
            SNPs = append(SNPs, SNP{ID: fields[0], CHR: fields[1], MAP: fields[2], POS: pos, REF: fields[4], ALT: fields[5]})
        }
    }

    if err := scanner.Err(); err != nil {
        fmt.Fprintln(os.Stderr, "Reading Stdin: ", err)

        os.Exit(1)
    }

    return SNPs
}

func max(x, y int) int {
    if x < y {
        return y
    }
    return x
}

// internal function to read a plain text EIGENSTRAT .geno file
func readEIG(inFile io.Reader, sampleCount, markerCount int) [][]uint8 {
    typing := make([][]uint8, 0, markerCount)

    scanner := bufio.NewScanner(inFile)
    scanner.Split(bufio.ScanLines) 

    var line []byte
    var lineCount int

    // one line per marker, sample number entries in each line
    for scanner.Scan() {
        line = scanner.Bytes()

        if len(line) != sampleCount {
            fmt.Printf("Invalid line: %s\n", string(line))
            os.Exit(1)
        }

        data := make([]uint8, sampleCount, sampleCount)

        // 9 - missing
        // 0 - HOM_ALT (0 REF ALLELE)
        // 1 - HET     (1 REF ALLELE)
        // 2 - HOM_REF (2 REF ALLELES)
        for i, b := range line {
            switch b {
                case '9':
                    data[i] = 3
                // skip esoteric character to uint conversion no clues if everywhere code tab is linear 0-1-2
                case '0':
                    data[i] = 0
                case '1':
                    data[i] = 1
                case '2':
                    data[i] = 2
            }
        }

        lineCount++

        typing = append(typing, data)
    }

    if err := scanner.Err(); err != nil {
        fmt.Fprintln(os.Stderr, "Reading Stdin: ", err)

        os.Exit(1)
    }

    if lineCount != markerCount {
        fmt.Printf("Invalid marker count in EIGENSTRAT geno line: %d vs %d\n", lineCount, markerCount)
        os.Exit(1)
    }

    return typing
}

// internal function to read a PACKEDANCESTRYMAP .geno file
func readPackedEIG(reader *bufio.Reader, byteCount, sampleCount, markerCount int) [][]uint8 {
    var typing [][]uint8

    buf := make([]byte, byteCount, byteCount)

    var b byte
    var s uint

    // file size is correct, so we don't have to test whether lineCount/MarkerCount is equal
    for V := 0; V < markerCount; V++ {
        _, err := io.ReadFull(reader, buf)

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(2)
        }

        data := make([]uint8, sampleCount, sampleCount)

        for S := 0; S < sampleCount; S++ {
            if (S % 4) == 0 {
                b = buf[S/4]
                s = 0
            }

            // packed from upper to lower 
            // 11 = missing
            // 00 = hom ALT (0 copies of hom REF)
            // 01 = het     (1 copies of hom REF)
            // 10 = hom REF (2 copies of hom REF)
            data[S] = (b >> (6 - s * 2)) & 0x3

            s++
        }

        typing = append(typing, data)
    }

    return typing
}

// wrapper for reading a .geno file that is either in plain text EIGSTRAT or a
// binary PACKEDANCESTRYMAP data file and return
// the header (that is nil in case of EGIENSTRAT data file or the header of the binary PACKEDANCESTRYMAP file)
// and 2d uint8 slice containing the EIG style genotype data in flat format
func ReadEIG(fn string, sampleCount, markerCount int)  ([]byte, [][]uint8) {
    inFile, err := os.Open(fn)
    if err != nil {
        // handle error
        fmt.Fprintln(os.Stderr, err)
        os.Exit(2)
    }
    defer inFile.Close()

    fileInfo, err := inFile.Stat()
    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        os.Exit(1)
    }

    fileSize := int(fileInfo.Size())

    byteCount := sampleCount / 4

    if (sampleCount % 4) > 0 {
        byteCount++
    }

    headerLen  := max(byteCount, 48);

    // check if file size conforms with PACKEDANCESTRYMAP expected size
    if fileSize == headerLen + byteCount * markerCount {
        reader := bufio.NewReader(inFile)
        header := make([]byte, headerLen, headerLen)

        _, err := io.ReadFull(reader, header)

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(1)
        }

        var indHash, snpHash, indCount, snpCount int
        _, err = fmt.Sscanf(string(header), "GENO %d %d %x %x", &indCount, &snpCount, &indHash, &snpHash)

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(1)
        }

        if indCount != sampleCount {
            fmt.Fprintf(os.Stderr, "Sample count in geno file (%d) is not the same as in .ind file (%d)\n",
                indCount, sampleCount)
            os.Exit(1)
        }
        if snpCount != markerCount {
            fmt.Fprintf(os.Stderr, "Marker count in geno file (%d) is not the same as in .snp file (%d)\n",
                snpCount, markerCount)
        }

        return header, readPackedEIG(reader, byteCount, sampleCount, markerCount)
    // otherwise read it as a plain text EIGENSTRAT data file
    } else {
        return nil, readEIG(inFile, sampleCount, markerCount)
    }
}

// write genotype data (.geno) in EIGENSTRAT plaint text format
func WriteEIG(outFn string, TYPING [][]uint8) {
    outFile, err := os.Create(outFn)
    if err != nil {
        // handle error
        fmt.Fprintln(os.Stderr, err)
        os.Exit(1)
    }
    defer outFile.Close()

    writer := bufio.NewWriter(outFile)

    // sample data is in rows
    sampleCount := len(TYPING[0])

    // buffer to hold the plain text EIG data for one marker (without new line)
    buf := make([]byte, sampleCount, sampleCount)

    // no header, just plain TEXT, however we have to remap 3 to 9 for MISSING
    // and recode uint to byte
    for _, markerDat := range TYPING {
        for j, geno := range markerDat {
            switch geno {
                case 3:
                    buf[j] = '9'
                case 0:
                    buf[j] = '0'
                case 1:
                    buf[j] = '1'
                case 2:
                    buf[j] = '2'
            }
        }

        // hopefully this is OS agnostic for new lines
        _, err = writer.WriteString(fmt.Sprintf("%s\n", buf))

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(1)
        }
    }

    writer.Flush()
}

// write genotype data (.geno) in binary PACKEDANCESTRYMAP format
func WritePackedEIG(outFn string, header []byte, TYPING [][]uint8) {
    outFile, err := os.Create(outFn)
    if err != nil {
        // handle error
        fmt.Fprintln(os.Stderr, err)
        os.Exit(1)
    }
    defer outFile.Close()

    writer := bufio.NewWriter(outFile)

    _, err = writer.Write(header)

    if err != nil {
        // handle error
        fmt.Fprintln(os.Stderr, err)
        os.Exit(1)
    }

    // sample data is in rows
    sampleCount := len(TYPING[0])

    // number of bytes to represent a marker for each samples in binary PED
    byteCount := sampleCount / 4

    if (sampleCount % 4) > 0 {
        byteCount++
    }

    // buffer to hold the BED data of one marker for each samples
    buf := make([]byte, byteCount, byteCount)

    var b byte

    for _, dat := range TYPING {
        for j, geno := range dat {
            // we have EIG binary format bitwise, so no need to recode

            // however PACKEDANCESTRYMAP sampes are high bits->low bits order
            b |= (geno << (6 - (j % 4) * 2))

            // pack byte
            if (j % 4) == 3 {
                buf[j/4] = b
                b = 0
            }
        }

        // last byte might not exactly divideable by 4 so we must set the lower bits we have so far
        if (sampleCount % 4) > 0 {
            buf[len(buf)-1] = b
            b = 0
        }

        // write out one marker
        _, err = writer.Write(buf)

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(1)
        }
    }

    // flush writer
    writer.Flush()
}
