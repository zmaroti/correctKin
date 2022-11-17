package correctKin

import (
    "os"
    "io"
    "bufio"
    "strings"
    "strconv"
    "fmt"
)

func ReadBIM(fn string, flipMajorMinor bool) []SNP {
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
            fmt.Fprintf(os.Stderr, "Invalid BIM line %s at %d\n", line, lineno)

            os.Exit(1)
        }

        pos, err = strconv.Atoi(fields[3])

        if err != nil {
            fmt.Fprintln(os.Stderr, err)

            os.Exit(1)
        }
        // PLINK format
        // CHR     ID              MAP(cm) POS     MINOR   MAJOR
        // 1       rs3094315       0.02013 752566  G       A
        // NOTE
        // in PLINK documentation field 5 is allele 1 (usually minor) and field 6 is allele 2 (usually major)
        // HOWEVER convertf will not flip alleles when converting EIGENSTRAT to BED or PED, the snp file REF (field 5)
        // is kept at field 5 of the resulting .bim or .pedsnp file, while ALT (field 6) is kept at field 6 of the PLINK format
        // so we have the flipMajorMinor option to use depending on whether your PLINK file is coming from a PLINK data set
        // or converted from EIGENSTRAT data set using convertf of EIGENSTRAT package
        // We have to note that it is not an issue in any population genetic analysis when one single data set is used, and
        // that data set is not created from many datasets where the minor/major SNPs are flipped differently
        if flipMajorMinor {
            SNPs = append(SNPs, SNP{ID: fields[1], CHR: fields[0], MAP: fields[2], POS: pos, REF: fields[5], ALT: fields[4]})
        } else {
            SNPs = append(SNPs, SNP{ID: fields[1], CHR: fields[0], MAP: fields[2], POS: pos, REF: fields[4], ALT: fields[5]})
        }
    }

    if err := scanner.Err(); err != nil {
        fmt.Fprintln(os.Stderr, "Reading Stdin: ", err)

        os.Exit(1)
    }

    return SNPs    
}

// read a PLINK FAM file and return the sampleIDs as a string slice
func ReadFAM(fn string) []string {
    inFile, err := os.Open(fn)
    defer inFile.Close()

    if err != nil {
        // handle error
        fmt.Println(err)
        os.Exit(2)
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

        if len(fields) != 6 {
            fmt.Fprintf(os.Stderr, "Invalid FAM line %s at %d\n", line, lineno)

            os.Exit(1)
        }

        // FAM fields:
        // FAMILYID SAMPLEID FATHERID MOTHERID SEX DISEASESTATUS
        samples = append(samples, fields[1])
    }

    if err := scanner.Err(); err != nil {
        fmt.Fprintln(os.Stderr, "Reading Stdin: ", err)

        os.Exit(1)
    }

    return samples
}

// read a BINARY PLINK .bed data file and return its 
// header (first 3 bytes) and the unpacked data in a 2d unit8 slice
func ReadBED(fn string, sampleCount, markerCount int) ([]byte, [][]uint8) {
    var genotypes [][]uint8

    byteCount := sampleCount / 4

    if (sampleCount % 4) > 0 {
        byteCount++
    }

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
        os.Exit(2)
    }

    fileSize := int(fileInfo.Size())

    // header is 3 bytes + the data for each indiv/marker
    expectedSize := 3 + (byteCount * markerCount)

    if fileSize != expectedSize {
        fmt.Fprintf(os.Stderr, "invalid number of entries: corrupted BED? (expected %d vs %d)\n", expectedSize, fileSize)
        os.Exit(2)
    }

    reader := bufio.NewReader(inFile)

    header := make([]byte, 3, 3)

    _, err = io.ReadFull(reader, header)

    if err != nil {
        fmt.Fprintln(os.Stderr, err)
        os.Exit(2)
    }

    // The first three bytes should be 0x6c, 0x1b, and 0x01 in that order. 
    if header[0] != 0x6c || header[1] != 0x1b {
        fmt.Fprintf(os.Stderr, "not a valid bed file: %s\n", fn)
        os.Exit(2)
    }

    // Third byte is 00000001 (SNP-major) or 00000000 (individual-major)
    // should be 00000001 for plink1.9
    if header[2] != 1 {
        fmt.Fprintf(os.Stderr, "%s is in individual-major format, please recode to SNP-major format\n", fn)
        os.Exit(2)
    }
    
    buf := make([]byte, byteCount, byteCount)

    var b byte

    // we tested the size so we have
    // marker count of blocks
    // each block consist of bytecount of bytes (data of all samples for 1 marker)
    for V := 0; V < markerCount; V++ {
        _, err := io.ReadFull(reader, buf)

        if err != nil {
            fmt.Fprintln(os.Stderr, err)
            os.Exit(2)
        }

        markerData := make([]uint8, sampleCount, sampleCount)

        // https://www.cog-genomics.org/plink2/formats#bed
        // sample order in the bytes are first lower bits, and last the upper bits
        // sample1 0b000000XX
        // sample2 0b0000XX00
        // sample3 0b00XX0000
        // sample4 0bXX000000
        for S := 0; S < sampleCount; S++ {
            if (S % 4) == 0 {
                b = buf[S/4]
            }

            // recoding PLINK to our internal uniform "EIG" format
            // GT         PLINK  EIG
            // missing    0b01  0b11
            // HET        0b10  0b01
            // hom REF    0b11  0b00
            // hom ALT    0b00  0b10
            switch b & 0x3 {
                case 1:
                    markerData[S] = 3
                case 2:
                    markerData[S] = 1
                case 3:
                    markerData[S] = 0
                case 0:
                    markerData[S] = 2
            }

            b >>= 2
        }

        genotypes = append(genotypes, markerData)
    }

    return header, genotypes
}

// write PLINK binary format
func WriteBED(outFn string, header []byte, typeData [][]uint8) {
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
    sampleCount := len(typeData[0])

    // number of bytes to represent a marker for each samples in binary PED
    byteCount := sampleCount / 4

    if (sampleCount % 4) > 0 {
        byteCount++
    }

    // buffer to hold the BED data of one marker for each samples
    buf := make([]byte, byteCount, byteCount)

    var b byte

    for _, dat := range typeData {
        for j, geno := range dat {
            // recoding from internal EIG (0b00 HOMREF, 0b01 HET, 0b10 HOMALT, 0b11 missing)
            // to PLINK                   (0b11 HOMREF, 0b10 HET, 0b00 HOMALT, 0b01 missing)
            switch geno {
                case 0: // HOMREF
                    geno = 3
                case 1: // HET
                    geno = 2
                case 2: // HOMALT
                    geno = 0
                case 3: // missing
                    geno = 1
            }

            // PLINK has low bits -> high bits order
            b |= (geno << ((j % 4) * 2))

            // pack byte
            if (j % 4) == 3 {
                buf[j/4] = b
                b = 0
            }
        }

        // last byte might not exactly divideable by 4 so we must set the last byte with its lower bits
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
