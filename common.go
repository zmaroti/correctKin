package correctKin

import (
    "os"
    "bufio"
    "fmt"
)

// struct to store data from .snp (EIG) or .bim (PLINK) file
type SNP struct {
    ID  string  // usually rs number ID
    CHR string  // chromosome, supposed to be number only but other organizm could have different naming or scaffolds
    MAP string  // genome position in centimorgan, should be float but use string to keep missing data what else
    POS int     // can be only INT
    REF string  // expected REF allele (major for plink)
    ALT string  // expected ALT allele (minor for plink)
}

// gets the linecount of a text file
func LineCount(fn string) int {
    inFile, err := os.Open(fn)
    defer inFile.Close()

    if err != nil {
        // handle error
        fmt.Println(err)
        os.Exit(1)
    }

    scanner := bufio.NewScanner(inFile)
    scanner.Split(bufio.ScanLines) 

    var lineCount int

    for scanner.Scan() {
        lineCount++
    }

    if err := scanner.Err(); err != nil {
        fmt.Fprintf(os.Stderr, "Reading %s: %v", fn, err)

        os.Exit(1)
    }

    return lineCount
}

