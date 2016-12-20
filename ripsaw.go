package main

import (
	"fmt"
	"math"
	"os"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"

	cli "gopkg.in/urfave/cli.v1"
)

type dinucCount [5][5]int
type dinucDistr [14]float64

func (c *dinucCount) GC() float64 {
	var gc, at int
	for i := 0; i < len(c); i++ {
		at += c[i][0]
		gc += c[i][1]
		gc += c[i][2]
		at += c[i][3]
	}
	return (float64(gc) / (float64(gc + at)))
}

// Very basic RIP index calculation.
// Returns the ratio of (ApC+GpT) to (CpA+TpG).
func (c *dinucCount) RipIndex() float64 {
	return (float64(c[0][1]+c[2][3]) / float64(c[1][0]+c[3][2]))
}

func (c *dinucCount) Distribution() *dinucDistr {
	var s [14]int
	s[0] = c[0][0] + c[3][3]  // AA,TT
	s[1] = c[0][1] + c[2][3]  // AC,GT
	s[2] = c[0][2] + c[1][3]  // AG,CT
	s[3] = c[0][3] + c[0][3]  // AT
	s[4] = c[1][0] + c[3][2]  // CA,TG
	s[5] = c[1][1] + c[2][2]  // CC,GG
	s[6] = c[1][2] + c[1][2]  // CG
	s[7] = c[2][1] + c[2][1]  // GC
	s[8] = c[3][0] + c[3][0]  // TA
	s[9] = c[0][4] + c[4][3]  // AN,NT
	s[10] = c[1][4] + c[4][2] // CN,NG
	s[11] = c[2][4] + c[4][1] // GN,NC
	s[12] = c[3][4] + c[4][0] // TN,NA
	s[13] = c[4][4] + c[4][4] // NN

	total := 0.0
	for _, count := range s {
		total += float64(count)
	}

	var p dinucDistr
	for i, count := range s {
		p[i] = float64(count) / total
	}

	return (&p)
}

func (p *dinucDistr) ShannonEntropy() (entropy float64) {
	entropy = 0
	for _, prob := range p {
		if prob > 0 {
			entropy -= (prob * math.Log2(prob))
		}
	}
	return (entropy)
}

func main() {
	app := cli.NewApp()
	app.Name = "ripsaw"
	app.Usage = "Segment the genome according to dinucleotide frequencies"
	app.Action = func(c *cli.Context) error {
		var filename string
		if c.NArg() > 0 {
			filename = c.Args().Get(0)
		} else {
			fmt.Printf("No filename supplied.\n\n")
			cli.ShowAppHelp(c)
			os.Exit(2)
		}
		f, err := os.Open(filename)
		check(err)

		t := &linear.Seq{}
		t.Alpha = alphabet.DNA
		sc := seqio.NewScanner(fasta.NewReader(f, t))
		for sc.Next() {
			s := sc.Seq().(*linear.Seq)

			it := s.Alphabet().LetterIndex()
			var baseStart, ncounter int
			baseStart = 0
			ncounter = 0
			for i, nuc := range s.Seq {
				// Are we in a region of 'n's?
				if it[nuc] < 0 {
					// If so, increment the count.
					ncounter++
				} else {
					// We're in good sequence
					if ncounter > 5 {
						if ncounter != i {
							contig := linear.NewSeq(s.Name(), s.Seq[baseStart:i-ncounter], alphabet.DNA)
							contig.Offset = baseStart
							AnalyseContig(contig, 0)
						}
						baseStart = i
						ncounter = 0
					} else {
						ncounter = 0
					}
				}
			}

			i := s.Seq.Len()
			contig := linear.NewSeq(s.Name(), s.Seq[baseStart:i-ncounter], alphabet.DNA)
			contig.Offset = baseStart
			AnalyseContig(contig, 0)
		}
		return nil
	}
	app.Run(os.Args)
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func AnalyseContig(contig *linear.Seq, depth int) {
	it := contig.Alphabet().LetterIndex()
	encoded := make([]int, contig.Len())
	// Switch from base to integer value:
	// A => 0, C => 1, G => 2, T => 3, other => 4
	for i, nucleotide := range contig.Seq {
		if it[nucleotide] >= 0 {
			encoded[i] = it[nucleotide]
		} else {
			encoded[i] = 4
		}
	}

	var lFrac, rFrac float64
	var lCounts, rCounts dinucCount
	for i := 1; i < len(encoded); i++ {
		rCounts[encoded[i-1]][encoded[i]]++
	}

	totalEntropy := rCounts.Distribution().ShannonEntropy()

	maxEntropy := 0.0
	maxEntropyIndex := 0
	for i := 1; i < contig.Len()-1; i++ {
		lCounts[encoded[i-1]][encoded[i]]++
		rCounts[encoded[i-1]][encoded[i]]--
		lDistr := lCounts.Distribution()
		rDistr := rCounts.Distribution()
		lFrac = (float64(i) + 1) / float64(contig.Len())
		rFrac = 1.0 - lFrac
		djs := totalEntropy - (lFrac*lDistr.ShannonEntropy() + rFrac*rDistr.ShannonEntropy())
		if djs > maxEntropy {
			maxEntropy = djs
			maxEntropyIndex = i
		}
	}

	lContig := linear.NewSeq(contig.Name(), contig.Seq[0:maxEntropyIndex], alphabet.DNA)
	rContig := linear.NewSeq(contig.Name(), contig.Seq[maxEntropyIndex:], alphabet.DNA)

	if lContig.Len() > 5000 && rContig.Len() > 5000 {
		lContig.Offset = contig.Offset
		AnalyseContig(lContig, depth+1)
		rContig.Offset = contig.Offset + maxEntropyIndex
		AnalyseContig(rContig, depth+1)
	} else {
		fmt.Printf("%s\t%d\t%d\td_%02d\t%.2f\t%.2f",
			contig.Name(),
			contig.Annotation.Offset,
			contig.Annotation.Offset+contig.Len(),
			depth,
			lCounts.RipIndex()*50,
			lCounts.GC())
		for _, prob := range lCounts.Distribution() {
			fmt.Printf("\t%.4f", prob)
		}
		fmt.Println()
	}
}
