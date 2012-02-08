#!/usr/bin/env python3

'''
This script, given a motif and an input FASTA sequence file,
calculates the motif densities for all 3 reading frames in each gene of the sequence file.

Usage: analysis SEQUENCE INPUT_FILE OUTPUT_FILE

e.x. analysis STP orf_coding.fasta results.csv
'''

import numpy, sys
from Bio import SeqIO, Motif, Seq, Alphabet

# Global Settings:
OUT_HANDLE = None
WINDOW_SIZE = 25
DNA_ALPHA = Alphabet.IUPAC.unambiguous_dna
AA_ALPHA = Alphabet.IUPAC.protein
KERNEL = numpy.ones(WINDOW_SIZE) / WINDOW_SIZE

# Motif config:
MOTIF_INP = Motif.Motif(alphabet=AA_ALPHA)

# Analyzes a record. Then writes the gene (record id) and densities to OUT_HANDLE
def analyze_record(record):
    sequences = [record.seq, Seq.Seq(str(record.seq)[1:], DNA_ALPHA), Seq.Seq(str(record.seq)[2:], DNA_ALPHA)]
    aaSequences = [seq.translate() for seq in sequences]
    densities = []

    for aaSeq in aaSequences:
        aaMotifPos = []

        # Find motif locations
        for pos, seq in MOTIF_INP.search_instances(aaSeq):
            aaMotifPos.append(pos)

        # calculate moving average using convolution
        signal = numpy.zeros(len(aaSeq))
        signal[aaMotifPos] = 1.0
        densities.append(numpy.convolve(signal, KERNEL, mode='valid'))

    OUT_HANDLE.write(bytes(record.id + '\n', 'UTF-8'))
    for array in densities:
        numpy.savetxt(OUT_HANDLE, numpy.reshape(array, (1, len(array))),
                            fmt='%5.5f', delimiter=',')

def main():
    motifInput = sys.argv[1]
    inputPath = sys.argv[2]
    outputPath = sys.argv[3]

    MOTIF_INP.add_instance(Seq.Seq(motifInput, alphabet=AA_ALPHA))

    global OUT_HANDLE
    OUT_HANDLE = open(outputPath, 'wb')

    with open(inputPath, 'r') as in_handle:
        for record in SeqIO.parse(in_handle, 'fasta', alphabet=DNA_ALPHA):
            analyze_record(record)
    
    OUT_HANDLE.close()

if __name__ == '__main__':
    main()
