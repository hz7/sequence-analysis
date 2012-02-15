#!/usr/bin/env python3

'''
This script, given a motif and an input FASTA sequence file,
calculates the motif densities for all 3 reading frames in each gene of the sequence file.

Usage: analysis INPUT_FILE OUTPUT_FILE MOTIFS...

e.x. analysis orf_coding.fasta results.csv SP TP

or if importing from another module/ calling from the interpreter:
MotifAnalyzer(orf_coding.fasta, results.csv, ['SP, 'TP'], 'csv')
MotifAnalyzer(orf_coding.fasta, results.h5, ['SP, 'TP'], 'hdf5')
'''

import h5py, numpy, sys
from Bio import SeqIO, Motif, Seq, Alphabet

class MotifAnalyzer:

    # Global Settings:
    WINDOW_SIZE = 25
    DNA_ALPHA = Alphabet.IUPAC.unambiguous_dna
    AA_ALPHA = Alphabet.IUPAC.protein
    KERNEL = numpy.ones(WINDOW_SIZE) / WINDOW_SIZE
    INPUT_FORMAT = 'fasta'

    def writeCsv(self, gene, densities):
        '''Write output as csv'''
        self.out_handle.write(bytes(gene + '\n', 'UTF-8'))
        for array in densities:
            numpy.savetxt(self.out_handle, numpy.reshape(array, (1, len(array))),
                                fmt='%5.5f', delimiter=',')

    def writeHdf5(self, gene, densities):
        '''Write output as hdf5'''
        i = 0
        group = self.out_handle.create_group('/{0}'.format(gene))
        for array in densities:
            group[str(i)] = numpy.reshape(array, (1, len(array)))
            i += 1

    def analyze_record(self, record):
        '''Does density calculations and writing of output for a given DNA sequence in record.'''
        sequences = [Seq.Seq(str(record.seq)[i:], self.DNA_ALPHA) for i in range(3)]
        aaSequences = [seq.translate() for seq in sequences]
        densities = []

        for aaSeq in aaSequences:
            # Mark motif locations with 1, else 0
            signal = numpy.zeros(len(aaSeq))
            for pos, seq in self.motif_inp.search_instances(aaSeq):
                signal[pos] = 1.0

            # calculate moving average using convolution
            densities.append(numpy.convolve(signal, self.KERNEL, mode='valid'))

        self.writer(record.id, densities)

    def __init__(self, inputPath, outputPath, motifInput, outputMode):
        '''Calculates motif densities.
        File inputPath contains FASTA DNA sequences.
        The density used will be for for all motifs listed in motifInput.
        Results are written to outputPath using outputMode = 'csv' or 'hdf5'
        '''
        self.motif_inp = Motif.Motif(alphabet=self.AA_ALPHA)
        for motif in motifInput:
            self.motif_inp.add_instance(Seq.Seq(motif, alphabet=self.AA_ALPHA))  

        if outputMode == 'csv':
            self.writer = self.writeCsv
            self.out_handle = open(outputPath, 'wb')
        elif outputMode == 'hdf5':
            self.writer = self.writeHdf5
            self.out_handle = h5py.File(outputPath, 'w')
        else:
            raise Exception('unknown mode ' + outputMode)

        with open(inputPath, 'r') as in_handle:
            for record in SeqIO.parse(in_handle, self.INPUT_FORMAT, alphabet=self.DNA_ALPHA):
                self.analyze_record(record)

        self.out_handle.close()

def main():
    inputPath = sys.argv[1]
    outputPath = sys.argv[2]
    motifInput = sys.argv[3:]
    MotifAnalyzer(inputPath, outputPath, motifInput, 'csv')

if __name__ == '__main__':
    main()
