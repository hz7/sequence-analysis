#!/usr/bin/env python3

'''
This script, given motifs and an input FASTA sequence file,
finds the positions of the motifs in the 3 frames.

Usage: locator INPUT_FILE OUTPUT_FILE OUTPUT_MODE MOTIFS...

e.x. locator orf_coding.fasta results.csv csv SP TP

or if importing from another module/ calling from the interpreter:
Locator(orf_coding.fasta, results.csv, ['SP, 'TP'], 'csv')
Locator(orf_coding.fasta, results.h5, ['SP, 'TP'], 'hdf5')
'''

import h5py, numpy, sys
from Bio import SeqIO, Motif, Seq, Alphabet

class Locator:

    # Global Settings:
    DNA_ALPHA = Alphabet.IUPAC.unambiguous_dna
    AA_ALPHA = Alphabet.IUPAC.protein
    INPUT_FORMAT = 'fasta'

    def writeCsv(self, gene, all_locs):
        '''Write output as csv'''
        self.out_handle.write(gene)
        self.out_handle.write('\n')
        for loc in all_locs:
            self.out_handle.write(','.join(str(num) for num in loc))
            self.out_handle.write('\n')

    def writeHdf5(self, gene, all_locs):
        '''Write output as hdf5'''
        i = 0
        group = self.out_handle.create_group('/{0}'.format(gene))
        for locs in all_locs:
            group[str(i)] = numpy.array(locs, dtype=numpy.int32)
            i += 1

    def analyze_record(self, record):
        '''Does density calculations and writing of output for a given DNA sequence in record.'''
        sequences = [Seq.Seq(str(record.seq)[i:], self.DNA_ALPHA) for i in range(3)]
        aaSequences = [seq.translate() for seq in sequences]

        all_locs = []
        for aaSeq in aaSequences:
            # Encoded as length, position1, position2, ...
            locations = [len(aaSeq)]
            for pos, seq in self.motif_inp.search_instances(aaSeq):
                locations.append(pos)
            all_locs.append(locations)

        self.writer(record.id, all_locs)

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
            self.out_handle = open(outputPath, 'w')
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
    mode = sys.argv[3]
    motifInput = sys.argv[4:]
    Locator(inputPath, outputPath, motifInput, mode)

if __name__ == '__main__':
    main()
