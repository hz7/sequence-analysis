#!/usr/bin/env python3

'''
This script, given motifs and an input FASTA sequence file,
finds the positions of the motifs in the 3 frames, then gets min. distances.

Usage: distance INPUT_FILE OUTPUT_FILE OUTPUT_MODE MOTIFS...

e.x. distance orf_coding.fasta results.csv csv SP TP

or if importing from another module/ calling from the interpreter:
DistFinder(orf_coding.fasta, results.csv, ['SP, 'TP'], 'csv')
DistFinder(orf_coding.fasta, results.h5, ['SP, 'TP'], 'hdf5')
'''

import h5py, numpy, sys
from Bio import SeqIO, Motif, Seq, Alphabet

def min_dist(v1, v2):
    ''' 
    Assume v1, v2 are sorted.
    Returns a list containing the minimum distance for ea. element in v2 to
    some element in v1.
    Because v1 and v2 are sorted, this function is linear because we know
    if x < y in v2 and r < s in v1, if d(x, s) is minimal, r < x <= s
    and d(r, y) < d(s, y)
    '''
    result = []
    if (len(v1) == 0):
        return result
    i = 0 #index of current candidate for closest element
    for x in v2:
        d1 = float("inf") #current distance
        d2 = abs(v1[i] - x) #next distance
        while d2 <= d1:
            d1 = d2
            if i < len(v1) - 1:
                i += 1
                d2 = abs(v1[i] - x)
            else:
                break
        i -= 1 #backtrack
        result.append(d1)
    return result

def get_dists(all_locs):
    diff1 = (x2 - x1 for x1,x2 in zip(all_locs[0], all_locs[0][1:]))
    diff2 = min_dist(all_locs[0], all_locs[1])
    diff3 = min_dist(all_locs[0], all_locs[2])
    return [diff1, diff2, diff3]

class DistFinder:

    # Global Settings:
    DNA_ALPHA = Alphabet.IUPAC.unambiguous_dna
    AA_ALPHA = Alphabet.IUPAC.protein
    INPUT_FORMAT = 'fasta'

    def writeCsv(self, gene, all_dists):
        '''Write output as csv'''
        self.out_handle.write(gene)
        self.out_handle.write('\n')
        for dists in all_dists:
            self.out_handle.write(','.join(str(num) for num in dists))
            self.out_handle.write('\n')

    def writeHdf5(self, gene, all_dists):
        '''Write output as hdf5'''
        i = 0
        group = self.out_handle.create_group('/{0}'.format(gene))
        for dists in all_dists:
            group[str(i)] = numpy.array(dists, dtype=numpy.int32)
            i += 1

    def analyze_record(self, record):
        '''Does density calculations and writing of output for a given DNA sequence in record.'''
        sequences = [Seq.Seq(str(record.seq)[i:], self.DNA_ALPHA) for i in range(3)]
        aaSequences = [seq.translate() for seq in sequences]

        all_locs = []
        for aaSeq in aaSequences:
            # Encoded as length, position1, position2, ...
            locations = []
            for pos, seq in self.motif_inp.search_instances(aaSeq):
                locations.append(pos)
            all_locs.append(locations)

        self.writer(record.id, get_dists(all_locs)) #record distance calculation

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
    DistFinder(inputPath, outputPath, motifInput, mode)

if __name__ == '__main__':
    main()
