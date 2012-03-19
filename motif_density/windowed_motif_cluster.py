'''
Calculate a moving window min distance score.

What I'm doing:
1. At each position in the amino acid sequence, I calculate the number of appearance of
the motif in the window. This is done for all 3 frames.

2. I make a mask for when the # of appearances exceeds a certain threshold in the original frame

3. I take the max. of appearances in frames 2 and 3 in a given window, but only at locations
allowed by the mask.

input: fasta dna sequence
output: a score for each gene
'''
import h5py, numpy, sys
from Bio import SeqIO, Motif, Seq, Alphabet

class ClusterScorer:

    # Global Settings:
    WINDOW_SIZE = 50 #Window size for density calculatioons
    DNA_ALPHA = Alphabet.IUPAC.unambiguous_dna
    AA_ALPHA = Alphabet.IUPAC.protein
    KERNEL = numpy.ones(WINDOW_SIZE) #basically counts occurences of motif falling in the window
    THRESHOLD = 2 #minimum number of appearances in the original frame to be interesting
    INPUT_FORMAT = 'fasta'

    def writeCsv(self, gene, score, position):
        '''Write output as csv'''
        self.out_handle.write(bytes('{0},{1},{2}\n'.format(gene, score, position), 'UTF-8'))

    def writeHdf5(self, gene, score, position):
        '''Write output as hdf5'''
        self.out_handle[gene] = score

    def analyze_record(self, record):
        '''Does density calculations and writing of output for a given DNA sequence in record.'''
        sequences = (Seq.Seq(str(record.seq)[i:], self.DNA_ALPHA) for i in range(3))
        aaSequences = (seq.translate() for seq in sequences)
        windowed_counts = []

        for aaSeq in aaSequences:
            # Mark motif locations with 1, else 0
            locations = numpy.zeros(len(aaSeq), dtype=int)
            for pos, seq in self.motif_inp.search_instances(aaSeq):
                locations[pos] = 1
            windowed_counts.append(numpy.convolve(locations, self.KERNEL, mode='valid'))

        # Windowing is hard-coded for 3 frames
        last = min((len(V) for V in windowed_counts))
        windowed_counts = [V[:last] for V in windowed_counts] #trim all to same size
        select = windowed_counts[0] > self.THRESHOLD
        if numpy.any(select):
            score = numpy.max(numpy.maximum(windowed_counts[1][select], windowed_counts[2][select]))
            position = numpy.argmax(numpy.maximum(windowed_counts[1][select], windowed_counts[2][select]))
        else:
            score = -1
            position = -1

        self.writer(record.id, score, position / last)

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
    mode = sys.argv[3]
    motifInput = sys.argv[4:]
    ClusterScorer(inputPath, outputPath, motifInput, mode)

if __name__ == '__main__':
    main()
