'''
Calculate a moving window min distance score.

What I'm doing:
1. At each position in the amino acid sequence, calculate the number of appearance of
the motif in the window. This is done for all 3 frames.

2. In the original frame, mark where the # of appearances in a window exceeds the threshold.

3. In each of frames 2 and 3, find the locations where the density is highest.
Record these highest locations, but only at locations meeting the threshold requirement found in 2.

arguments: input fasta file, output directory, output format, motifs...
'''
import csv, h5py, numpy, os, sys
from Bio import SeqIO, Motif, Seq, Alphabet
from matplotlib import pyplot

class ClusterScorer:

    # Global Settings:
    WINDOW_SIZE = 50 #Window size for density calculatioons
    DNA_ALPHA = Alphabet.IUPAC.unambiguous_dna
    AA_ALPHA = Alphabet.IUPAC.protein
    KERNEL = numpy.ones(WINDOW_SIZE) #basically counts occurences of motif falling in the window
    THRESHOLD = 2 #minimum number of appearances in the original frame to be interesting
    INPUT_FORMAT = 'fasta'

    def writeCsv(self, gene, score, positions):
        '''Write output as csv'''
        row = [gene, score] + positions
        self.csv_writer.writerow(row)

    def writeHdf5(self, gene, score, positions):
        '''Write output as hdf5'''
        self.out_handle[gene] = (score, positions)

    def writeDistribution(self):
        with open(os.path.join(self.output_dir, 'all_locs.csv'), 'w') as out_handle:
            for num in self.master_distr:
                out_handle.write('{0:03.2f}\n'.format(num))
        pyplot.hist(self.master_distr, bins=10)
        pyplot.title('Where Max Densities Occur in a Sequence')
        pyplot.xlabel('Relative Location')
        pyplot.savefig(os.path.join(self.output_dir, 'histogram.png'))

    def analyze_record(self, record):
        '''Pick locations over a threshold and with max density'''
        sequences = (Seq.Seq(str(record.seq)[i:], self.DNA_ALPHA) for i in range(3))
        aaSequences = (seq.translate() for seq in sequences)
        all_locations = []
        windowed_counts = []

        for aaSeq in aaSequences:
            # Mark motif locations with 1, else 0
            locations = numpy.zeros(len(aaSeq), dtype=int)
            for pos, seq in self.motif_inp.search_instances(aaSeq):
                locations[pos] = 1
            all_locations.append(locations)
            windowed_counts.append(numpy.convolve(locations, self.KERNEL, mode='valid'))

        # Windowing is hard-coded for 3 frames        
        last = min((len(V) for V in windowed_counts))
        windowed_counts = [V[:last] for V in windowed_counts] #trim all to same size
        select = windowed_counts[0] > self.THRESHOLD #select sites meeting threshold
        if numpy.any(select):
            # trim first part of sequence since they're not covered by the window
            all_locations[1] = all_locations[1][(self.WINDOW_SIZE - 1):]
            all_locations[2] = all_locations[2][(self.WINDOW_SIZE - 1):]

            score1 = numpy.max(windowed_counts[1][select])
            score1_sel = numpy.logical_and(windowed_counts[1] == score1, select)
            hits1 = numpy.nonzero(numpy.logical_and(all_locations[1], score1_sel))[0]

            score2 = numpy.max(windowed_counts[2][select])
            score2_sel = numpy.logical_and(windowed_counts[2] == score2, select)
            hits2 = numpy.nonzero(numpy.logical_and(all_locations[2], score2_sel))[0]

            if score1 > score2:
                hit = hits1
                score = score1
            elif score2 > score1:
                hit = hits2
                score = score2
            elif score2 > 0:
                hit = numpy.hstack((hits1, hits2))
                score = score2
            else:
                hit = numpy.array([])
                score = -1

            normalized_locations = numpy.divide(hit, last).tolist()
            self.master_distr.extend(normalized_locations)
        else:
            score = -1
            normalized_locations = []

        self.writer(record.id, score, normalized_locations)

    def __init__(self, inputPath, output_dir, motifInput, outputMode):
        '''Calculates motif densities.
        File inputPath contains FASTA DNA sequences.
        The density used will be for for all motifs listed in motifInput.
        Results are written to output_dir using outputMode = 'csv' or 'hdf5'
        '''
        self.motif_inp = Motif.Motif(alphabet=self.AA_ALPHA)
        for motif in motifInput:
            self.motif_inp.add_instance(Seq.Seq(motif, alphabet=self.AA_ALPHA))

        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if outputMode == 'csv':
            self.writer = self.writeCsv
            self.out_handle = open(os.path.join(output_dir, 'out.csv'), 'w', newline='')
            self.csv_writer = csv.writer(self.out_handle)
        elif outputMode == 'hdf5':
            self.writer = self.writeHdf5
            self.out_handle = h5py.File(os.path.join(output_dir, 'out.h5'), 'w')
        else:
            raise Exception('unknown mode ' + outputMode)

        self.master_distr = [] # Combined relative locations for all genes

        with open(inputPath, 'r') as in_handle:
            for record in SeqIO.parse(in_handle, self.INPUT_FORMAT, alphabet=self.DNA_ALPHA):
                self.analyze_record(record)

        self.out_handle.close()
        self.writeDistribution()


def main():
    inputPath = sys.argv[1]
    outputDir = sys.argv[2]
    mode = sys.argv[3]
    motifInput = sys.argv[4:]
    ClusterScorer(inputPath, outputDir, motifInput, mode)

if __name__ == '__main__':
    main()
