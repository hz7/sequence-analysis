'''
Goes through the alignment folders

Outputs csv: gene, species, score
'''

import os, numpy, csv, re
from Bio import SeqIO, Motif, Seq, Alphabet

INPUT_FORMAT = 'fasta'
ALPHA = Alphabet.IUPAC.protein
WINDOW_SIZE = 50
KERNEL = numpy.ones(WINDOW_SIZE) / WINDOW_SIZE

motif = 'Q'
motif_inp = Motif.Motif(alphabet=ALPHA)
motif_inp.add_instance(Seq.Seq(motif, alphabet=ALPHA))

#only records with id matching this pattern get saved
SPEC_PATTERN = re.compile('\[(.*?)\]') 

def analyze_dir(gene, file_path):
    '''Does density calculations and writing of output for a given AA sequence in record.'''
    with open(file_path, 'r') as in_handle:
        for record in SeqIO.parse(in_handle, INPUT_FORMAT, alphabet=ALPHA):
            match = SPEC_PATTERN.search(record.id)
            if match:
                species = match.group(1)
            else:
                continue
            
            # process alignment
            aaSeq = Seq.Seq(str(record.seq).replace('-', ''), ALPHA)
            # Mark motif locations with 1, else 0
            signal = numpy.zeros(len(aaSeq))
            for pos, seq in motif_inp.search_instances(aaSeq):
                signal[pos] = 1.0
            # calculate moving average using convolution
            density = numpy.convolve(signal, KERNEL, mode='valid')

            score = density.max()
            csv_writer.writerow([gene, species, score])

csv_writer = csv.writer(open('c:/holt/out.csv', 'w', newline=''))
csv_writer.writerow(['gene', 'species', 'score']) # write headers

for dirpath, dirnames, filenames in os.walk('c:/holt/alignments'):
    for subdirname in dirnames:
        gene = subdirname
        subdir_path = os.path.join(dirpath, subdirname)
        subdir_files = os.listdir(subdir_path)
        subdir_file_path = os.path.join(subdir_path, subdir_files[0])
        analyze_dir(gene, subdir_file_path)

del(csv_writer)
