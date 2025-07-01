#!/usr/bin/env python3

import sys
from Bio import SeqIO

class ReferenceReader :
    genome_seq = {}
    def __init__(self, fasta_file) :
        for record in SeqIO.parse(fasta_file, 'fasta') :
            self.genome_seq[record.name] = record.seq
    
    def get_sequence(self, chrom, start, end=None, is_one_based_coord = True) :
        assert(chrom in self.genome_seq)
        if end == None :
            end = start + 1
        if is_one_based_coord :
            # SeqRecord is zero-based so position must be decremented
            start -= 1
            end -= 1
        
        if start < 0 or len(self.genome_seq[chrom]) <= end :
            return None
        return self.genome_seq[chrom][start:end]

if __name__ == '__main__' :
    import pax_io
    pax_result_file = sys.argv[1]
    fasta_file = sys.argv[2]
    pax_reader = pax_io.PaxReader(pax_result_file)
    reference_reader = ReferenceReader(fasta_file)
    for variant in pax_reader :
        reference_sequence = reference_reader.get_sequence(variant.chrom, variant.start)
        print(variant.ref, '\t', reference_sequence)
        assert(variant.ref == reference_sequence)
    print("OK")
