#!/usr/bin/env python3

import gzip
import sys
from typing import Protocol

class VariantLike(Protocol) :
    @property
    def variant_information(self) -> dict :
        ...

class PaxVariant(VariantLike) :
    def __init__(self, items) :
        self.chrom = items['Chr']
        self.start = int(items['Start'])
        self.end   = int(items['End'])
        self.ref   = items['Ref']
        self.alt   = items['Alt']
        if self.ref == '-' :
            self.mut_type = 'Ins'
            self.mut_length = len(self.alt)
        elif self.alt == '-' :
            self.mut_type = 'Del'
            self.mut_length = len(self.ref)
        else :
            raise ValueError(f'Invalid mutation type ref: {self.ref}, alt: {self.alt}')
        self._variant_information = items

    @property
    def variant_information(self) -> dict :
        return self._variant_information

class PaxReader :
    essential_colnames = ['Chr', 'Start', 'End', 'Ref', 'Alt']
    def __init__(self, file_name, header = None) :
        if file_name.endswith('.gz') :
            self.file_handler = gzip.open(file_name, 'rt')
        else :
            self.file_handler = open(file_name, 'r')
        if header is not None :
            self.header = header
        else :
            self.header = self.file_handler.readline().rstrip('\n').split('\t')
        for colname in self.essential_colnames :
            assert(colname in self.header)
    
    def __iter__(self) :
        return self

    def __next__(self) :
        line = self.file_handler.readline()
        if line == '' :
            raise StopIteration
        cols = line.rstrip('\n').split('\t')
        assert(len(self.header) == len(cols))
        items = {self.header[i]: cols[i] for i in range(len(self.header))}
        return(PaxVariant(items))

class PaxWriter :
    def __init__(self, file_name, header) :
        self.header = header
        self.file_handler = open(file_name, 'wt')
    
    def write_head(self) :
        self.file_handler.write('\t'.join(self.header) + '\n')
    
    def write(self, pax_variant: VariantLike) :
        output_cols = [pax_variant.variant_information[k] for k in self.header]
        output_line = '\t'.join([str(x) for x in output_cols]) + '\n'
        self.file_handler.write(output_line)

if __name__ == '__main__' :
    pax_reader = PaxReader(sys.argv[1])
    pax_writer = PaxWriter(sys.argv[2], pax_reader.header)
    pax_writer.write_head()
    for l in pax_reader :
        pax_writer.write(l)
