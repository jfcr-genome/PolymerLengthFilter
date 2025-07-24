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
    
    def _calc_mutect2_vaf(self, sample_info) -> dict :
        allelic_depth = [int(x) for x in sample_info['AD'].split(',')]
        alt_depth = sum(allelic_depth[1:])
        ref_depth = allelic_depth[0]
        total_depth = int(sample_info['DP'])
        vaf = float(alt_depth) / float(total_depth)
        return {'ref_reads':   ref_depth, 
                'alt_reads':   alt_depth,
                'total_reads': total_depth,
                'vaf':         vaf}
    
    def _calc_strelka2_vaf(self, sample_info) -> dict :
        tir_list = [int(x) for x in sample_info['TIR'].split(',')]
        tar_list = [int(x) for x in sample_info['TAR'].split(',')]
        alt_depth = tir_list[0]
        ref_depth = tar_list[0]
        total_depth = ref_depth + alt_depth
        vaf = float(alt_depth) / float(total_depth)
        return {'ref_reads':   ref_depth, 
                'alt_reads':   alt_depth,
                'total_reads': total_depth,
                'vaf':         vaf}

    def calc_vaf(self) :
        if self.variant_information['mutect2.FORMAT'] != 'NA' :
            format_header =    self.variant_information['mutect2.FORMAT'].split(':')
            normal_info_list = self.variant_information['mutect2.NORMAL'].split(':')
            tumor_info_list  = self.variant_information['mutect2.TUMOR'].split(':')
            assert(len(format_header) == len(normal_info_list))
            assert(len(format_header) == len(tumor_info_list))
            normal_info = {format_header[i]: normal_info_list[i] for i in range(len(format_header))}
            tumor_info  = {format_header[i]: tumor_info_list[i]  for i in range(len(format_header))}
            normal_reads_info = self._calc_mutect2_vaf(normal_info)
            for k, v in normal_reads_info.items() :
                self._variant_information['normal_' + k] = v
            tumor_reads_info  = self._calc_mutect2_vaf(tumor_info)
            for k, v in tumor_reads_info.items() :
                self._variant_information['tumor_' + k] = v
        elif self.variant_information['strelka2.FORMAT'] != 'NA' :
            format_header =    self.variant_information['strelka2.FORMAT'].split(':')
            normal_info_list = self.variant_information['strelka2.NORMAL'].split(':')
            tumor_info_list  = self.variant_information['strelka2.TUMOR'].split(':')
            assert(len(format_header) == len(normal_info_list))
            assert(len(format_header) == len(tumor_info_list))
            normal_info = {format_header[i]: normal_info_list[i] for i in range(len(format_header))}
            tumor_info  = {format_header[i]: tumor_info_list[i]  for i in range(len(format_header))}
            normal_reads_info = self._calc_strelka2_vaf(normal_info)
            for k, v in normal_reads_info.items() :
                self._variant_information['normal_' + k] = v
            tumor_reads_info  = self._calc_strelka2_vaf(tumor_info)
            for k, v in tumor_reads_info.items() :
                self._variant_information['tumor_' + k] = v
        else :
            sys.exit("No Mutect2 or Strelka2 Information")
        return

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
