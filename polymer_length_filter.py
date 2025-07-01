#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO
import argparse
from pax_io import PaxReader, PaxWriter, PaxVariant, VariantLike
from context_annotator import ContextAnnotator, AnnotatedVariant
from reference_reader import ReferenceReader

parser = argparse.ArgumentParser()
parser.add_argument('pax_result_file')
parser.add_argument('output_file')
parser.add_argument('fasta_file')
parser.add_argument('-v', '--vaf_threshold', default = 0.25)
parser.add_argument('-r', '--repeat_length', default = 6)
parser.add_argument('-l', '--mutation_length', default = 2)
parser.add_argument('-s', '--compute_spectrum', action = 'store_true')
args = parser.parse_args()


pax_reader = PaxReader(args.pax_result_file)
if 'tumor_vaf' not in pax_reader.header :
    sys.exit("There is no vaf column.")
ref_reader = ReferenceReader(args.fasta_file)
output_header = pax_reader.header + \
    ['annotation_label', 'rep_seq', 'rep_length', 'mh_seq', 'mh_length']
pax_writer = PaxWriter(args.output_file, output_header)
pax_writer.write_head()
cont_annot = ContextAnnotator(ref_reader)
for variant in pax_reader :
    if 'str_contraction' in variant.variant_information['mutect2.FILTER'].split(';') :
        continue
    annotated_variant = cont_annot.annotate(variant)
    # print(annotated_variant.label + '\t' + annotated_variant.variant_information['tumor_vaf'] + '\t' + str(annotated_variant.mut_len))
    if (annotated_variant.annot_type in ['C', 'T', 'R'] and 
        float(annotated_variant.variant_information['tumor_vaf']) < args.vaf_threshold and 
        annotated_variant.mut_len <= args.mutation_length and 
        annotated_variant.annot_len >= args.repeat_length)  :
        continue
    pax_writer.write(annotated_variant)

