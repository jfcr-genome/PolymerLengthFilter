#!/usr/bin/env python3

import os
import os.path
import sys
from Bio import SeqIO
import argparse
from pax_io import PaxReader, PaxWriter, PaxVariant, VariantLike
from context_annotator import ContextAnnotator, AnnotatedVariant
from reference_reader import ReferenceReader
import spectrum_context

GENOME_BUILD = 'GRCh38'

parser = argparse.ArgumentParser()
parser.add_argument('pax_result_file')
parser.add_argument('output_file')
parser.add_argument('fasta_file')
parser.add_argument('-v', '--vaf_threshold', default = 0.25)
parser.add_argument('-r', '--repeat_length', default = 6)
parser.add_argument('-l', '--mutation_length', default = 2)
parser.add_argument('--compute_spectrum', action = 'store_true')
parser.add_argument('-s', '--sample', default = "Sample")
args = parser.parse_args()

calculate_vaf_flag = False

context_83_count = {}
if args.compute_spectrum :
    for annotation_label in spectrum_context.context_83 :
        context_83_count[annotation_label] = 0

pax_reader = PaxReader(args.pax_result_file)
output_header = pax_reader.header + \
    ['annotation_label', 'rep_seq', 'rep_length', 'mh_seq', 'mh_length']
if 'tumor_vaf' not in pax_reader.header :
    calculate_vaf_flag = True
    output_header += ['normal_total_reads', 'normal_ref_reads', 'normal_alt_reads', 'normal_vaf',
                      'tumor_total_reads',  'tumor_ref_reads',  'tumor_alt_reads',  'tumor_vaf']
    # sys.exit("There is no vaf column.")
ref_reader = ReferenceReader(args.fasta_file)
pax_writer = PaxWriter(args.output_file, output_header)
pax_writer.write_head()
cont_annot = ContextAnnotator(ref_reader)
for variant in pax_reader :
    if 'str_contraction' in variant.variant_information['mutect2.FILTER'].split(';') :
        continue
    if calculate_vaf_flag :
        variant.calc_vaf()
    annotated_variant = cont_annot.annotate(variant)
    # print(annotated_variant.label + '\t' + annotated_variant.variant_information['tumor_vaf'] + '\t' + str(annotated_variant.mut_len))
    if (annotated_variant.annot_type in ['C', 'T', 'R'] and 
        float(annotated_variant.variant_information['tumor_vaf']) < args.vaf_threshold and 
        annotated_variant.mut_len <= args.mutation_length and 
        annotated_variant.annot_len >= args.repeat_length)  :
        continue
    pax_writer.write(annotated_variant)
    if args.compute_spectrum :
        context_83_count[annotated_variant.label] += 1

if args.compute_spectrum :
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
    import SigProfilerMatrixGenerator.scripts.MutationMatrixGenerator
    from SigProfilerAssignment import Analyzer as Analyze
    import pandas as pd
    context_count_file = os.path.splitext(args.output_file)[0] + ".context_83_count.txt"
    with open(context_count_file, 'w') as f :
        f.write('MutationType' + '\t' + args.sample + '\n')
        for annotation_label in spectrum_context.context_83 :
            f.write(annotation_label + '\t' + str(context_83_count[annotation_label]) + '\n')

    output_dir = os.path.abspath(os.path.dirname(args.output_file))

    Analyze.cosmic_fit(context_count_file, output_dir,
                       collapse_to_SBS96 = False, genome_build=GENOME_BUILD)

    signature_dat = \
        pd.read_csv(output_dir + \
                    "/Assignment_Solution/Activities/Assignment_Solution_Activities.txt", \
                    sep="\t")

    sig_cols = signature_dat.filter(like = 'ID').columns
    signature_dat['total'] = signature_dat[sig_cols].sum(axis=1)
    signature_dat[sig_cols] = signature_dat[sig_cols].div(signature_dat['total'], axis=0)

    signature_dat.to_csv(output_dir + \
                         "/Assignment_Solution/Activities/Assignment_Solution_Activities_Relative.txt",
                         sep="\t", index = False)

