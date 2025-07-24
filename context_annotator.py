#!/usr/bin/env python3

import reference_reader
from pax_io import VariantLike, PaxVariant

class AnnotatedVariant(VariantLike): 
    """
    Represents a genomic variant with sequence context annotation

    This class stores the repeat length and microhomology length calculated from the reference 
    genome sequence aroud the variant, as well as the variant class assigned based on these features.

    Attributes:
        chrom (str) : Chromosome name
        start (int) : Start position of the variant (1-based)
        end (int)   : End position of the variant
        ref (str)   : Reference allele
        alt (str)   : Alternative allele
        variant (PaxVariant) : Variant information loaded from Pax result file.
        annotation (dict) : Repeat length and microhomology information annotated by ContextAnnotator
        label (str) : Variant class assigned based on repeat, microhomology, and variant length.

    
    """
    def __init__(self, variant: PaxVariant, annotation: dict) :
        self.chrom = variant.chrom
        self.start = variant.start
        self.end   = variant.end
        self.ref   = variant.ref
        self.alt   = variant.alt
        self.variant = variant
        self.annotation = annotation
        self.label = self._summarize_annotation()
        self.annotation['annotation_label'] = self.label

    def _summarize_annotation(self) -> str :
        """
        Returns a label representing the indel mutation type and its contextual classification.
        
        The label consists of four fields separated by colons:
            <mutation_length>:<mutation_type>:<annotatype_type>:<annotation_length>

        The format follows conventions inspired by SigProfileMatrixGenerator.
        """
        if self.variant.mut_type == 'Del' :
            mut_seq = self.variant.ref
        elif self.variant.mut_type == 'Ins' :
            mut_seq = self.variant.alt
        else :
            raise ValueError(f'Invalid mutation type: {self.variant.mut_type}')
        
        mut_len = len(mut_seq)
        if mut_len == 1 :
            if mut_seq in ['C', 'G'] :
                annotation_type = 'C'
            elif mut_seq in ['T', 'A'] :  
                annotation_type = 'T'
            else :
                raise ValueError('Invalid mutation nucleotide: ' + self.ref + ' -> ' + self.alt)
            annotation_length = self.annotation['rep_length']
        else :
            if self.annotation['mh_length'] == 0 or self.annotation['rep_length'] > 0 :
                annotation_type = 'R'
                annotation_length = int(self.annotation['rep_length'] / mut_len)
            else :
                annotation_type = 'M'
                annotation_length = self.annotation['mh_length']

        # For insertions labeld as Microhomology ('M'), 
        # reclassify as Regular ('R') with repeat length zero.
        # This mirrors behavior in SigProfilerMatrixGenerator (line 1788), 
        # though the rationale is not clearly documented. 
        if self.variant.mut_type == 'Ins' and annotation_type == 'M' :
            annotation_type = 'R'
            annotation_length = 0

        self.mut_len = mut_len
        self.annot_len = annotation_length
        self.annot_type = annotation_type
        annotation_length = min(annotation_length, 5)
        mut_len = min(mut_len, 5)
        return ':'.join([str(mut_len), self.variant.mut_type, annotation_type, str(annotation_length)])

    @property
    def variant_information(self) -> dict:
        return {**self.variant.variant_information, **self.annotation}


class ContextAnnotator :
    def __init__(self, reference_reader:reference_reader.ReferenceReader) :
        self.reference_reader = reference_reader

    def annotate(self, variant) -> AnnotatedVariant:
        rep_annot = self._annotate_repeat_length(variant)
        mh_annot  = self._annotate_microhomology(variant)
        annot_dict = {**rep_annot, **mh_annot}
        return AnnotatedVariant(variant, annot_dict)
    
    def _annotate_repeat_length(self, variant) :
        genome_length = len(self.reference_reader.genome_seq[variant.chrom])
        if variant.mut_type == 'Ins' :
            mut_seq = variant.alt
            mut_length = len(variant.alt)
            forward_rep_unit_start = variant.start + 1
            reverse_rep_unit_start = variant.start + 1 - mut_length
        elif variant.mut_type == 'Del' :
            mut_seq = variant.ref
            mut_length = len(variant.ref)
            forward_rep_unit_start = variant.start + mut_length
            reverse_rep_unit_start = variant.start - mut_length
        else :
            raise ValueError(f'Invalid mutation type: {variant.mut_type}')

        forward_rep = ''
        while (forward_rep_unit_start + mut_length < genome_length and
               mut_seq == 
               self.reference_reader.get_sequence(variant.chrom, 
                                                  forward_rep_unit_start, 
                                                  forward_rep_unit_start + mut_length)) :
            forward_rep += mut_seq
            forward_rep_unit_start += mut_length

        reverse_rep = ''
        while (0 < reverse_rep_unit_start and
               mut_seq == 
               self.reference_reader.get_sequence(variant.chrom, 
                                                  reverse_rep_unit_start, 
                                                  reverse_rep_unit_start + mut_length)) :
            reverse_rep += mut_seq
            reverse_rep_unit_start -= mut_length
        
        return {'forward_rep': forward_rep, 'reverse_rep': reverse_rep, 
                'rep_seq': forward_rep + reverse_rep, 'rep_length': len(forward_rep) + len(reverse_rep)}

    def _annotate_microhomology(self, variant) :
        # "mh" stands for microhomology
        genome_length = len(self.reference_reader.genome_seq[variant.chrom])
        if variant.mut_type == 'Ins' :
            mut_seq = variant.alt
            mut_length = len(variant.alt)
            forward_mh_start = variant.start + 1
            reverse_mh_end   = variant.start + 1
        elif variant.mut_type == 'Del' :
            mut_seq = variant.ref
            mut_length = len(variant.ref)
            forward_mh_start = variant.start + mut_length
            reverse_mh_end   = variant.start
        else :
            raise ValueError(f'Invalid mutation type: {variant.mut_type}')
        if mut_length == 1 :
            return {'forward_mh': '', 'reverse_mh': '', 'mh_seq': '', 'mh_length': 0}
        forward_mh = mut_seq
        forward_mh_length = mut_length
        while forward_mh_length > 0 :
            if (forward_mh == 
                self.reference_reader.get_sequence(variant.chrom, forward_mh_start, 
                                                   forward_mh_start + forward_mh_length)) :
                break
            forward_mh_length -= 1
            forward_mh = forward_mh[:forward_mh_length]

        reverse_mh = mut_seq
        reverse_mh_length = mut_length
        while reverse_mh_length > 0 :
            if (reverse_mh == 
                self.reference_reader.get_sequence(variant.chrom, reverse_mh_end - reverse_mh_length, 
                                                   reverse_mh_end)) :
                break
            reverse_mh_length -= 1
            reverse_mh = reverse_mh[1:]

        assert(forward_mh_length >= 0)
        assert(reverse_mh_length >= 0)
        if forward_mh_length > reverse_mh_length :
            mh_seq = forward_mh
            mh_length = forward_mh_length
        else :
            mh_seq = reverse_mh
            mh_length = reverse_mh_length
        return {'mh_seq': mh_seq, 'forward_mh': forward_mh, 'reverse_mh': reverse_mh, 'mh_length': mh_length}

    
if __name__ == '__main__': 
    import sys
    import pax_io
    import reference_reader
    pax_result_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]
    pax_reader = pax_io.PaxReader(pax_result_file)
    ref_reader = reference_reader.ReferenceReader(fasta_file)
    pax_writer = pax_io.PaxWriter(output_file, pax_reader.header + ['annotation_label'])
    pax_writer.write_head()
    cont_annot = ContextAnnotator(ref_reader)
    for variant in pax_reader :
        annotated_variant = cont_annot.annotate(variant)
        pax_writer.write(annotated_variant)
        if variant.variant_information['Spectrum'] == 'NA' :
            continue
        spectrum_simple = ':'.join(variant.variant_information['Spectrum'].split(':')[1:])
        assert(spectrum_simple == annotated_variant.label)
        
        """
        print('\t'.join([variant.chrom, str(variant.start), str(variant.end),
                         variant.ref, variant.alt,
                         variant.variant_information['Spectrum'],
                         annotated_variant.annotation['label']]))
        """    
        """ 
        rep_annot = cont_annot._annotate_repeat_length(variant)
        mh_annot  = cont_annot._annotate_microhomology(variant)
        print('\t'.join([variant.chrom, str(variant.start), str(variant.end),
                         variant.ref, variant.alt,
                         variant.variant_information['Spectrum'],
                         rep_annot['forward_rep'], rep_annot['reverse_rep'],
                         str(rep_annot['rep_length']),
                         mh_annot['reverse_mh'], mh_annot['forward_mh'],
                         str(mh_annot['mh_length'])]))
         """   
    


