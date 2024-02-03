#!/usr/bin/env python
# encoding: utf-8
# AUTHOR: Simon Michau (smichau@ukaachen.de)

import primertool_v1.config as config
import datetime
import pandas as pd


class Primer(object):
    """
    General class describing Primers of any type (Exon, Posedit, Genomic Position)
    """

    def __init__(self, order_notation: str, primer: dict = None, gene_info: dict = None, direction=None):
        # melting_temperature, sequence, gc, start, length, end, order_notation
        # etc. TODO
        self.direction = direction  # 'forwards' or 'reverse'
        self.mt = None
        self.seq = None
        self.gc = None
        self.start = None
        self.length = None
        self.end = None
        self.order_notation = None  # do later?

    def __str__(self):
        return self.order_notation


class PrimerPair(object):
    def __init__(self):
        self.forward = None
        self.reverse = None
        self.length = None  # length in bp

        # additional information (genename,exonnumber,index|genename,posedit,index|chromosome,start,end,index) (header1)
        self.info = dict()
        self.index = None

        # check on init if both forward and reverse primers are of the same subclass
        if not isinstance(self.forward, type(self.reverse)):  # TODO: test this
            raise TypeError(f'Both primers must be of the same primer subclass! '
                            f'Got {type(self.forward)} and {type(self.reverse)}')

    def __str__(self):
        pass

    def __repr__(self):
        pass


class ExonPrimerPair(PrimerPair):
    def __init__(self, gene_info: dict, primer: dict, exon_number: int, index: int):
        super().__init__()
        self.gene_info = gene_info  # Gene information
        self.gene_name = gene_info['name']
        self.chromosome = gene_info['chromosome']
        self.nm_number = gene_info['nm_number']
        self.strand = gene_info['strand']

        self.primer = primer  # primer3 output
        self.mt = (self.primer['PRIMER_RIGHT_0_TM'] + self.primer['PRIMER_LEFT_0_TM']) / 2
        self.bp = self.primer['PRIMER_PAIR_0_PRODUCT_SIZE']

        self.exon_number = exon_number
        self.index = index

        # Output
        header, header_2, primer_pairs, info, primer_forwards, primer_reverse = self.generate_output()
        self.output = [header, header_2, primer_pairs, info]
        self.primer_forwards = primer_forwards
        self.primer_reverse = primer_reverse

    def __str__(self):
        return self.output

    def __repr__(self):
        return self.output

    def get_ordertable(self):
        """ Returns the order table for the primer pair as a string.
        Format:
        date | kuerzel | primer | gene/fragment | transcript | melting temp | bp
        """
        date = datetime.datetime.now().strftime("%d.%m.%Y")
        kuerzel = config.KUERZEL

        # order table as string/csv
        order_table = (
            f'Datum\tAuftraggeber\tPrimer zur Bestellung\tGen/Fragment\tTranskript\tSchmelztemperatur\tbp\n'
            f'{date}\t{kuerzel}\t{self.primer_forwards}\t{self.gene_name}\t{self.nm_number}\t{self.mt}\t{self.bp}\n'
            f'{date}\t{kuerzel}\t{self.primer_reverse}\t{self.gene_name}\t{self.nm_number}\t{self.mt}\t{self.bp}\n'
        )

        # order table as pandas dataframe
        df_order_table = pd.DataFrame({
            'Datum': [date, date],
            'Auftraggeber': [kuerzel, kuerzel],
            'Primer zur Bestellung': [self.primer_forwards, self.primer_reverse],
            'Gen/Fragment': [self.gene_name, self.gene_name],
            'Transkript': [self.nm_number, self.nm_number],
            'Schmelztemperatur': [self.mt, self.mt],
            'bp': [self.bp, self.bp]
        })

        return df_order_table

    def generate_output(self):
        header = f'{self.gene_name}, Exon: {self.exon_number}, Primerpaar: {self.index}\n'
        header_2 = """{0:7}\t{1:<5}\t{2:<6}\t{3:4}\t{4:4}\t{5:35}\n""".format('', 'Start', 'Length', 'Tm', 'GC%', 'Sequence')
        exon_str = str(self.exon_number).zfill(2)

        if self.strand == '+':
            primer_pairs = (
                "Forward\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35}\n"
                "Reverse\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}\n"
                "Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}\n\n").format(**self.primer)

            info = (f"{self.gene_name}-E{exon_str}F;{self.primer['PRIMER_LEFT_0_SEQUENCE']}\n"
                    f"{self.gene_name}-E{exon_str}R;{self.primer['PRIMER_RIGHT_0_SEQUENCE']}\n"
                    f"{self.gene_name}; {round(self.mt)} °C; {self.bp}bp; {self.nm_number}\n")

            primer_forwards = f'{self.gene_name}-E{exon_str}F;{self.primer["PRIMER_LEFT_0_SEQUENCE"]}'
            primer_reverse = f'{self.gene_name}-E{exon_str}R;{self.primer["PRIMER_RIGHT_0_SEQUENCE"]}'
        else:
            primer_pairs = (
                "Forward\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}\n"
                "Reverse\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35}\n"
                "Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}\n\n").format(**self.primer)

            info = (f"{self.gene_name}-E{exon_str}F;{self.primer['PRIMER_RIGHT_0_SEQUENCE']}\n"
                    f"{self.gene_name}-E{exon_str}R;{self.primer['PRIMER_LEFT_0_SEQUENCE']}\n"
                    f"{self.gene_name}; {round(self.mt)} °C; {self.bp}bp; {self.nm_number}\n")

            primer_forwards = f'{self.gene_name}-E{exon_str}F;{self.primer["PRIMER_RIGHT_0_SEQUENCE"]}'
            primer_reverse = f'{self.gene_name}-E{exon_str}R;{self.primer["PRIMER_LEFT_0_SEQUENCE"]}'

        return header, header_2, primer_pairs, info, primer_forwards, primer_reverse


class GenomicPositionPrimerPair(PrimerPair):
    def __init__(self, primer: dict, chromosome, start, end, index: int):
        super().__init__()
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.index = index
        self.primer = primer  # primer3 output
        self.mt = (self.primer['PRIMER_RIGHT_0_TM'] + self.primer['PRIMER_LEFT_0_TM']) / 2
        self.bp = self.primer['PRIMER_PAIR_0_PRODUCT_SIZE']

        # Output
        header, header_2, primer_pairs, info, primer_forwards, primer_reverse = self.generate_output()
        self.output = [header, header_2, primer_pairs, info]
        self.primer_forwards = primer_forwards
        self.primer_reverse = primer_reverse

    def __str__(self):
        return self.output

    def __repr__(self):
        return self.output

    def get_ordertable(self):
        date = datetime.datetime.now().strftime("%d.%m.%Y")
        kuerzel = config.KUERZEL

        # order table as string/csv
        #order_table = (
        #    f'Datum\tAuftraggeber\tPrimer zur Bestellung\tGen/Fragment\tTranskript\tSchmelztemperatur\tbp\n'
        #    f'{date}\t{kuerzel}\t{self.primer_forwards}\t{self.gene_name}\t{self.nm_number}\t{self.mt}\t{self.bp}\n'
        #    f'{date}\t{kuerzel}\t{self.primer_reverse}\t{self.gene_name}\t{self.nm_number}\t{self.mt}\t{self.bp}\n'
        #)

        # TODO: find out what information should be in the ordertable in case of a GP primer
        # order table as pandas dataframe
        df_order_table = pd.DataFrame({
            'Datum': [date, date],
            'Auftraggeber': [kuerzel, kuerzel],
            'Primer zur Bestellung': [self.primer_forwards, self.primer_reverse],
            'Gen/Fragment': ['chromosome', 'chromosome'],
            'Transkript': ['transcript', 'transcript'],
            'Schmelztemperatur': [self.mt, self.mt],
            'bp': [self.bp, self.bp]
        })

        return df_order_table

    def generate_output(self):
        header = f'Chromosome: {self.chromosome}, Start: {self.start}, Ende: {self.end}, Primerpaar: {self.index}\n'
        header_2 = """{0:7}\t{1:<5}\t{2:<6}\t{3:4}\t{4:4}\t{5:35}\n""".format('', 'Start', 'Length', 'Tm', 'GC%', 'Sequence')

        primer_pairs = (
            "Forward\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35}\n"
            "Reverse\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}\n"
            "Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}\n\n").format(**self.primer)

        info = (f"{self.chromosome}-{self.start}F;{self.primer['PRIMER_LEFT_0_SEQUENCE']}\n"
                f"{self.chromosome}-{self.end}R;{self.primer['PRIMER_RIGHT_0_SEQUENCE']}\n"
                f"{self.chromosome}-{self.start}-{self.end}; {round(self.mt)} °C; {self.bp}bp;\n")

        primer_forwards = f'{self.chromosome}-{self.start}F;{self.primer["PRIMER_LEFT_0_SEQUENCE"]}'
        primer_reverse = f'{self.chromosome}-{self.end}R;{self.primer["PRIMER_RIGHT_0_SEQUENCE"]}'

        return header, header_2, primer_pairs, info, primer_forwards, primer_reverse


class PoseditPrimerPair(PrimerPair):
    def __init__(self, gene_info: dict, primer: dict, posedit, index: int):
        super().__init__()
        self.gene_name = gene_info['name']
        self.chromosome = gene_info['chromosome']
        self.nm_number = gene_info['nm_number']
        self.strand = gene_info['strand']
        self.index = index
        self.primer = primer  # primer3 output
        self.mt = (self.primer['PRIMER_RIGHT_0_TM'] + self.primer['PRIMER_LEFT_0_TM']) / 2
        self.bp = self.primer['PRIMER_PAIR_0_PRODUCT_SIZE']
        self.posedit = posedit

        # Output
        header, header_2, primer_pairs, info, primer_forwards, primer_reverse = self.generate_output()
        self.output = [header, header_2, primer_pairs, info]
        self.primer_forwards = primer_forwards
        self.primer_reverse = primer_reverse

    def __str__(self):
        return self.output

    def __repr__(self):
        return self.output

    def get_ordertable(self):
        # TODO
        pass

    def generate_output(self):
        header = f'{self.gene_name}, Position: {self.posedit}, Primerpaar: {self.index}\n'
        header_2 = """{0:7}\t{1:<5}\t{2:<6}\t{3:4}\t{4:4}\t{5:35}\n""".format('', 'Start', 'Length', 'Tm', 'GC%', 'Sequence')
        exon_str = str(self.posedit).zfill(2)

        if self.strand == '+':
            primer_pairs = (
                "Forward\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35}\n"
                "Reverse\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}\n"
                "Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}\n\n").format(**self.primer)

            info = (f"{self.gene_name}-{exon_str}F;{self.primer['PRIMER_LEFT_0_SEQUENCE']}\n"
                    f"{self.gene_name}-{exon_str}R;{self.primer['PRIMER_RIGHT_0_SEQUENCE']}\n"
                    f"{self.gene_name}; {round(self.mt)} °C; {self.bp}bp; {self.nm_number}\n")

            primer_forwards = f'{self.gene_name}-{exon_str}F;{self.primer["PRIMER_LEFT_0_SEQUENCE"]}'
            primer_reverse = f'{self.gene_name}-{exon_str}R;{self.primer["PRIMER_RIGHT_0_SEQUENCE"]}'
        else:
            primer_pairs = (
                "Forward\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}\n"
                "Reverse\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35})\n"
                "Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}\n\n").format(**self.primer)

            info = (f"{self.gene_name}-E{exon_str}F;{self.primer['PRIMER_RIGHT_0_SEQUENCE']}\n"
                    f"{self.gene_name}-E{exon_str}R;{self.primer['PRIMER_LEFT_0_SEQUENCE']})\n"
                    f"{self.gene_name}; {round(self.mt)} °C; {self.bp}bp; {self.nm_number}\n")

            primer_forwards = f'{self.gene_name}-{exon_str}F;{self.primer["PRIMER_RIGHT_0_SEQUENCE"]}'
            primer_reverse = f'{self.gene_name}-{exon_str}R;{self.primer["PRIMER_LEFT_0_SEQUENCE"]}'

        return header, header_2, primer_pairs, info, primer_forwards, primer_reverse
