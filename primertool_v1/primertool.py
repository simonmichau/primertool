#!/usr/bin/env python
# encoding: utf-8
# AUTHOR: Daniela Dey (ddey@ukaachen.de)

import genomepy
import genomepy.exceptions
import os
import primer3
import math
import requests
import json
import pandas as pd

import primertool_v1.functions as functions
import primertool_v1.exceptions as exceptions
import primertool_v1.primer as primer
import primertool_v1.config as config

import logging

logger = logging.getLogger('Primertool')
logging.basicConfig(level=logging.INFO)


class Primertool(object):

    def __init__(self, reference='hg38', max_insert=800, min_insert=200, dist_exon_borders=40, **kwargs):
        """ Includes all necessary functions for all types of input to generate primers.

        Args:
            reference: str (hg19/hg38)
            max_insert: int
            min_insert: int
            dist_exon_borders: int
        """
        self.reference = reference
        self.output_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'output')
        self.genome_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'genomes')
        self.genome = None

        # general ucsc sql database config
        self.ucsc_config = dict(user='genome',
                                password='',
                                host='genome-euro-mysql.soe.ucsc.edu',
                                database=self.reference,
                                raise_on_warnings=True,
                                )
        self.mutalyzer_url = 'https://mutalyzer.nl/api'
        # basic primer config as dictionary
        self.primer_config = dict(PRIMER_OPT_SIZE=20,
                                  PRIMER_MIN_SIZE=20,
                                  PRIMER_MAX_SIZE=22,
                                  PRIMER_OPT_TM=60,
                                  PRIMER_MIN_TM=58,
                                  PRIMER_MAX_TM=62,
                                  PRIMER_MAX_POLY_X=5,
                                  GCCLAMP=1,
                                  )
        self.max_insert = max_insert
        self.min_insert = min_insert
        self.dist_exon_borders = dist_exon_borders

        # create output directory if not already exists
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def check_genome(self):
        """Check for the reference genome.

        Check if the reference genome is downloaded via genome-py. If not, download the reference genome.
        """
        if self.reference not in ['hg38', 'hg19']:
            logger.exception('This tool is designed to create primer for the human reference genome. '
                             'Please choose hg19 or hg38 as reference.')
            raise exceptions.PrimertoolInputError('Reference can only be hg19 or hg38.')
        else:
            try:
                genome = genomepy.Genome(self.reference, genomes_dir=self.genome_dir)
                self.genome = genome
            except FileNotFoundError:
                try:
                    genomepy.install_genome(self.reference, "UCSC", genomes_dir=self.genome_dir)
                    genome = genomepy.Genome(self.reference, genomes_dir=self.genome_dir)
                    self.genome = genome
                except genomepy.exceptions.GenomeDownloadError as msg:
                    raise exceptions.PrimertoolGenomeError(f'Cannot download given reference {self.reference} '
                                                           f'from UCSC. Please check your input- {msg}')

    def get_gene_information(self, nm_number):
        """ Retrieve gene information from RefSeq database.

        Query to the RefSeq database using the  transcript number (nm_number) to get the gene name, strand, exons etc.
        Only takes the coding transcript number without a version number!.

        Returning the query result as dict, dismissing alternative chromosomes by checking the chromosome name
        is smaller than 6 characters.

        Args:
            nm_number: str

        Returns: dict

        """

        logging.info('Retrieving gene information from RefSeq database')
        query = f'SELECT chrom, strand, name2, exonCount, cdsStart, cdsEnd, exonStarts, exonEnds ' \
                f'FROM refGene WHERE name="{nm_number}"'
        query_results = functions.query_database(self.ucsc_config, query)

        if not query_results:
            raise exceptions.PrimertoolInputError(f'Could not find gene information for {nm_number} in RefSeq database')

        if len(query_results) > 1:
            for entry in query_results:
                if len(entry[0]) < 6:
                    result = entry
        else:
            result = query_results[0]

        gene = dict(nm_number=nm_number,
                    chromosome=result[0],
                    strand=result[1],
                    name=result[2],
                    exoncount=result[3],
                    cds_start=result[4],
                    cds_end=result[5],
                    exon_starts=[int(x) for x in result[6].decode("utf-8").split(',') if x],
                    exon_ends=[int(x) for x in result[7].decode("utf-8").split(',') if x]
                    )
        return gene

    def check_insert_size(self, start, end):
        """ Check insert size is between min/max insert size.

        The insert size for the primer needs to be between min and max insert size, given the sequencing method.
        If the insert size is smaller, the difference to the min insert size is added/subtracted from start/end
        position.
        If the insert size is bigger than the max insert size, the range is split into multiple chunks,


        Args:
            start: int
            end: int

        Returns: list

        """
        insert_len = end - start
        positions = []
        if insert_len < self.min_insert:
            exon_start = int(start - (self.min_insert - insert_len) / 2)
            exon_end = int(end + (self.min_insert - insert_len) / 2)
            positions.append([exon_start, exon_end])
        elif insert_len > self.max_insert:
            num = math.ceil(insert_len / (math.ceil(insert_len / self.max_insert)))
            lst = list(range(start, end))
            chunks = [lst[i:i + num] for i in range(0, len(lst), num)]
            for i in chunks:
                positions.append([i[0], i[-1]])
        else:
            start = start - self.dist_exon_borders
            end = end + self.dist_exon_borders
            positions.append([start, end])
        return positions

    def iterate_positions(self, positions, chromosome):
        """ Generate primers with the given positions.

        If primer3 does not return a result, increase the sequence length in which primers can be generated.

        Args:
            positions: list
            chromosome: int

        Returns: list

        """
        output = []
        for index, sub_exon in enumerate(positions):
            pos_start = sub_exon[0]
            pos_end = sub_exon[1]

            primer_bases = 100
            primers, target_size = self.design_primer(chromosome, pos_start, pos_end, primer_bases)
            while primers['PRIMER_PAIR_NUM_RETURNED'] == 0:
                logging.info('No primers were found yet, increasing sequence size')
                primer_bases = primer_bases + 100

                if target_size <= self.max_insert:
                    primers, target_size = self.design_primer(chromosome, pos_start, pos_end, primer_bases)

                else:
                    # if target_size is bigger than max_insert: split into chunks and try for primers again?
                    # new_positions = self.check_insert_size(pos_start + primer_bases, pos_end + primer_bases)
                    # primers = self.iterate_positions(new_positions, chromosome)
                    # logging.info('No primers found and target size größer max insert')
                    break

            if primers['PRIMER_PAIR_NUM_RETURNED'] == 0:
                logging.info('No primers could be found')
            else:
                output.append([primers, pos_start, pos_end, index + 1])
        return output

    def design_primer(self, chromosome, start, end, primer_bases):
        """ Design a primer pair using primer3.

        Firstly the targets are defined. Then the genomic sequence is retrieved and the most common SNPs are masked
        in the sequence.

        Args:
            chromosome: str
            start: int
            end: int
            primer_bases: int

        Returns: tuple (primer3 results as dict and insert size as int)

        """
        target_info = functions.calculate_targets(start, end, primer_bases)
        genomic_sequence = functions.mask_snps(self.genome, chromosome, target_info['seq_start'],
                                               target_info['seq_end'], self.ucsc_config)
        seq_dict = dict(SEQUENCE_TEMPLATE=genomic_sequence, SEQUENCE_TARGET=[target_info['target_base'],
                                                                             target_info['target_length']])
        self.primer_config['PRIMER_PRODUCT_SIZE_RANGE'] = target_info['size_range']
        primer_out = primer3.bindings.designPrimers(seq_dict, self.primer_config)
        return primer_out, target_info['size_range'][1]


class PrimerMutation(Primertool):

    def __init__(self, mutation, reference, **kwargs):
        """ Generating primer for a mutation in hgvs nomenclature.

        Args:
            mutation: str
            reference: str (hg19/hg38)
        """
        Primertool.__init__(self, reference, **kwargs)
        self.mutation = mutation

    def check_mutation(self):
        """ Checking the given mutation using the mutalyzer api and converting into a genomic position.

        Firstly running the mutalyzer name checker (runMutalyzer) to check the given mutation nomenclature while trying
        to catch any possible errors. Then the mutation is parsed using the hgvs parser and converted into a genomic
        position using the mutalyzer api.

        Returns: hgvs parser object

        """
        api_request = requests.get(f'{self.mutalyzer_url}/normalize/{self.mutation}?only_variants=false')

        # check if mutalyzer api request was successful
        if api_request.status_code == 200:
            response = api_request.json()
        else:
            raise exceptions.PrimertoolMutalyzerError(f'Could not connect to mutalyzer api. '
                                                      f'Status code: {api_request.status_code}')

        # write mutalyzer response to json file [for debugging purposes, TODO: remove]
        with open(f'{self.mutation}.json', 'w') as outfile:
            json.dump(response, outfile)

        functions.mutalyzer_error_handler(response)  # Error handling for mutalyzer response

        coding_mutation = functions.parse_mutation(self.mutation)  # TODO: fix behavior when given non-coding mutations

        # This was a workaround to handle missing version numbers in the mutation input and automatically find a
        # corrected accession number (ac). Currently, this case is handled by the mutalyzer_error_handler, which makes a
        # logging info when it was able to infer the most recent version
        if 'infos' in response:
            logging.info(response['infos'][0]['details'])
            coding_mutation.ac = response['corrected_model']['reference']['id']

        if coding_mutation.type is not 'c':
            logging.warning(f'The input mutation is valid but not in coding reference!')
            raise exceptions.PrimertoolInputError('Input is not in coding reference!')

        if 'equivalent_descriptions' in response:
            genomic_mutation = functions.parse_mutation(response['equivalent_descriptions']['g'][0]['description'])
        elif 'chromosomal_descriptions' in response:
            genomic_mutation = functions.parse_mutation(response['chromosomal_descriptions'][0]['g'])
        else:
            raise exceptions.PrimertoolInputError('Could not resolve genomic mutation from mutalyzer response')

        return coding_mutation, genomic_mutation

    def write_outfile(self, gene_info, nm_number, posedit, list_primers):
        """ Writing a file with the found primers.

        Args:
            gene_info: list
            nm_number: str
            posedit:
            list_primers:

        Returns: str

        """
        outfile = '{0}/{1}_exon{2}_primer.txt'.format(self.output_dir, gene_info['name'], posedit)
        primer_strings = []
        primers_forward, primers_reverse = [], []
        for item in list_primers:
            epp = primer.ExonPrimerPair(gene_info, item[0], posedit, item[3])
            # TODO: is this correct? I dont think primer_output_exon takes the posedit as an argument. Also there is primer_output_posedit, which is not used at all
            #primer_string, primer_forward, primer_reverse = functions.primer_output_exon(gene_info['name'],
            #                                                                             nm_number, item[0],
            #                                                                             gene_info['strand'],
            #                                                                             posedit, item[3])
            primer_strings.append(epp.output)
            primers_forward.append(epp.primer_forwards)
            primers_reverse.append(epp.primer_reverse)
        functions.write_output_file(outfile, primer_strings)
        return outfile

    def create_primer(self):
        """ Create primers for a mutation in hgvs nomenclature.

        The following steps are taken:
            - checking if genome is downloaded
            - retrieving mutation information
            - generating primers based on position of mutation in gene/exon
        """
        # 1. check if genome is available
        self.check_genome()
        # 1b. check if sequence identifier is for a coding transcript
        self.check_sequence_identifier()
        # 2. check mutation with mutalyzer, get genomic position
        coding_mutation, genomic_mutation = self.check_mutation()
        nm_number, version = functions.split_nm(coding_mutation.ac)
        gene_info = self.get_gene_information(nm_number)
        mutation_position = functions.find_sequence_positions(gene_info['exon_starts'], gene_info['exon_ends'],
                                                              gene_info['exoncount'], gene_info['strand'],
                                                              genomic_mutation.posedit.pos)
        # 3. generate primers based on mutation position in gene/exon
        if mutation_position['is_in_exon'] and mutation_position['exon_len'] <= self.max_insert:
            logging.info('Mutation in exon, running PrimerExon')
            x = PrimerExon(nm_number, mutation_position['exon_number'], self.reference)
            primer_output = x.create_primer()
        elif mutation_position['is_in_exon'] and mutation_position['exon_len'] > self.max_insert:
            # TODO: test this case, it is not properly implemented yet
            logging.info('Mutation is in an exon, but the exon length is bigger than the max insert size')
            x = PrimerGenomicPosition(gene_info['chromosome'], mutation_position['mut_start'],
                                      mutation_position['mut_end'], self.reference)
            primer_output = x.create_primer(write_file=False)
            self.write_outfile(gene_info, nm_number, coding_mutation.posedit, primer_output)
        else:
            logging.info('Mutation is not in an exon')
            x = PrimerGenomicPosition(gene_info['chromosome'], mutation_position['mut_start'],
                                      mutation_position['mut_end'], self.reference)
            primer_output = x.create_primer(write_file=False)
            self.write_outfile(gene_info, mutation_position['mut_start'], mutation_position['mut_end'], primer_output)

        return primer_output

    def check_sequence_identifier(self):
        # check if sequence identifier is for a RefSeq coding transcript (denoted by NM_)
        if not self.mutation.startswith('NM_'):
            raise exceptions.PrimertoolInputError(f'Given input "{self.mutation}" is not a RefSeq coding transcript! '
                                                  f'Make sure to provide an NM_ identifier.')
        else:
            logging.info(f'Valid Sequence Identifier in {self.mutation}')


class PrimerExon(Primertool):

    def __init__(self, nm_number, exon, reference, **kwargs):
        """ Generating primer for a specific exon of a transcript.

        Args:
            nm_number: str
            exon: int
            reference: str (hg19/hg38)
        """
        Primertool.__init__(self, reference, **kwargs)

        self.nm_number = nm_number
        self.exon = int(exon)

    def check_exon(self, gene):
        """ Checking if exon exists in gene.

        Args:
            gene: list
        """
        logging.info('number of exons ' + str(gene['exoncount']))
        if self.exon > gene['exoncount']:
            logging.exception(f'This exon number is not in number of exons in gene {self.exon}')
            raise exceptions.PrimertoolInputError(f'Exon number not in number of exons in gene {self.exon}')
        else:
            logging.info('This exon number exists in this gene.')

    def write_outfile(self, gene_info, list_primers):
        """ Write primers to file.

        Args:
            gene_info: list
            list_primers: list

        """
        outfile = '{0}/{1}_exon{2}_primer.txt'.format(self.output_dir, gene_info['name'], self.exon)
        epp_list = []
        for item in list_primers:
            epp = primer.ExonPrimerPair(gene_info, item[0], self.exon, item[3])
            epp_list.append(epp)

        primer_list = [pp.output for pp in epp_list]

        primer_list = []
        for pp in epp_list:
            primer_list.append(pp.output)

        functions.write_output_file(outfile, primer_list)

        # Ordertable
        ordertable = [epp.get_ordertable() for epp in epp_list]
        if len(ordertable) > 1:
            ordertable = pd.concat(ordertable, ignore_index=True)
        elif len(ordertable) == 1:
            ordertable = ordertable[0]
        else:
            ordertable = pd.DataFrame()
        return ordertable

    def create_primer(self, return_ordertable=True):
        """ Create primers for a specific exon.

        The following steps are taken:
            - checking if genome is downloaded
            - retrieving mutation information
            - generating primers based on position of mutation in gene/exon

        Args:
            write_file: boolean

        """
        # 1. Check if genome build is available
        self.check_genome()
        # 2. Get RefSeq Information about the gene
        gene_info = self.get_gene_information(self.nm_number)
        # 3. Check if this exon number exists in the gene
        self.check_exon(gene_info)
        # 4. Define exon boundaries
        if gene_info['strand'] == '+':
            exon_start = gene_info['exon_starts'][self.exon - 1]
            exon_end = gene_info['exon_ends'][self.exon - 1]
        else:
            exon_start = gene_info['exon_starts'][len(gene_info['exon_starts']) - self.exon]
            exon_end = gene_info['exon_ends'][len(gene_info['exon_starts']) - self.exon]

        exon_positions = self.check_insert_size(exon_start, exon_end)

        list_primers = self.iterate_positions(exon_positions, gene_info['chromosome'])
        logging.info(f'Primers generated. {str(len(list_primers))} primers found. ')

        ordertable = self.write_outfile(gene_info, list_primers)

        if not return_ordertable:
            return list_primers

        return ordertable


class PrimerGen(Primertool):

    def __init__(self, nm_number, reference, **kwargs):
        """ Generating primer for all exons of specific transcript.

        Args:
            nm_number: str
            reference: str (hg19/hg38)
        """

        Primertool.__init__(self, reference, **kwargs)

        self.nm_number = nm_number

    def write_outfile(self, gene_info, out_strings):
        """ Write primers to file.

        Args:
            gene_info: list
            out_strings: list

        Returns:

        """
        outfile = f'{self.output_dir}/{gene_info["name"]}_primer.txt'
        functions.write_output_file(outfile, out_strings)
        return outfile

    def create_primer(self):
        """ Create primers for all exons in a transcript."""
        # 1. Check if genome build is available
        self.check_genome()
        # 2. Get RefSeq Information about the gene
        gene_info = self.get_gene_information(self.nm_number)

        out_strings = []
        for exon in range(0, gene_info['exoncount']):
            logging.info(f'Creating primers for exon {exon + 1}')
            x = PrimerExon(self.nm_number, exon + 1, self.reference)
            list_primers = x.create_primer(return_ordertable=False)
            for primer_pair in list_primers:
                epp = primer.ExonPrimerPair(gene_info, primer_pair[0], exon + 1, primer_pair[3])
                out_strings.append(epp.output)
                # out_strings.append(functions.primer_output_exon(gene_info['name'], self.nm_number, primer_pair[0],
                #                                                gene_info['strand'], exon + 1, primer_pair[3]))

        ordertable = self.write_outfile(gene_info, out_strings)
        return ordertable


class PrimerGenomicPosition(Primertool):

    def __init__(self, chromosome, start, end, reference, **kwargs):
        """ Primer for a genomic position.

        Args:
            chromosome: str
            start: int
            end: int
            reference: str (hg19/hg38)
        """
        Primertool.__init__(self, reference, **kwargs)
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def write_outfile(self, list_primers):
        """ Write primers to file.

        Args:
            list_primers: list

        """
        outfile = '{0}/{1}_{2}_{3}_primer.txt'.format(self.output_dir, self.chromosome, self.start, self.end)
        gppp_list = []
        for item in list_primers:
            gppp = primer.GenomicPositionPrimerPair(item[0], self.chromosome, item[1], item[2], item[3])
            gppp_list.append(gppp.output)
            #primer_strings.append(functions.primer_output_genomic(item[0], self.chromosome, item[1], item[2], item[3]))
        functions.write_output_file(outfile, gppp_list)

        #primer_list = []
        #for pp in gppp_list:
        #    primer_list.append(pp.output)
        #functions.write_output_file(outfile, primer_list)

        # Ordertable
        ordertable = [gppp.get_ordertable() for gppp in gppp_list]
        if len(ordertable) > 1:
            ordertable = pd.concat(ordertable, ignore_index=True)
        return ordertable

    def create_primer(self, write_file=True):
        """ Create primers for a specific genomic position.

        Args:
            write_file: boolean

        """
        # 1. Check if genome build is available
        self.check_genome()

        positions = self.check_insert_size(self.start, self.end)
        list_primers = self.iterate_positions(positions, self.chromosome)
        if write_file:
            self.write_outfile(list_primers)
        else:
            return list_primers
