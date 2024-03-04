import genomepy
import genomepy.exceptions
import primer3
import os
import math
import datetime
import requests
import pandas as pd
from typing import Tuple

from . import exceptions
from . import functions
from . import logger
from .insilicopcr import InSilicoPCR

logger = logger.init_logger(level=logger.logging.INFO)


class PrimerPair(object):

    def __init__(self, primer_info):
        self.primer_info = primer_info[0]  # primer3 output

        self.mt = (self.primer_info['PRIMER_RIGHT_0_TM'] + self.primer_info['PRIMER_LEFT_0_TM']) / 2
        self.bp = self.primer_info['PRIMER_PAIR_0_PRODUCT_SIZE']

        self.chromosome = None
        self.ordertable = pd.DataFrame()
        self.orderprimer_forwards = None
        self.orderprimer_reverse = None
        self.primer_forwards = None
        self.primer_reverse = None


class ExonPrimerPair(PrimerPair):

    def __init__(self, gene_info, primer_info, exon_number: int):
        super().__init__(primer_info)
        self.gene_info = gene_info  # Gene information
        self.gene_name = gene_info['name']
        self.chromosome = gene_info['chromosome']
        self.nm_number = gene_info['nm_number']
        self.strand = gene_info['strand']

        self.exon_number = exon_number
        self.orderprimer_forwards, self.orderprimer_reverse, self.primer_forwards, self.primer_reverse = self.get_order_primers()

        # create pandas DataFrame for order table
        self.ordertable = self.make_order_table()

    def get_order_primers(self) -> Tuple[str, str, str, str]:
        """ Returns the forward and reverse primers in Gene-E(Exonnumber)F;Sequence and Gene-E(Exonnumber)R;Sequence format for ordering. """
        exon_str = str(self.exon_number).zfill(2)

        if self.strand == '+':
            orderprimer_forwards = f'{self.gene_name}-E{exon_str}F;{self.primer_info["PRIMER_LEFT_0_SEQUENCE"]}'
            orderprimer_reverse = f'{self.gene_name}-E{exon_str}R;{self.primer_info["PRIMER_RIGHT_0_SEQUENCE"]}'
            primer_forwards = self.primer_info["PRIMER_LEFT_0_SEQUENCE"]
            primer_reverse = self.primer_info["PRIMER_RIGHT_0_SEQUENCE"]
        else:
            orderprimer_forwards = f'{self.gene_name}-E{exon_str}F;{self.primer_info["PRIMER_RIGHT_0_SEQUENCE"]}'
            orderprimer_reverse = f'{self.gene_name}-E{exon_str}R;{self.primer_info["PRIMER_LEFT_0_SEQUENCE"]}'
            primer_forwards = self.primer_info["PRIMER_RIGHT_0_SEQUENCE"]
            primer_reverse = self.primer_info["PRIMER_LEFT_0_SEQUENCE"]
        return orderprimer_forwards, orderprimer_reverse, primer_forwards, primer_reverse

    def make_order_table(self) -> pd.DataFrame:
        date = datetime.datetime.now().strftime("%d.%m.%Y")
        df_order_table = pd.DataFrame({
            'date': [date, date],
            'person': [None, None],
            'primer': [self.orderprimer_forwards, self.orderprimer_reverse],
            'gene': [self.gene_info['name'], self.gene_info['name']],
            'nm_number': [self.nm_number, self.nm_number],
            'mt': [self.mt, self.mt],
            'bp': [self.bp, self.bp]
        })
        return df_order_table


class GenomicPositionPrimerPair(PrimerPair):

    def __init__(self, primer_info, chromosome: str, start: int, end: int):
        super().__init__(primer_info)
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.orderprimer_forwards, self.orderprimer_reverse, self.primer_forwards, self.primer_reverse = self.get_order_primers()

        # create pandas DataFrame for order table
        self.ordertable = self.make_order_table()

    def get_order_primers(self) -> Tuple[str, str, str, str]:
        """ Returns the forward and reverse primers in ChrStartF;Sequence and ChrStartR;Sequence format for ordering."""
        orderprimer_forwards = f'{self.chromosome}-{self.start}F;{self.primer_info["PRIMER_LEFT_0_SEQUENCE"]}'
        orderprimer_reverse = f'{self.chromosome}-{self.end}R;{self.primer_info["PRIMER_RIGHT_0_SEQUENCE"]}'
        primer_forwards = self.primer_info["PRIMER_LEFT_0_SEQUENCE"]
        primer_reverse = self.primer_info["PRIMER_RIGHT_0_SEQUENCE"]
        return orderprimer_forwards, orderprimer_reverse, primer_forwards, primer_reverse

    def make_order_table(self) -> pd.DataFrame:
        date = datetime.datetime.now().strftime("%d.%m.%Y")
        df_order_table = pd.DataFrame({
            'date': [date, date],
            'person': [None, None],
            'primer': [self.orderprimer_forwards, self.orderprimer_reverse],
            'gene': [None, None],
            'nm_number': [None, None],
            'mt': [self.mt, self.mt],
            'bp': [self.bp, self.bp]
        })
        return df_order_table


class PrimerGenerator(object):

    def __init__(self, genome_assembly: str, kuerzel: str, **kwargs):
        self.genome_assembly = self.check_genome_assembly(genome_assembly)

        self.genome_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'genomes')
        self.output_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'output')
        # create output directory if not already exists
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # check if genome directory exists and try to create it if not
        self.genome = self.fetch_genome()

        self.kuerzel = kuerzel

        self.max_insert = 800
        self.min_insert = 200
        self.dist_exon_borders = 40

        self.chromosome = None  # if chromosome is known later from gene information

    @staticmethod
    def check_genome_assembly(genome_assembly: str) -> str or exceptions.PrimertoolInputError:
        if genome_assembly not in ['hg38', 'hg19']:
            msg = f'Given genome assembly {genome_assembly} is invalid. Only hg19 and hg38 are accepted.'
            logger.exception(msg)
            raise exceptions.PrimertoolInputError(msg)
        else:
            return genome_assembly

    def fetch_genome(self) -> genomepy.Genome or exceptions.PrimertoolInputError or exceptions.PrimertoolGenomeError:
        """
        Check if the given genome assembly is available in the genome directory. If not, download it from UCSC.

        Returns: genomepy.Genome object

        """
        try:
            return genomepy.Genome(self.genome_assembly, genomes_dir=self.genome_dir)
        except FileNotFoundError:
            try:
                logger.info(f'Downloading genome {self.genome_assembly} from UCSC')
                genomepy.install_genome(self.genome_assembly, "UCSC", genomes_dir=self.genome_dir)
                return genomepy.Genome(self.genome_assembly, genomes_dir=self.genome_dir)
            except genomepy.exceptions.GenomeDownloadError as error_msg:
                msg = (f'Cannot download given reference {self.genome_assembly} from UCSC. '
                       f'Please check your input: {error_msg}')
                raise exceptions.PrimertoolGenomeError(msg)

    def check_insert_size(self, start: int, end: int) -> list:
        """ Check if insert size is between min/max insert size.

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

    def iterate_positions(self, positions: list) -> list:
        """ Generate primers with the given positions.

        If primer3 does not return a result, increase the sequence length in which primers can be generated.

        Args:
            positions: list

        Returns: list

        """
        output = []
        for index, sub_exon in enumerate(positions):
            pos_start, pos_end = sub_exon[0], sub_exon[1]

            primer_bases = 100
            primers_found = 0
            target_size = 0
            primers = None
            while not primers_found and target_size <= self.max_insert:
                # Design primers using primer3
                primers, target_size = self.design_primer(pos_start, pos_end, primer_bases)

                if target_size > self.max_insert:
                    primers = None
                    logger.warning('Stopping search (target size > max insert). No primers were found.')
                    break

                # if target_size is bigger than max_insert: split into chunks and try for primers again?
                # new_positions = self.check_insert_size(pos_start + primer_bases, pos_end + primer_bases)
                # primers = self.iterate_positions(new_positions, self.chromosome)
                # logging.info('No primers found and target size größer max insert')

                # Filter out all primer pairs that are not uniquely binding
                primers, invalid_flag = functions.filter_unique_primers(primers)

                # If primers were found, but all of them are not unique
                if invalid_flag:
                    # update pos_start and pos_end for next iteration (a.k.a. widen the target)
                    # note that the update size here is arbitrary
                    pos_start = pos_start + 100
                    pos_end = pos_end + 100
                    # reset primer_bases to init value
                    primer_bases = 0

                # increment loop variables
                primer_bases = primer_bases + 100
                primers_found = primers['PRIMER_PAIR_NUM_RETURNED']

                if not primers_found:  # log message if no primers were found
                    msg = (f'No primers <{primer_bases}primer bases around target position found yet, increasing '
                           f'allowed distance from target position to {primer_bases + 100} primer bases')
                    logger.debug(msg)

            if not primers['PRIMER_PAIR_NUM_RETURNED'] == 0:
                logger.debug(f'Primers found for position {pos_start}-{pos_end}.')
                output.append([primers, pos_start, pos_end, index + 1])

        logger.info(f'Found {len(output)} primers.')
        return output

    def design_primer(self, start: int, end: int, primer_bases: int) -> tuple:
        """ Design a primer pair using primer3.

        Firstly the targets are defined. Then the genomic sequence is retrieved and the most common SNPs are masked
        in the sequence.

        Args:
            start: int
            end: int
            primer_bases: int

        Returns: tuple (primer3 results as dict and insert size as int)

        """
        target_info = functions.calculate_targets(start, end, primer_bases)
        genomic_sequence = functions.mask_snps(self.genome, self.chromosome, target_info['seq_start'],
                                               target_info['seq_end'], self.genome_assembly)
        seq_dict = dict(SEQUENCE_TEMPLATE=genomic_sequence, SEQUENCE_TARGET=[target_info['target_base'],
                                                                             target_info['target_length']])
        primer3_config = dict(PRIMER_OPT_SIZE=20,
                              PRIMER_MIN_SIZE=20,
                              PRIMER_MAX_SIZE=22,
                              PRIMER_OPT_TM=60,
                              PRIMER_MIN_TM=58,
                              PRIMER_MAX_TM=62,
                              PRIMER_MAX_POLY_X=5,
                              GCCLAMP=1)
        primer3_config['PRIMER_PRODUCT_SIZE_RANGE'] = target_info['size_range']
        primer_out = primer3.bindings.designPrimers(seq_dict, primer3_config)
        return primer_out, target_info['size_range'][1]


class ExonPrimerGenerator(PrimerGenerator):

    def __init__(self, nm_number: str, exon_number: int, genome_assembly: str = 'hg38', kuerzel: str = None):
        super().__init__(genome_assembly, kuerzel)
        self.nm_number = nm_number
        self.exon_number = exon_number

        # # # create primer # # #
        # 1. Verify that the NM number is valid
        self.check_nm_number()
        # 2. Get RefSeq Information about the gene
        self.gene_info = functions.get_gene_information(self.genome_assembly, self.nm_number)
        self.chromosome = self.gene_info['chromosome']
        # 3. Check if this exon number exists in the gene
        self.check_exon()
        # 4. Get exon boundaries
        exon_start, exon_end = self.get_exon_boundaries()
        # 5.
        exon_positions = self.check_insert_size(exon_start, exon_end)
        # 6. Iterate over positions and generate primers
        list_primers = self.iterate_positions(exon_positions)
        # 7. Create ordertable
        self.ordertable = self.get_ordertable(self.gene_info, list_primers)

    def check_nm_number(self):
        """ Check if nm_number is valid. """
        if not self.nm_number.startswith('NM_'):
            msg = f'Given NM number {self.nm_number} is invalid. NM number should start with "NM_"'
            logger.exception(msg)
            raise exceptions.PrimertoolInputError(msg)

    def check_exon(self):
        """ Checking if exon exists in gene. """
        logger.info(f'Number of exons in gene {self.gene_info["name"]}: {self.gene_info["exoncount"]}')
        if self.exon_number > self.gene_info['exoncount']:
            msg = f'The given exon number ({self.exon_number}) is larger than the number of exons in gene {self.gene_info["name"]}'
            logger.exception(msg)
            raise exceptions.PrimertoolInputError(msg)

    def get_exon_boundaries(self):
        if self.gene_info['strand'] == '+':
            start = self.gene_info['exon_starts'][self.exon_number - 1]
            end = self.gene_info['exon_ends'][self.exon_number - 1]
        else:
            start = self.gene_info['exon_starts'][len(self.gene_info['exon_starts']) - self.exon_number]
            end = self.gene_info['exon_ends'][len(self.gene_info['exon_starts']) - self.exon_number]
        return start, end

    def get_ordertable(self, gene_info: dict, list_primers: list) -> pd.DataFrame:
        """ Get ordertable for exon primers.

        Args:
            gene_info: list
            list_primers: list

        """
        outfile = f'{self.output_dir}/{gene_info["name"]}_exon{self.exon_number}_primer.txt'

        exon_primer_pairs = []
        for primers in list_primers:
            exon_primer_pairs.append(ExonPrimerPair(gene_info, primers, self.exon_number))

        if len(exon_primer_pairs) == 0:
            logger.warning(f'No primers found for exon {self.exon_number} in gene {gene_info["name"]}')
            df = pd.DataFrame()
        else:
            df = pd.concat([pair.ordertable for pair in exon_primer_pairs], ignore_index=True)
            df['person'] = self.kuerzel
        return df


class GenomicPositionPrimerGenerator(PrimerGenerator):

    def __init__(self, chromosome: str, start: int, end: int, genome_assembly: str = 'hg38', kuerzel: str = None):
        super().__init__(genome_assembly, kuerzel)
        self.chromosome = chromosome
        self.start = start
        self.end = end

        # # # create primer # # #
        # 1. Check if the given chromosome is valid
        # self.check_chromosome()
        # 2. Check if the given start and end position are valid
        # self.check_positions()
        # 3. Check insert size
        positions = self.check_insert_size(self.start, self.end)
        # 4. Iterate over positions and generate primers
        list_primers = self.iterate_positions(positions)
        # 5. Create ordertable
        self.ordertable = self.get_ordertable(list_primers)

    def get_ordertable(self, list_primers: list) -> pd.DataFrame:
        """ Get ordertable for genomic position primers.

        Args:
            list_primers: list

        """
        outfile = f'{self.output_dir}/{self.chromosome}_{self.start}_{self.end}_primer.txt'

        position_primer_pairs = []
        for primers in list_primers:
            position_primer_pairs.append(GenomicPositionPrimerPair(primers, self.chromosome, self.start, self.end))

        df = pd.concat([pair.ordertable for pair in position_primer_pairs], ignore_index=True)
        df['person'] = self.kuerzel
        return df


class VariantPrimerGenerator(PrimerGenerator):

    def __init__(self, variant: str, genome_assembly: str = 'hg38', kuerzel: str = None):
        super().__init__(genome_assembly, kuerzel)

        # # # create primer # # #
        # 1. check if sequence identifier is for a coding transcript
        self.variant = self.check_variant(variant)
        # 2. check mutation with mutalyzer, get genomic position
        coding_mutation, genomic_mutation = self.check_mutation()
        self.nm_number, version = functions.split_nm(coding_mutation.ac)
        gene_info = functions.get_gene_information(self.genome_assembly, self.nm_number)
        mutation_position = functions.find_sequence_positions(gene_info['exon_starts'], gene_info['exon_ends'],
                                                              gene_info['exoncount'], gene_info['strand'],
                                                              genomic_mutation.posedit.pos)
        # 3. generate primers based on mutation position in gene/exon
        if mutation_position['is_in_exon'] and mutation_position['exon_len'] <= self.max_insert:
            logger.info('Mutation in exon, running PrimerExon')
            variant_primer = ExonPrimerGenerator(self.nm_number, mutation_position['exon_number'], self.genome_assembly,
                                                 self.kuerzel)
        elif mutation_position['is_in_exon'] and mutation_position['exon_len'] > self.max_insert:
            logger.info('Mutation is in an exon, but the exon length is bigger than the max insert size')
            raise exceptions.PrimertoolExonLengthError('Exon length is bigger than the max insert size. ')
            # variant_primer = GenomicPositionPrimerGenerator(gene_info['chromosome'], mutation_position['mut_start'],
            #                                                mutation_position['mut_end'], self.genome_assembly,
            #                                                self.kuerzel)
        else:
            logger.info('Mutation is not in an exon')
            variant_primer = GenomicPositionPrimerGenerator(gene_info['chromosome'], mutation_position['mut_start'],
                                                            mutation_position['mut_end'], self.genome_assembly,
                                                            self.kuerzel)
        # 4. set ordertable
        self.ordertable = variant_primer.ordertable

    @staticmethod
    def check_variant(variant: str):
        # correct case
        if variant.startswith('Chr') or variant.startswith('chr'):
            msg = f'Only NM numbers are supported by VariantPrimerGenerator. Use GenomicPositionPrimerGenerator instead'
            logger.exception(msg)
            raise exceptions.PrimertoolInputError(msg)
        # check if variant is valid
        if not variant.startswith('NM_'):
            msg = f'Given variant {variant} is invalid. Variant should start with "NM_"'
            logger.exception(msg)
            raise exceptions.PrimertoolInputError(msg)
        else:
            return variant

    def check_mutation(self):
        """ Checking the given mutation using the mutalyzer api and converting into a genomic position.

        Firstly running the mutalyzer name checker (runMutalyzer) to check the given mutation nomenclature while trying
        to catch any possible errors. Then the mutation is parsed using the hgvs parser and converted into a genomic
        position using the mutalyzer api.

        Returns: hgvs parser object

        """
        api_response = requests.get(f'https://mutalyzer.nl/api/normalize/{self.variant}?only_variants=false')

        # check if mutalyzer api request was successful
        if api_response.status_code == 200:
            response = api_response.json()
        else:
            errors = [f"{error['code']}: {error['details']}" for error in api_response.json()['custom']['errors']]
            raise exceptions.PrimertoolMutalyzerError(
                f'Status code: {api_response.status_code}. Mutalyzer API request failed. Errors: {errors}')

        # Error handling for mutalyzer response
        functions.mutalyzer_error_handler(response)

        coding_mutation = functions.parse_mutation(self.variant)

        if 'infos' in response:
            logger.info(response['infos'][0]['details'])
            coding_mutation.ac = response['corrected_model']['reference']['id']

        if coding_mutation.type != 'c':
            logger.warning(f'The input mutation is valid but not in coding reference!')
            raise exceptions.PrimertoolInputError('Input is not in coding reference!')

        if 'equivalent_descriptions' in response:
            genomic_mutation = functions.parse_mutation(response['equivalent_descriptions']['g'][0]['description'])
        elif 'chromosomal_descriptions' in response:
            genomic_mutation = functions.parse_mutation(response['chromosomal_descriptions'][0]['g'])
        else:
            raise exceptions.PrimertoolInputError('Could not resolve genomic mutation from mutalyzer response')

        return coding_mutation, genomic_mutation


class GenePrimerGenerator(PrimerGenerator):

    def __init__(self, nm_number: str, genome_assembly: str = 'hg38', kuerzel: str = None):
        super().__init__(genome_assembly, kuerzel)
        self.nm_number = nm_number

        # # # create primer # # #
        # 1. Verify that the NM number is valid
        self.check_nm_number()
        # 2. Get RefSeq Information about the gene
        self.gene_info = functions.get_gene_information(self.genome_assembly, self.nm_number)
        self.chromosome = self.gene_info['chromosome']
        # 3. Get number of exons in gene
        num_exons = self.gene_info['exoncount']
        # 4. Iterate over exons and create primers for each exon
        exon_primers = []
        for exon in range(1, num_exons + 1):
            exon_primers.append(
                ExonPrimerGenerator(self.nm_number, exon, self.genome_assembly, self.kuerzel).ordertable)

        self.ordertable = pd.concat(exon_primers, ignore_index=True)

    def check_nm_number(self):
        """ Check if nm_number is valid. """
        if not self.nm_number.startswith('NM_'):
            msg = f'Given NM number {self.nm_number} is invalid. NM number should start with "NM_"'
            logger.exception(msg)
            raise exceptions.PrimertoolInputError(msg)
