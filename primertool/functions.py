from typing import Tuple
import genomepy
import re
import hgvs.parser
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.dataproviders.uta
import hgvs.exceptions
import mysql.connector
from mysql.connector import errorcode
import logging

from . import exceptions
from .insilicopcr import InSilicoPCR

logger = logging.getLogger(__package__)


def query_ucsc_database(genome_assembly: str, query: str):
    """ Query a SQL database using a config and specific query.

    Using mysql.connector with a dict as config, see
    https://dev.mysql.com/doc/connector-python/en/connector-python-example-connecting.html

    Args:
        genome_assembly:
        query: str

    Returns: list

    """
    # general ucsc sql database config
    ucsc_config = dict(user='genome',
                       password='',
                       host='genome-euro-mysql.soe.ucsc.edu',
                       database=genome_assembly,
                       raise_on_warnings=True,
                       )
    query_result = None
    try:
        cnx = mysql.connector.connect(**ucsc_config)
        cursor = cnx.cursor()
        cursor.execute(query)
        query_result = cursor.fetchall()
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            logger.error("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            logger.error("Database does not exist")
        else:
            logger.error(err)

    if query_result:
        return query_result
    else:
        logger.debug('This database query did not return any results. Please check your input.')
        return None


def get_gene_information(genome_assembly: str, nm_number: str):
    """ Retrieve gene information from RefSeq database for a given NM number. """

    logger.info('Collecting gene information from RefSeq database')
    query = f'SELECT chrom, strand, name2, exonCount, cdsStart, cdsEnd, exonStarts, exonEnds ' \
            f'FROM refGene WHERE name="{nm_number}"'

    query_results = query_ucsc_database(genome_assembly, query)

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


def calculate_targets(target_start: int, target_end: int, primer_bases: int) -> dict:
    """ Defining the sequence start and end positions.

    Primer bases are the number of bases to each side of the target sequence in which primer3 looks for a possible
    primer. These are added to sequence start and end respectively to calculate the sequence start and end positions.
    The size range is a list with the length of the target sequence and target sequence with the primer bases added.
    Primer3 should design primers in this size range.


    Args:
        target_start: int
        target_end: int
        primer_bases: int

    Returns: dict

    """
    seq_start = target_start - primer_bases
    seq_end = target_end + primer_bases
    target_base = primer_bases
    target_length = target_end - target_start
    size_range = [target_length, target_length + int(primer_bases / 2)]

    target_info = dict(target_start=target_start,
                       target_end=target_end,
                       seq_start=seq_start,
                       seq_end=seq_end,
                       target_base=target_base,
                       target_length=target_length,
                       size_range=size_range)
    return target_info


def mask_snps(genome: genomepy.Genome, chromosome: str, seq_start: int, seq_end: int, genome_assembly: str):
    """ Mask common SNP positions with an N in the sequence.

    Using the get_snps() function to retrieve common SNPs in the given sequence from UCSC. Extract the genomic sequence
    between the given positions via genomepy and replace the bases at common SNP positions with an N.

    Args:
        genome: genomepy object
        chromosome: str
        seq_start: int
        seq_end: int
        genome_assembly: str

    Returns: list

    """
    sequence = str(genome[chromosome][seq_start:seq_end]).upper()

    snps = get_snps(chromosome, seq_start, seq_end, genome_assembly)
    if snps is None:
        seq_snps = sequence
    else:
        snp_seq = list(sequence)
        for j in snps:
            snp_seq[j - 2] = 'N'
        seq_snps = ''.join(snp_seq)

    return seq_snps


def get_snps(chromosome: str, seq_start: int, seq_end: int, genome_assembly: str):
    """Retrieve common SNPs from the UCSC database.

    Query the UCSC database with the chromosome and position data to find common snps in the sequence, which need to
    be masked. Snps are returned as list.

    Args:
        chromosome: str
        seq_start: int
        seq_end: int
        genome_assembly: str

    Returns: list

    """
    snpquery = f"SELECT chromStart FROM snp150Common WHERE chrom='{chromosome}' AND class='single' AND chromEnd BETWEEN'{seq_start}'AND'{seq_end}'"
    snp = query_ucsc_database(genome_assembly, snpquery)
    snps = []
    if snp is not None:
        for idx, i in enumerate(snp):
            snps.append(seq_end - int(str(snp[idx])[1:-2]))

    return snps


def mutalyzer_error_handler(response):
    """ Checks for errors in the mutalyzer response and raises an exception if there is an error. """

    if 'message' in response and 'custom' in response:
        # If the entries at the top level of the response are message and custom, there is a problem with the input
        logger.info(response['message'])

        # Handle infos and errors
        if 'infos' in response['custom']:
            for info in response['custom']['infos']:
                logger.info(f'{info["code"]}: {info["details"]}')  # print all infos
        if 'errors' in response['custom']:
            for error in response['custom']['errors']:
                logger.error(f'{error["code"]}: {error["details"]}')  # print all errors

            error_code = response['custom']['errors'][0]['code']
            error_message = response['custom']['errors'][0]['details']

            if error_code == 'EPARSE':
                raise exceptions.PrimertoolInputError('There is an error in the given mutation', error_code,
                                                      error_message)
            elif error_code == 'ERETR':
                raise exceptions.PrimertoolInputError('The given NM number has an error and could not be found',
                                                      error_code,
                                                      error_message)
            elif error_code == 'ENOINTRON':
                raise exceptions.PrimertoolInputError('The given NM number has an error and could not be found',
                                                      error_code,
                                                      error_message)
            elif error_code == 'ESYNTAXUC':
                raise exceptions.PrimertoolInputError(error_code, error_message)
            else:
                raise exceptions.PrimertoolInputError('There was a problem with the input. ', error_code, error_message)


def parse_mutation(mutation):
    """ Parse mutation with hgvs parser.

    Genenames in brackets are removed from the variant (eg: eg: NM_003165.6(STXBP1):c.1702G>A).
    The mutation is then parsed using hgvs.parser and a hgvs tree object (see
    https://hgvs.readthedocs.io/en/stable/key_concepts.html#variant-object-representation) is returned.
    Parses coding and genomic variants.

    Args:
        mutation: str

    Returns: hgvs.parser tree object

    """
    mutation = re.sub(r'\([^)]*\)', '', mutation)

    try:
        hp = hgvs.parser.Parser()
        hgvs_mutation = hp.parse_hgvs_variant(mutation)
    except hgvs.exceptions.HGVSParseError as msg:
        raise exceptions.PrimertoolInputError(
            f'Could not parse the input. There is a problem with the hgvs nomenclature: {msg}')
    return hgvs_mutation


def split_nm(nm_number):
    """ Split variant.ac from hgvs object into transcript and version number

    If the transcript number includes a version number, split at '.'. If not, set the version number to 1.

    Args:
        nm_number: str

    Returns: tuple

    """
    nm_split = nm_number.split('.')
    transcript = nm_split[0]

    if len(nm_split) == 1:
        version = 1
    else:
        version = int(nm_split[1])

    return transcript, version


def find_sequence_positions(exon_starts, exon_ends, exon_count, strand, mutation_position):
    """ Establish if variant position is in an exon.

    The hgvs.parser object has the start and end position of the variant (in case of an indel rather than a snp). Using
    these positions and the lists of exon starts and ends, determine if the variant is located in an exon.

    The information is returned as a dictionary. If the gene is based on the - strand, the exon number needs to be
    inverted.

    Args:
        exon_starts: list
        exon_ends: list
        exon_count: int
        strand: str
        mutation_position: hgvs.location.interval

    Returns: dict

    """
    mut_start = int(str(mutation_position.start))
    mut_end = int(str(mutation_position.end))
    exon_number = 0
    exon_len = 0
    is_in_exon = False

    for exon in range(exon_count):
        if exon_starts[exon] <= mut_start and mut_end <= exon_ends[exon]:
            is_in_exon = True
            exon_len = exon_ends[exon] - exon_starts[exon]
            if strand == '-':
                exon_number = exon_count - exon
            else:
                exon_number = exon + 1

    return dict(exon_number=exon_number, mut_start=mut_start, mut_end=mut_end, mut_length=mut_end - mut_start,
                is_in_exon=is_in_exon, exon_len=exon_len)


def filter_unique_primers(primer3_dict: dict) -> Tuple[dict, bool]:
    """
    Check primer uniqueness and remove all primers with multiple binding sites, so that only the uniquely binding ones
    remain.
    Also returns a flag if all primers were invalid (i.e. not uniquely binding), but only if there was at least one
    primer to begin with.
    """
    pre_filter_primer_count = primer3_dict['PRIMER_PAIR_NUM_RETURNED']

    idx = 0
    while idx < primer3_dict['PRIMER_PAIR_NUM_RETURNED']:  # iterate over all generated primers
        forward_primer = primer3_dict[f'PRIMER_LEFT_{idx}_SEQUENCE']
        reverse_primer = primer3_dict[f'PRIMER_RIGHT_{idx}_SEQUENCE']
        if not InSilicoPCR(forward_primer, reverse_primer).is_uniquely_binding():
            # delete primerpair with index idx from primers
            primer3_dict = purge_primer_pair(primer3_dict, idx)  # TODO: test
            # don't increment idx because primer idx was removed and replaced by primer idx+1
        else:
            idx = idx + 1

    post_filter_primer_count = primer3_dict['PRIMER_PAIR_NUM_RETURNED']

    # If >0 primers were found, but all were invalid (i.e. not uniquely binding), set flag
    all_primers_invalid_flag = pre_filter_primer_count > 0 and post_filter_primer_count == 0

    return primer3_dict, all_primers_invalid_flag


def purge_primer_pair(primer3_dict: dict, index: int) -> dict:
    """
    Removes primer pair of given index from primer3_dict and updates the dictionary so that the remaining data stays
    consistent. (I.e. reset indices so enumeration does not have gaps and update)
    """
    # Remove primerpair with given index (remove every item from the dict where the index appears in the key)
    for key, value in primer3_dict:
        if re.match(f'.*{index}.*', key):
            del primer3_dict[key]

    # Reduce index of all primers with larger index than the one we removed
    for idx in range(index + 1, primer3_dict['PRIMER_PAIR_NUM_RETURNED']):
        for key, value in primer3_dict:
            # if the corresponding index is found
            if re.findall(str(idx), key):
                # Rename keys
                new_key = reduce_numbers_in_string(key)
                primer3_dict[new_key] = primer3_dict.pop(key)

    # Count down number of returned primers
    primer3_dict['PRIMER_LEFT_NUM_RETURNED'] -= 1
    primer3_dict['PRIMER_RIGHT_NUM_RETURNED'] -= 1
    primer3_dict['PRIMER_PAIR_NUM_RETURNED'] -= 1

    return primer3_dict


def reduce_numbers_in_string(input_string: str) -> str:
    """
    Takes an input string and reduces any (positive) integer in it by one.
    """
    # Define a regex pattern to match numbers
    pattern = r"\d+"

    # Find all occurrences of the pattern in the input string
    matches = re.findall(pattern, input_string)

    # Decrement each matched number by 1
    decremented_values = [str(int(match) - 1) for match in matches]

    # Replace the original numbers with the decremented values
    for i, match in enumerate(matches):
        input_string = input_string.replace(match, decremented_values[i], 1)

    return input_string
