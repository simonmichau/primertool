import re
import hgvs.parser
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.dataproviders.uta
import hgvs.exceptions
import logging
import mysql.connector
from mysql.connector import errorcode

from primertool_v1.exceptions import PrimertoolInputError, PrimertoolHGVSError, PrimertoolGenomeError


def query_database(config, query):
    """ Query a SQL database using a config and specific query.

    Using mysql.connector with a dict as config, see
    https://dev.mysql.com/doc/connector-python/en/connector-python-example-connecting.html

    Args:
        config: dict with user/password/host&database
        query: str

    Returns: list

    """
    try:
        cnx = mysql.connector.connect(**config)
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)

    cursor = cnx.cursor()
    cursor.execute(query)
    query_result = cursor.fetchall()
    if query_result:
        return query_result
    else:
        logging.info('This database query did not return any results. Please check your input.')
        return None


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
        raise PrimertoolInputError('Could not parse the input. There is a problem with the hgvs nomenclature: '
                                   '{0}'.format(msg))
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

    return dict(exon_number=exon_number,
                mut_start=mut_start,
                mut_end=mut_end,
                mut_length=mut_end - mut_start,
                is_in_exon=is_in_exon,
                exon_len=exon_len)


def calculate_targets(target_start, target_end, primer_bases):
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

    target_base = target_start - seq_start
    target_length = target_end - target_start
    size_range = [target_length, int(target_length + primer_bases / 2)]
    target_info = dict(target_start=target_start,
                       target_end=target_end,
                       seq_start=seq_start,
                       seq_end=seq_end,
                       target_base=target_base,
                       target_length=target_length,
                       size_range=size_range
                       )
    return target_info


def get_snps(chromosome, seq_start, seq_end, ucsc_config):
    """Retrieve common SNPs from the UCSC database.

    Query the UCSC database with the chromosome and position data to find common snps in the sequence, which need to
    be masked. Snps are returned as list.

    Args:
        chromosome: str
        seq_start: int
        seq_end: int
        ucsc_config: dict

    Returns: list

    """
    snpquery = f"SELECT chromStart FROM snp150Common WHERE chrom='{chromosome}' AND class='single' AND chromEnd BETWEEN'{seq_start}'AND'{seq_end}'"
    snp = query_database(ucsc_config, snpquery)
    if snp is None:
        return None
    else:
        snps = []
        for idx, i in enumerate(snp):
            snps.append(seq_end - int(str(snp[idx])[1:-2]))
        return snps


def mask_snps(genome, chromosome, seq_start, seq_end, ucsc_config):
    """ Mask common SNP positions with an N in the sequence.

    Using the get_snps() function to retrieve common SNPs in the given sequence from UCSC. Extract the genomic sequence
    between the given positions via genomepy and replace the bases at common SNP positions with an N.

    Args:
        genome: genomepy object
        chromosome: str
        seq_start: int
        seq_end: int
        ucsc_config: dict

    Returns: list

    """
    snps = get_snps(chromosome, seq_start, seq_end, ucsc_config)
    sequence = str(genome[chromosome][seq_start:seq_end]).upper()
    if snps is None:
        seq_snps = sequence
    else:
        snp_seq = list(sequence)
        for j in snps:
            snp_seq[j - 2] = 'N'
        seq_snps = ''.join(snp_seq)

    return seq_snps


def write_output_file(outfile, primer_strings):
    """ Write primer pairs to outfile.

    Args:
        outfile: str
        primer_strings: list of strings

    """
    logging.info(f'Writing output to {outfile}')
    with open(outfile, 'w') as f:
        for item in primer_strings:
            for x in item:
                f.write(x)


def mutalyzer_error_handler(response):
    """ Checks for errors in the mutalyzer response and raises an exception if there is an error. """

    if 'message' in response and 'custom' in response:
        # If the entries at the top level of the response are message and custom, there is a problem with the input
        logging.info(response['message'])

        # Handle infos and errors
        if 'infos' in response['custom']:
            for info in response['custom']['infos']:
                logging.info(f'{info["code"]}: {info["details"]}')  # print all infos
        if 'errors' in response['custom']:
            for error in response['custom']['errors']:
                logging.error(f'{error["code"]}: {error["details"]}')  # print all errors

            error_code = response['custom']['errors'][0]['code']
            error_message = response['custom']['errors'][0]['details']

            if error_code == 'EPARSE':
                raise PrimertoolInputError('There is an error in the given mutation', error_code, error_message)
            elif error_code == 'ERETR':
                raise PrimertoolInputError('The given NM number has an error and could not be found', error_code,
                                           error_message)
            elif error_code == 'ENOINTRON':
                raise PrimertoolInputError('The given NM number has an error and could not be found', error_code,
                                           error_message)
            elif error_code == 'ESYNTAXUC':
                raise PrimertoolInputError(error_code, error_message)
            else:
                raise PrimertoolInputError('There was a problem with the input. ', error_code, error_message)
