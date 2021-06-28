import re
import hgvs.parser
import hgvs.assemblymapper
import hgvs.variantmapper
import hgvs.dataproviders.uta
import hgvs.exceptions
import logging
import zeep
import math
import mysql.connector
from mysql.connector import errorcode

from primertool.exceptions import PrimertoolInputError, PrimertoolHGVSError, PrimertoolGenomeError


def query_database(config, query):
    # INFO: Query a SQL database using a config and specific query
    #   Using dict for connection
    #   (https://dev.mysql.com/doc/connector-python/en/connector-python-example-connecting.html)

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
    # INFO: parse mutation to hgvs.parser, returns tree object with ac, posedit, type etc
    # logging.info('Parsing the given mutation into hgvs')
    mutation = re.sub(r'\([^)]*\)', '', mutation)
    try:
        hp = hgvs.parser.Parser()
        hgvs_mutation = hp.parse_hgvs_variant(mutation)
    except hgvs.exceptions.HGVSParseError as msg:
        raise PrimertoolInputError('Could not parse the input. There is a problem with the hgvs nomenclature: '
                                   '{0}'.format(msg))

    return hgvs_mutation


def convert_variant_notation(mutation, reference, url):
    # INFO: takes hgvs.parser object, converts to genomic reference based on in put reference genome
    # TODO: Converting to genomic ref
    #  1: Does mutation.ac have a version number? -> how to somehow find the 'newest' version number
    # logging.info('Coverting the given mutation into a genomic position')
    client = zeep.Client(url)
    var_g = client.service.numberConversion(reference, mutation)

    if not var_g[0]:
        # logging.info('First conversion did not work, trying with lower version number')
        nm, version = split_nm(str(mutation.ac))
        if version != 1:
            version = int(version) - 1
            new_nm = nm + '.' + str(version)
            mutation.ac = new_nm
            # logging.info(f'running again with {new_nm}')
            genomic = convert_variant_notation(mutation, reference, url)
        else:
            # logging.info('Mutation could not be converted into a genomic position')
            raise PrimertoolInputError('Mutation could not be converted into a genomic position')
    else:
        genomic = parse_mutation(var_g[0])
        logging.info('The genomic position of {0} is: {1}'.format(mutation, genomic))
    return genomic


def split_nm(nm_number):
    # INFO: split mutation.ac into nm_number and version number.
    #    if no version number in string, set as 0?

    nm_split = nm_number.split('.')
    number = nm_split[0]

    if len(nm_split) == 1:
        version = 1
    else:
        version = int(nm_split[1])

    return number, version


def find_sequence_positions(exon_starts, exon_ends, exon_count, strand, mutation_position):
    # INFO: check if position of mutation in an exon
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


def calculate_targets(start, end, primer_bases):
    # 1. Add bases to exon borders
    target_start = start
    target_end = end
    # 2. Define sequence start and end (area in which primers will be generated)
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


def get_snps(chromosome, seq_start, seq_end, uscsc_config):
    snpquery = f"SELECT chromStart FROM snp150Common WHERE chrom='{chromosome}' AND class='single' AND chromEnd BETWEEN'{seq_start}'AND'{seq_end}'"
    snp = query_database(uscsc_config, snpquery)
    if snp is None:
        return None
    else:
        snps = []
        for idx, i in enumerate(snp):
            snps.append(seq_end - int(str(snp[idx])[1:-2]))
        return snps


def mask_snps(genome, chromosome, seq_start, seq_end, ucsc_config):
    snps = get_snps(chromosome, seq_start, seq_end, ucsc_config)
    sequence = str(genome[chromosome][seq_start:seq_end]).upper()
    if snps is None:
        seq_snps = sequence
    else:
        snp_seq = list(sequence)
        for j in snps:
            snp_seq[j-2] = 'N'
        seq_snps = ''.join(snp_seq)

    return seq_snps


def primer_output_exon(genname, nm_number, primer, strand, exon_number, index):
    # INFO: write generated primers to file
    #   switch forward/reverse based on strand information
    header = f'{genname}, Exon: {exon_number}, Primerpaar: {index}\n'
    header_2 = """{0:7}\t{1:<5}\t{2:<6}\t{3:4}\t{4:4}\t{5:35}""".format('', 'Start', 'Length', 'Tm', 'GC%', 'Sequence')
    exon_str = str(exon_number).zfill(2)
    temp = (primer['PRIMER_RIGHT_0_TM'] + primer['PRIMER_LEFT_0_TM']) / 2

    if strand == '+':
        primer_pairs = """
Forward\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35}
Reverse\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}
Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}
""".format(**primer)

        info = f"""
{genname}-E{exon_str}F;{primer['PRIMER_LEFT_0_SEQUENCE']}
{genname}-E{exon_str}R;{primer['PRIMER_RIGHT_0_SEQUENCE']}
{genname}; {round(temp)} °C; {primer['PRIMER_PAIR_0_PRODUCT_SIZE']}bp; {nm_number}

"""
    else:
        primer_pairs = """

Forward\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}
Reverse\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35}
Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}
""".format(**primer)

        info = f"""
{genname}-E{exon_str}F;{primer['PRIMER_RIGHT_0_SEQUENCE']}
{genname}-E{exon_str}R;{primer['PRIMER_LEFT_0_SEQUENCE']}
{genname}; {round(temp)} °C; {primer['PRIMER_PAIR_0_PRODUCT_SIZE']}bp; {nm_number}

"""
    output = [header, header_2, primer_pairs, info]

    return output


def primer_output_posedit(genname, nm_number, primer, strand, posedit, index):
    # INFO: write generated primers to file
    #   switch forward/reverse based on strand information
    header = f'{genname}, POsition: {posedit}, Primerpaar: {index}\n'
    header_2 = """{0:7}\t{1:<5}\t{2:<6}\t{3:4}\t{4:4}\t{5:35}""".format('', 'Start', 'Length', 'Tm', 'GC%', 'Sequence')
    exon_str = str(posedit).zfill(2)
    temp = (primer['PRIMER_RIGHT_0_TM'] + primer['PRIMER_LEFT_0_TM']) / 2

    if strand == '+':
        primer_pairs = """
Forward\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35}
Reverse\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}
Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}
""".format(**primer)

        info = f"""
{genname}-{exon_str}F;{primer['PRIMER_LEFT_0_SEQUENCE']}
{genname}-{exon_str}R;{primer['PRIMER_RIGHT_0_SEQUENCE']}
{genname}; {round(temp)} °C; {primer['PRIMER_PAIR_0_PRODUCT_SIZE']}bp; {nm_number}

"""
    else:
        primer_pairs = """

Forward\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}
Reverse\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35}
Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}
""".format(**primer)

        info = f"""
{genname}-E{exon_str}F;{primer['PRIMER_RIGHT_0_SEQUENCE']}
{genname}-E{exon_str}R;{primer['PRIMER_LEFT_0_SEQUENCE']}
{genname}; {round(temp)} °C; {primer['PRIMER_PAIR_0_PRODUCT_SIZE']}bp; {nm_number}

"""
    output = [header, header_2, primer_pairs, info]

    return output


def primer_output_genomic(primer, chromosome, start, end, index):
    # INFO: write generated primers to file
    #   switch forward/reverse based on strand information
    # TODO: is this simplifiable? how to handle strand orientation?
    #  Do we need this, if this is put into a database?

    header = f'Chromosome: {chromosome}, Start: {start}, Ende: {end}, Primerpaar: {index}\n'
    header_2 = """{0:7}\t{1:<5}\t{2:<6}\t{3:4}\t{4:4}\t{5:35}""".format('', 'Start', 'Length', 'Tm', 'GC%', 'Sequence')

    primer_pairs = """
Forward\t{PRIMER_LEFT_0[0]:<5}\t{PRIMER_LEFT_0[1]:<6}\t{PRIMER_LEFT_0_TM:2.2f}\t{PRIMER_LEFT_0_GC_PERCENT:2.2f}\t{PRIMER_LEFT_0_SEQUENCE:35}
Reverse\t{PRIMER_RIGHT_0[0]:<5}\t{PRIMER_RIGHT_0[1]:<6}\t{PRIMER_RIGHT_0_TM:2.2f}\t{PRIMER_RIGHT_0_GC_PERCENT:2.2f}\t{PRIMER_RIGHT_0_SEQUENCE:35}
Product Length\t{PRIMER_PAIR_0_PRODUCT_SIZE}
""".format(**primer)

    temp = (primer['PRIMER_RIGHT_0_TM'] + primer['PRIMER_LEFT_0_TM']) / 2
    info = f"""
{chromosome}-{start}F;{primer['PRIMER_LEFT_0_SEQUENCE']}
{chromosome}-{end}R;{primer['PRIMER_RIGHT_0_SEQUENCE']}
{chromosome}-{start}-{end}; {round(temp)} °C; {primer['PRIMER_PAIR_0_PRODUCT_SIZE']}bp;

"""
    output = [header, header_2, primer_pairs, info]

    return output


def write_output_file(outfile, primer_strings):
    print(primer_strings)
    with open(outfile, 'w') as f:
        for item in primer_strings:
            for x in item:
                f.write(x)