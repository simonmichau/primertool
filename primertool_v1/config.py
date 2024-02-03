#!/usr/bin/env python
# encoding: utf-8
# AUTHOR: Simon Michau (smichau@ukaachen.de)
import logging

KUERZEL = None  # Kuerzel of the user for the order form
GENOME = None  # Genome assembly (hg19, hg38)
VARIANT = None  # HGVS variant [only used for PrimerMutation]
NM_NUMBER = None  # NM number [used for PrimerExon and PrimerGen]
EXON_NUMBER = None  # Exon number [only used for PrimerExon]
CHROMOSOME = None  # Chromosome [only used for PrimerGenomicPosition]
START = None  # Start position [only used for PrimerGenomicPosition]
END = None  # End position [only used for PrimerGenomicPosition]

# Logging config
LOGGING_LEVEL = logging.INFO  # Logging level


class Config:
    """ Class to manage global variables and configs for the primertool package """

    def __init__(self, **kwargs):
        self.KUERZEL = None  # Kuerzel of the user for the order form
        self.GENOME = None  # Genome assembly (hg19, hg38)
        self.VARIANT = None  # HGVS variant [only used for PrimerMutation]
        self.NM_NUMBER = None  # NM number [used for PrimerExon and PrimerGen]
        self.EXON_NUMBER = None  # Exon number [only used for PrimerExon]
        self.CHROMOSOME = None  # Chromosome [only used for PrimerGenomicPosition]
        self.START = None  # Start position [only used for PrimerGenomicPosition]
        self.END = None  # End position [only used for PrimerGenomicPosition]
