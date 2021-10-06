Generating Primers
==================

General requirements for primers
---------------------------------

Since the primers are used for sanger sequencing, these settings are important to create usable primers. Of course the
standard settings can be changed if needed.

- maximum insert size for a primer: 700 bp
- minimum insert size for a primer: 200 bp
- distance of primers from exon borders: 40 bp


These are settings used for primer3:

- optimal primer size: 20 bp 
- minimum primer size: 20 bp
- maximum primer size: 22 bp
- optimal temperature: 60°C
- minimum temperature: 58°C
- maximum temperature: 62°C
- maximum poly x (maximum allowable length of a mononucleotide repeat): 5
- gcc clamp (Require the specified number of consecutive Gs and Cs at the 3' end of both the left and right primer): 1


PrimerMutation
---------------

A variant in HGVS nomenclature is used to generate primers for the specific position. If the variant is located in an exon, primers are generated for the whole exon.

1. Check if given reference genome is downloaded and if not, download the specific file.
2. Check the given mutation with mutalyzer and convert from coding into the genomic position.
3. Extract transcript number (NM) and transcript specific information from the NCBI Refseq database.



PrimerGenomic
-------------

This is the easiest version to generate primers, since it only requires a start and end position, chromosome and the name of the reference genome.

1. Check if given reference genome is downloaded and if not, download the specific file.
2. Are start and stop position in the min/max insert size.