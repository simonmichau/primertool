Examples


NM_014946.3:c.131C>T: no primers can be found using PrimerExon, but x = primertool.PrimerGenomicPosition('chr2', 32063820, 32064250, 'hg38')
 works

NM_001353345:c.1960G>A: give error of incorrect version number using name checker (version number is 2), but conversion is not possible. Mutalyzer web version gives following error: The accession number NM_001353345 could not be found in our database (or is not a transcript). Example of inconsistency between name checker and position converter.


NM_003165:c.1702+1G>A: primertool.exceptions.PrimertoolInputError: ('There was a problem with the input. ', 'WNOVER', 'No version number is given, using NM_003165.6. Please use this number to reduce downloading overhead.')
And +/- variants which are very close to exon borders should still be first tried to generate primers for the exon.


