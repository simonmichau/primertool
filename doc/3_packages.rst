Needed packages and their documentation
=======================================



Primer3
-------

Primer3 exists both as online_ version and as python package_. An additional manual_  exists, that contains explanations for all variables.

This package generates primers with a given genomic sequence and a dictionary containing all needed variables.


genomepy
--------

genomepy_ allows easy installation and use of genomes in python. Here it is used to download the needed version of the human genome (hg19 or hg38) and extract the sequences in which the primers are generated.


hgvs.parser
-----------

The hgvs_ parser is a package to handle biological sequence variants based on the HGVS nomenclature. Here it is used to parse a variant in


.. Links
.. _online: https://primer3.ut.ee/
.. _package: https://libnano.github.io/primer3-py/index.html
.. _manual:  http://htmlpreview.github.io/?https://github.com/libnano/primer3-py/master/primer3/src/libprimer3/primer3_manual.htm
.. _genomepy: https://github.com/vanheeringen-lab/genomepy
.. _hgvs: https://hgvs.readthedocs.io/en/stable/index.html