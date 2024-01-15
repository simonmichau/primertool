Primertool
==========

The premise of this package is to generate primers for PCR/Sanger sequencing for either:

- a specific variant (hgvs_ variant nomenclature), either the whole exon or if not in an exon for the genomic position
- an exon (transcript number and exon )
- all exons of a transcript (transcript number)
- around a genomic position (chromosome and start/stop position)


This tool allows for primers based the human reference genome version hg19 or hg38.


Installation
--------------

The Primertool requires Python 3.6+. All required packages are listed in the requirements file.

Local installation is possible by running

.. code:: bash

    pip install .

in the main directory.

You will also need gcc installed on your system, so check if you already have it by running

.. code:: bash

    gcc --version

If you don't have gcc installed, you can install it by running

.. code:: bash

    sudo apt update
    sudo apt install build-essential

Basic usage examples
---------------------

1. Generating primers for a mutation using PrimerMutation:

.. code:: Python

    import primertool
    x = primertool.PrimerMutation('NM_003165.6:c.1702G>A', 'hg38')
    x.create_primer()

2. Generating primer for a single exon using PrimerExon:

.. code:: Python

    import primertool
    x = primertool.PrimerExon('NM_000451', '6', 'hg38')
    x.create_primer()

3. Generating primers for every exon in a gene using PrimerGen:

.. code:: Python

    import primertool
    x = primertool.PrimerGen('NM_000451', 'hg38')
    x.create_primer()

4. Generating primer for a genomic position using PrimerGenomicPosition:

.. code:: Python

    import primertool
    x = primertool.PrimerGenomicPosition('chr12', 121814175, 121814175, 'hg38')
    x.create_primer()





.. Links
.. _hgvs: https://varnomen.hgvs.org/
