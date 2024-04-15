**************************************
Usage
**************************************

Quickstart (for the impatient)
================================

To type polysaccharide loci from genome assemblies::

   kaptive assembly <database> <path_to_assemblies> > kaptive_results.tsv


This will run ``kaptive assembly`` with the default parameters, and produce a table detailing the best match locus, predicted phenotype, confidence score and detailed typing information for each input genome assembly in the file called results_table.txt (plus other default outputs described below). 


Detailed usage
================

We designed Kaptive 3 to be easier to use on the command-line than previous versions by structuring the program as a
series of sub-commands that follow the general pattern of ``kaptive <mode> <database> <input>``.
There are three modes:

* **assembly**: type assemblies
* **extract**: extract features from Kaptive databases in different formats
* **convert**: convert Kaptive results to different formats

To see the full list of commands and options, run ``kaptive --help``.


kaptive assembly
------------------

Given a Kaptive database and a bacterial genome assembly, ``kaptive assembly`` will perform 3 main tasks:

* Determines the most likely locus type of the genome assembly.
* Reconstructs the biosynthetic gene cluster from the assembly contig sequences.
* Predicts the corresponding serotype/phenotype of the genome assembly.

.. note::
 As of version 3, Kaptive no longer supports allelic (*wzi*, *wzc*) typing.


To perform K-locus typing on a directory of *Klebsiella pneumoniae* assemblies, you would run::

    kaptive assembly kpsc_k assemblies/*.fasta

Here we have told Kaptive to perform typing of assemblies with ``assembly`` and used the database keyword
``kpsc_k`` to specify the *Klebsiella pneumoniae* K-locus database. All other parameters are set to the default.


Database keywords are a handy short-cut for using the databases distributed with Kaptive and located in
the ``reference_databases`` directory. For available database keywords, see :ref:`Database keywords`.


Alternatively, you can specify the full path to your own database.

You may also want to specify the locations and/or filenames of the output files using the following options:: 

     -o , --out            Results table filename/path (default: stdout)
     -f [], --fasta []     Output locus sequence to "{input}_kaptive_results.fna"
                            - Optionally pass output directory (default: cwd)
     -j [], --json []      Output results to JSON lines
                            - Optionally pass file name (default: kaptive_results.json)
     -p [], --plot []      Output locus plots to "{input}_kaptive_results.{fmt}"
                            - Optionally pass output directory (default: cwd)
     --plot-fmt png,svg    Format for locus plots (default: png)
     --no-header           Suppress header line
     --debug               Append debug columns to table output


Advanced options
^^^^^^^^^^^^^^^^^^
Advanced users may wish to customise Kaptive's scoring options (for picking the best match locus), confidence options (for marking matches as 'untypeable', 'good' or 'low' confindence) or database parsing options. We recommend keeping the default options for standard typing using the *Klebsiella* and/or *A. baumanii* databases distributed with Kaptive.


Scoring options::

     --score-metric        Minimap alignment metric to use for scoring (default: AS)
     --weight-metric       Weighting for scoring metric (default: prop_genes_found)
                            - none: No weighting
                            - locus_length: length of the locus
                            - genes_expected: # of genes expected in the locus
                            - genes_found: # of genes found in the locus
                            - prop_genes_found: genes_found / genes_expected
     --min-cov             Minimum gene %coverage to be used for scoring (default: 50.0)

Confidence options::

     --gene-threshold      Minimum translated gene identity required for expected genes (default: database specific)
     --allow-below-threshold     Ignore minimum gene identity threshold
     --percent-expected-genes     Minimum percentage of expected genes required to consider a locus typeable (applies to fragmented loci only, default: 50)
     --max-other-genes     Maximum number of other genes in locus permitted for typeable loci (applies to fragmented loci only, default: 1)


Database options::

     --locus-regex         Pattern to match locus names in db source note
     --type-regex          Pattern to match locus types in db source note
     --filter              Pattern to select loci to include in the database


kaptive convert
----------------

The ``convert`` command allows you to convert the Kaptive results JSON file into a range of useful formats, including:

* :ref:`Tabular` output (tsv)
* :ref:`Locus nucleotide sequence(s) <Fasta>` (fna)
* Locus gene nucleotide sequences (ffn)
* Locus gene amino acid sequences (faa)
* Locus :ref:`plot <Plot>`

This means if you didn't want to or forgot to output these files during the initial run, we've got you covered!

Simply run ``kaptive convert <JSON file> <format>`` and the file will be output to the current directory.


API
------
Whilst Kaptive isn't designed to be a full API, it is possible to use it as a module in your own Python scripts.
For typing assemblies, you can use the ``kaptive.assembly.typing_pipeline`` function, which takes an assembly path and a
``kaptive.database.Database`` object as input and returns a ``kaptive.typing.TypingResult`` object.

.. code-block:: python

    from kaptive.database import Database, get_database
    from kaptive.assembly import typing_pipeline
    from pathlib import Path

    db = Database.from_genbank(database_path)
    results = [typing_pipeline(assembly, db, threads=8) for assembly in Path('assemblies').glob('*.fna')]

For example, if you wanted to perform K and O locus typing on a single assembly, you could do the following::

    k_db, o_db = get_database('kpsc_k'), get_database('kpsc_o')
    k_db, o_db = Database.from_genbank(k_db), Database.from_genbank(o_db)
    k_results, o_results = typing_pipeline(a, k_db, threads=8), typing_pipeline(a, o_db, threads=8)
    print(k_results.as_table(), o_results.as_table())

