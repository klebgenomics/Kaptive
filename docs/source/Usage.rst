**************************************
Usage
**************************************

Quickstart (for the impatient)
================================

To type polysaccharide loci from genome assemblies::

   kaptive assembly <database> <path_to_assemblies> > kaptive_results.tsv


This will run ``kaptive assembly`` with the default parameters, and produce a table detailing the best match locus,
predicted phenotype, confidence score and detailed typing information for each input genome assembly in the file
called ``results_table.txt``.


Detailed usage
================

We designed Kaptive 3 to be easier to use on the command-line than previous versions by structuring the program as a
series of sub-commands that follow the general pattern of ``kaptive <mode> <database> <input>``.
There are three modes:

* **assembly**: :ref:`type assemblies <kaptive-assembly>`
* **extract**: :ref:`extract <kaptive-extract>` features from Kaptive databases in different formats
* **convert**: :ref:`convert <kaptive-convert>` Kaptive results to different formats

.. note::
 To see the full list of commands and options, run ``kaptive --help``.

.. _kaptive-assembly:

kaptive assembly
------------------

Given a Kaptive database and a bacterial genome assembly, ``kaptive assembly`` will perform 3 main tasks:

* Determines the most likely locus type of the genome assembly.
* Reconstructs the biosynthetic gene cluster from the assembly contig sequences.
* Predicts the corresponding serotype/phenotype of the genome assembly.

.. note::
 As of version 3, Kaptive no longer supports allelic (*wzi*, *wzc*) typing.


To perform K locus typing on a directory of *Klebsiella pneumoniae* assemblies, you would run::

    kaptive assembly kpsc_k assemblies/*.fasta > kaptive_results.tsv

Here we have told Kaptive to perform typing of assemblies with ``assembly`` and used the database keyword
``kpsc_k`` to specify the *Klebsiella pneumoniae* K locus database. All other parameters are set to the default.


Database keywords are a handy short-cut for using the databases distributed with Kaptive and located in
the ``reference_databases`` directory. See available `keywords <Database-keywords>`.


Alternatively, you can specify the full path to your own database.

You may also want to specify the locations and/or filenames of the output files using the following options::

    -o file, --out file   Output file to write/append tabular results to (default: stdout)
    -f [dir], --fasta [dir]
                        Turn on fasta output, defaulting "{input}_kaptive_results.fna"
                         - Optionally choose the output directory (default: cwd)
    -j [file], --json [file]
                        Turn on JSON lines output
                         - Optionally choose file (can be existing) (default: kaptive_results.json)
    -p [dir], --plot [dir]
                        Turn on plot output, defaulting to "{input}_kaptive_results.{fmt}"
                         - Optionally choose the output directory (default: cwd)
    --plot-fmt png,svg    Format for locus plots (default: png)
    --no-header           Suppress header line
    --debug               Append debug columns to table output

Example::

    kaptive assembly kpsc_k assemblies/*.fasta -o kaptive_results.tsv -f -j -p

This will output a tabular file called ``kaptive_results.tsv``, a fasta file for each assembly called
``<assembly_name>_kaptive_results.fna``, a JSON lines file called ``kaptive_results.json`` and a plot for each assembly
called ``<assembly_name>_kaptive_results.png``.

Advanced options
^^^^^^^^^^^^^^^^^^
Advanced users may wish to customise Kaptive's scoring options (for picking the best match locus), confidence options
(for marking matches as 'Typeable' or 'Untypeable') or database parsing options. We recommend keeping the default
options for standard typing using the *Klebsiella* and/or *A. baumanii* databases distributed with Kaptive.


:ref:`Scoring options <Scoring-algorithm>`::

    --score-metric        Alignment metric to use for scoring (default: AS)
    --weight-metric       Weighting for scoring metric (default: prop_genes_found)
                         - none: No weighting
                         - locus_length: length of the locus
                         - genes_expected: # of genes expected in the locus
                         - genes_found: # of genes found in the locus
                         - prop_genes_found: genes_found / genes_expected
    --min-cov             Minimum gene %coverage to be used for scoring (default: 50.0)

:ref:`Confidence options <Confidence-score>`::

    --gene-threshold      Species-level locus gene identity threshold (default: database specific)
    --max-other-genes     Typeable if <= other genes (default: 1)
    --percent-expected-genes
                        Typeable if >= % expected genes (default: 50)
    --allow-below-threshold
                        Typeable if any genes are below threshold

See database options :ref:`here <Database-options>` and other options::

    -V, --verbose         Print debug messages to stderr
    -v , --version        Show version number and exit
    -h , --help           Show this help message and exit
    -t , --threads        Number of threads for alignment (default: maximum available CPUs)

.. _kaptive-convert:

kaptive convert
----------------
The ``convert`` command allows you to convert the Kaptive results JSON file into a range of useful formats, including:

* **tsv**: :ref:`Tabular` output (tsv)
* **json**: JSON lines format (same as input but optionally filtered)
* **loci**: :ref:`Locus nucleotide sequence(s) <Fasta>` in fasta (fna) format
* **genes**: Locus gene nucleotide sequences in fasta (ffn) format
* **proteins**: Locus gene amino acid sequences in fasta (faa) format
* **png**: Locus :ref:`plots <Plot>` as PNG
* **svg**: Locus :ref:`plots <Plot>` as SVG

This means if you didn't want to or forgot to output these files during the initial run, we've got you covered!

Simply run ``kaptive convert <JSON file> <format>`` and the file will be output to the current directory.

For example, to convert the JSON file to a tabular format, you would run::

    kaptive convert kaptive_results.json tsv > kaptive_results.tsv

OR if you have multiple JSON files::

        cat *.json | kaptive convert - tsv > kaptive_results.tsv

OR if I want to convert the results to a protein fasta of all the locus genes::

    kaptive convert kaptive_results.json proteins > proteins.faa

OR if I want to do the same as above but generate a separate file for each sample::

    kaptive convert kaptive_results.json proteins -d proteins

Where the ``-d`` option specifies the output directory for the converted results and each sample will have its own file
with the name ``<sample_name>_kaptive_results.faa``.

.. note::
 Plots will always be written to files, even if the output is set to stdout, and one file will be written per sample.

Inputs::

    db path/keyword       Kaptive database path or keyword
    json                  Kaptive JSON lines file or - for stdin
    format                Output format
                             - json: JSON lines format (same as input but optionally filtered)
                             - tsv: Tab-separated values (results table)
                             - loci: Locus nucleotide sequence in fasta format
                             - proteins: Locus proteins in fasta format
                             - genes: Locus genes in fasta format
                             - png: Locus plot in PNG format
                             - svg: Locus plot in SVG format

Filter options::

    -r , --regex          Python regular-expression to select JSON lines (default: All)
    -l  [ ...], --loci  [ ...]
                        Space-separated list to filter locus names (default: All)
    -s  [ ...], --samples  [ ...]
                        Space-separated list to filter sample names (default: All)

Output options::

    -o , --out            Output file to write/append results to (default: stdout)
                         - Note: Only for text formats, figures will be written to files
    -d , --outdir         Output directory for converted results
                         - Note: This forces the output to be written to files
                                 If used with locus, proteins or genes, one file will be written per sample

.. note::
 Filters take precedence in descending order


See database options :ref:`here <Database-options>` and other options::

    -V, --verbose         Print debug messages to stderr
    -v , --version        Show version number and exit
    -h , --help           Show this help message and exit



.. _api:

API
------
Whilst Kaptive isn't designed to be a full API, it is possible to use it as a module in your own Python scripts.
For typing assemblies, you can use the ``kaptive.assembly.typing_pipeline`` function, which takes an assembly path and a
``kaptive.database.Database`` object as input and returns a ``kaptive.typing.TypingResult`` object.

.. code-block:: python

    from kaptive.database import Database, get_database
    from kaptive.assembly import typing_pipeline
    from pathlib import Path

    db = Database.from_genbank(get_database('kp_k'))
    results = [typing_pipeline(assembly, db, threads=8) for assembly in Path('assemblies').glob('*.fna')]

For example, if you wanted to perform K and O locus typing on a single assembly, you could do the following:

.. code-block:: python

    k_results = typing_pipeline(a, Database.from_genbank(get_database('kp_k')), threads=8)
    o_results = typing_pipeline(a, Database.from_genbank(get_database('kp_o')), threads=8)
    print(k_results.as_table(), o_results.as_table())

.. note::
 By default the ``typing_pipeline`` function runs ``minimap2`` on a single thread, which is recommended for running
 multiple assemblies in parallel.

If you have lots of CPUs and know how many assemblies you have, it may be faster to use multiprocessing to type multiple
assemblies at once. Here's an example of how you could do that:

.. code-block:: python

    from kaptive.database import Database, get_database
    from kaptive.assembly import typing_pipeline
    from pathlib import Path
    from concurrent.futures import ProcessPoolExecutor as Executor

    db = Database.from_genbank(get_database('kp_k'))
    with Executor(max_workers=8) as executor:
        results = list(executor.map(lambda a: typing_pipeline(a, db), Path('assemblies').glob('*.fna')))

The ``TypingResult`` object is designed to only include information about the locus, and not the target
assembly, partially in an attempt to reduce its memory footprint. This means you can return the results from
multiples assemblies in parallel and then access the information you need from them.

.. note::
 This wasn't included in the Kaptive 3 CLI as writing the output in parallel is more complex, but we are open
 to suggestions if this dramatically improves performance!


