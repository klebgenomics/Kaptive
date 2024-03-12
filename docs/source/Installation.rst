====================================
Installation
====================================

Kaptive should work on both Python 2 and 3, but we run/test it on Python 3 and recommend you do the same.

Clone and run
=================

Kaptive is a single Python script, so you can simply clone (or download) from GitHub and run it. Kaptive depends on
`Biopython <http://biopython.org/wiki/Main_Page>`_, so make sure it's installed (either ``pip3 install biopython`` or
read detailed instructions `here <http://biopython.org/DIST/docs/install/Installation.html>`_.

.. code-block:: bash

 git clone https://github.com/klebgenomics/Kaptive.git
 Kaptive/kaptive.py -h

Install with pip
====================

Alternatively, you can install Kaptive using `pip <https://pip.pypa.io/en/stable/>`_.
This will take care of the Biopython requirement (if necessary) and put the ``kaptive.py`` script in your PATH for easy
access. Pip installing will *not* provide the reference databases, so you'll need to download them separately from
`here <https://github.com/katholt/Kaptive/tree/master/reference_database>`_.

.. code-block:: bash

 pip install kaptive
 kaptive.py -h

Other dependencies
======================

Regardless of how you download/install Kaptive, it requires that `BLAST+ <http://www.ncbi.nlm.nih.gov/books/NBK279690/>`_
is available on the command line (specifically the commands ``makeblastdb``, ``blastn`` and ``tblastn``).
BLAST+ can usually be easily installed using a package manager such as `Homebrew <http://brew.sh/>`_ (on Mac) or
`apt-get <https://help.ubuntu.com/community/AptGet/Howto>`_ (on Ubuntu and related Linux distributions). Some later
versions of BLAST+ have been associated with sporadic crashes when running tblastn with multiple threads; to avoid this
problem we recommend running Kaptive with BLAST+ v 2.3.0 or using the ``--threads 1`` option (see below for full command
argument details).
