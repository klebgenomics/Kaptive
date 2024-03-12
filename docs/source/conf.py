# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Kaptive'
copyright = '2024, Tom Stanton, Ryan Wick, Kathryn Holt, Kelly Wyres'
author = 'Tom Stanton'
release = '2.0.9'
version = '2.0.9'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autosectionlabel']
autosectionlabel_prefix_document = True
templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_theme_options = {
    "show_toc_level": 2
}
html_static_path = ['_static']
html_logo = 'https://github.com/klebgenomics/Kaptive/blob/master/extras/kaptive_logo.png?raw=true'
html_favicon = 'https://github.com/klebgenomics/Kaptive-Web/blob/master/static/images/favicon.png?raw=true'
