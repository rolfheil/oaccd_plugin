#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

# -*- coding: utf-8 -*-
#
# Psithon documentation build configuration file, created by
# sphinx-quickstart on Sun Feb 12 04:25:25 2012.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys, os

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('/home/rolf/psi4/objdir/stage//usr/local/psi4/lib//psi4/driver'))  # for qcdb for db
sys.path.insert(0, os.path.abspath('/home/rolf/psi4/objdir/stage//usr/local/psi4/lib/'))
#sys.path.insert(0, os.path.abspath('/home/rolf/psi4/objdir/stage//usr/local/psi4/lib/psi4'))
#sys.path.insert(0, os.path.abspath('/home/rolf/psi4/psi4/driver/'))
sys.path.insert(0, os.path.abspath('/home/rolf/psi4/psi4/share/psi4/databases/'))
sys.path.insert(0, os.path.abspath('/home/rolf/psi4/plugins/'))
sys.path.insert(0, os.path.abspath('/home/rolf/psi4/doc/sphinxman/source'))

# Import Sphinx themes
import cloud_sptheme as csp
import psi4doc as psp


# -- General configuration -----------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
# (1.3 needed by cloud. 1.3 works if unfix an error fixed for 1.4)
needs_sphinx = '1.4'

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
# 'psi4_sptheme.ext.autodoc_sections',
extensions = [
    # from Sphinx
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.extlinks',
    'sphinx.ext.graphviz',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    # from Jupyter
      # nbsphinx
    # from Cloud
    'cloud_sptheme.ext.index_styling',
    'cloud_sptheme.ext.escaped_samp_literals',
    # from Astropy
    'astropy_helpers.sphinx.ext.automodapi',
    'astropy_helpers.sphinx.ext.automodsumm',
    # from Psi4
    'psi4doc.ext.psidomain',
    'psi4doc.ext.relbar_toc',
]

autosummary_generate = True
autodoc_default_flags = ['members',
                         'undoc-members',
                         'inherited-members',  # disabled because there's a bug in sphinx
                         'show-inheritance',
                        ]
automodapi_toctreedirnm = 'source/api'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# Suppress warnings (sphinx 1.4.2)
suppress_warnings = ['image.nonlocal_uri']

# General information about the project.
project = u'Psi4'  # u'P\uA731\u026A4'
copyright = u'2017, The {0} Project'.format(project)

def get_version_info():
    import psi4
    version = psi4.__version__
    githash = psi4.version_formatter('{githash}')
    return version, githash

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
# The full version, including alpha/beta/rc tags.
version, release = get_version_info()

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
#exclude_patterns = []
exclude_patterns = ['template_index.rst',          'template_appendices.rst',
                    'autodoc_index_html.rst',      'autodoc_appendices_html.rst', 
                    'autodoc_index_htmluser.rst',  'autodoc_appendices_htmluser.rst', 
                    'autodoc_index_htmlprog.rst',  'autodoc_appendices_htmlprog.rst', 
                    'autodoc_index_latexuser.rst', 'autodoc_appendices_latexuser.rst', 
                    'autodoc_index_latexprog.rst', 'autodoc_appendices_latexprog.rst',
                    'abbr_accents.rst',
                    'autodoc_abbr_options_c.rst',
                    'autodoc_abbr_options_plugins.rst']


# The reST default role (used for this markup: `text`) to use for all documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False
show_authors = True

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []


# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'psi4doc'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {"roottarget": "index" }

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = [psp.get_theme_dir(), csp.get_theme_dir()]

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = ''

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "psi4square.png"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = 'favicon-psi4.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%A, %d %B %Y %I:%M%p'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {'**': ['localtoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html']}
html_sidebars = {'**': ['searchbox.html', 'localtoc.html'],
                 'index': ['searchbox.html', 'impressum.html']}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
#html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'Psi4doc'


# -- Options for LaTeX output --------------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
'preamble': '\\usepackage{amsfonts}\\usepackage{amssymb}',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
  ('index', 'Psi4.tex', u'Psi4 Documentation',
   u'Psi4 Project', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True

nbsphinx_allow_errors = True
nbsphinx_execute = 'never'
nbsphinx_timeout = 180

# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'psithon', u'Psithon Documentation',
     [u'Psi4 Project'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output ------------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'Psithon', u'Psithon Documentation',
   u'Psi4 Project', 'Psithon', 'One line description of project.',
   'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'

# Abbreviations
rst_epilog = """

.. color: #273896;">
.. |PSIfour| raw:: html

    <span style="font-family: Optima, sans-serif; text-transform: none;">P<span style="font-size: 82%;">SI</span>4</span>

.. |PSIfours| replace:: |PSIfour|\ 's
.. |psirc| replace:: :file:`~/.psi4rc`
.. |dl| replace:: :math:`\Rightarrow`
.. |dr| replace:: :math:`\Leftarrow`
.. |kcalpermol| replace:: kcal mol\ :sup:`-1`
.. |Angstrom| replace:: |AA|\ ngstr\ |o_dots|\ m
.. |MollerPlesset| replace:: M\ |o_slash|\ ller--Plesset
.. |--| unicode:: U+02013 .. en dash
   :trim:
.. |w--w| unicode:: U+02013 .. en dash
.. |---| unicode:: U+02014 .. em dash
   :trim:
.. |w---w| unicode:: U+02014 .. em dash
.. include:: /abbr_accents.rst
"""

#.. |PSIfourIM| image:: /psi4_blue.png
#   :height: 20px

#.. |scs| unicode:: U+A731
#.. |sci| unicode:: U+026A
#.. |PSIfourSC| replace:: P\ |scs|\ |sci|\ 4

#.. role:: sc
#.. raw:: html
#
#  <style type="text/css"><!--
#   .green {color: red;}
#   .sc {font-variant: small-caps;}
#   --></style>
#
#.. |PSIfour| replace:: :sc:`Psi4`

rst_prolog = """
.. highlight:: python
   :linenothreshold: 1
"""
# Logo at top of page
#.. image:: /psi4banner.png
#   :width: 100 %
#   :alt: Psi4 Project Logo

extlinks = {'source':    ('https://github.com/psi4/psi4/blob/master/%s', 'psi4/'),
            'srcsample': ('https://github.com/psi4/psi4/blob/master/samples/%s/input.dat', ''),
            'srcbasis':  ('https://github.com/psi4/psi4/blob/master/psi4/share/psi4/basis/%s.gbs', ''),
            'srcplugin': ('https://github.com/psi4/psi4/blob/master/plugins/%s', ''),
            'srcefpfrag':('https://github.com/ilyak/libefp/blob/master/fraglib/%s.efp', ''),
            'srcdb':     ('https://github.com/psi4/psi4/blob/master/psi4/share/psi4/databases/%s.py', '') }

