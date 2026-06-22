import sys
import os

sys.path.insert(0, os.path.abspath('..'))

# -- General configuration ------------------------------------------------

extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'myst_parser',
    'autoapi.extension',
]

# autoapi: scan the porems package and generate API pages automatically
autoapi_dirs = ['../porems']
autoapi_type = 'python'
autoapi_options = [
    'members',
    'undoc-members',
    'show-inheritance',
    'show-module-summary',
]
autoapi_root = 'autoapi'
autoapi_keep_files = True

source_suffix = {'.md': 'myst', '.rst': 'restructuredtext'}

master_doc = 'index'

project = 'PoreMS'
copyright = '2026, Hamzeh Kraus'
author = 'Hamzeh Kraus'

version = '0.5'
release = '0.5.0'

language = 'en'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

todo_include_todos = False


# -- HTML output ----------------------------------------------------------

html_theme = 'furo'

html_theme_options = {
    'light_logo': 'logo_text.svg',
    'dark_logo': 'logo_text.svg',
}

html_favicon = 'favicon.ico'

# _static holds favicon, style, etc.; pics holds logo SVGs (served at _static/ root after build)
html_static_path = ['_static', 'pics']

htmlhelp_basename = 'PoreMSdoc'


# -- LaTeX output ---------------------------------------------------------

latex_documents = [
    (master_doc, 'PoreMS.tex', 'PoreMS Documentation', 'Hamzeh Kraus', 'manual'),
]

man_pages = [
    (master_doc, 'porems', 'PoreMS Documentation', [author], 1)
]

texinfo_documents = [
    (master_doc, 'PoreMS', 'PoreMS Documentation',
     author, 'PoreMS', 'Pore Generator for Molecular Simulations.',
     'Miscellaneous'),
]
