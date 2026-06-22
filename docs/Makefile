# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
ifdef OS
   SPHINXBUILD   = python -msphinx
else
   SPHINXBUILD   = python3 -msphinx
endif
# SPHINXBUILD   = python3 -msphinx
SPHINXPROJ    = PoreMS
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@if [ "$@" = "html" ] && [ -d v0.3 ]; then \
		cp -r v0.3 $(BUILDDIR)/html/v0.3; \
	fi
