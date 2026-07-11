# Pore Generator for Molecular Simulations

:::{figure} /pics/pore.svg
:align: center
:width: 60%
:::

This python package generates pore structures to be used in molecular simulations.
For an overview of the programs operating principle, check out the corresponding
publication **Kraus et al.** (doi:[10.1080/08927022.2020.1871478](https://doi.org/10.1080/08927022.2020.1871478)).
Additionally, tutorials for generating [molecules](molecule.md) and [pores](pore.md) are provided.

Check out an exemplary [workflow](workflow.md) for using the PoreMS package to create a pore system
and run molecular dynamics simulation using [Gromacs](http://www.gromacs.org/).

To see the code, report a bug or contribute, please visit the
[GitHub repository](https://github.com/PoreMS/PoreMS).

A browser-based interface for PoreMS is available via the
[PoreMS App](https://github.com/PoreMS/PoreMS-App) — a step-by-step wizard for
configuring pore geometry, attaching surface molecules, and downloading a
GROMACS-ready structure without writing Python scripts.

## API Reference

Full API documentation is auto-generated from source: [API Reference](autoapi/index)

:::{toctree}
:hidden:
:maxdepth: 2

molecule
pore
workflow
shape_examples
autoapi/index
:::
