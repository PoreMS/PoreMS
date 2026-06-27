# v0.5.0

### Performance
* `Molecule`: positions stored as a NumPy array; bulk operations avoid per-atom Python loops
* `geometry.py`: scalar functions (`cross_product`, `length`, `dot_product`, `unit`, `angle`, `rotate`) rewritten with `math` module and manual formulas — eliminates NumPy Python-layer overhead (`normalize_axis_tuple`, `moveaxis`) for 3D vectors called ~16k times per pore; `rotate()` uses a fast list path for single vectors and NumPy `einsum` for arrays
* `Dice.find_parallel()`: replaced bound-method `apply_async(self.find_bond, ...)` with module-level `_find_bond_task()` and `starmap()` to avoid pickling the full `Molecule` object per worker; uses `fork` context on single-threaded Unix processes (zero-copy via copy-on-write) with `spawn` fallback on Windows and in test environments — **~2× geometric mean speedup vs 0.3.0 (up to 5× for small pores)**
* `Dice.find_bond()`: vectorized pairwise distance check with NumPy broadcasting; pre-computed `mol_box` array moved to constructor
* `Dice.neighbor()`: replaced fragile `pop(13)` magic index with explicit `n != cube_id` filter
* `Molecule.overlap()`: replaced O(n²) loop with SciPy `KDTree` (Chebyshev metric)
* `Store.gro()`, `Store.pdb()`: replaced per-atom string concatenation with list + `writelines`
* `Store.lmp()`: replaced `list.index()` per atom with pre-built dict for O(1) atom-type lookup
* `Pore.sites()`: `_oxygen_ex` converted to set for O(1) membership test
* `Pore.attach()`: pre-built `{si_id: idx}` dict replaces `list.index()` in placement loop
* `Pore.prepare()`: Counter over bound O atoms computed once per while-iteration, not twice

### Logic fixes
* `Pore.amorph()`: bond check now correctly evaluates the trial (displaced) position; atoms are tentatively moved, checked, and reverted on rejection — previously the pre-displacement position was always checked, making every displacement pass
* `PoreAmorphCylinder.attach_special()`: removed spurious `self._normal_in` positional argument that did not exist and would raise `AttributeError` at runtime; argument list now matches `Pore.attach()` signature
* `system.py:table()`: bare `except: pass` replaced with `except AttributeError`; magic literal `20` replaced with `_UNASSIGNED_KEY`; `sites_attach_mol[i]["SL"]`/`["SLG"]` accesses use `.get()` with default 0 to guard against incomplete site dictionaries
* `system.py:attach()`: shape index parsing changed from `int(shape[-1])` (broke for indices ≥ 10) to `int(shape.split("_")[1])`

### New features
* `PoreCone`: convenience class for conical pores (linearly varying diameter)
* `attach_special()` moved to `PoreKit` base class — removed identical duplicate methods from `PoreCylinder`, `PoreSlit`, and `PoreAmorphCylinder`
* `Hourglass` shape class added to `shape.py` with correct `is_in()` and `normal()` implementations

### CI / tooling
* GitHub Actions: added ruff linting workflow (`lint.yml`)
* GitHub Actions: added pip-audit security scan workflow (`security.yml`)
* GitHub Actions: publishing handled via existing `python-publish.yml` (manually triggered)
* Dependabot configuration added for automated dependency updates
* Python version matrix updated: 3.10–3.13; `python_requires` bumped to `>=3.10`

### Documentation
* RST source files migrated to MyST Markdown
* Sphinx theme updated to furo; API docs via sphinx-autoapi
* Copyright year updated to 2026

### Administrative
* `setup.py`: version 0.5.0, `python_requires='>=3.10'`, author email updated
* README: fixed `docsrc/` → `docs/` image paths, updated Python version, added testing commands


# v0.4.0
* New version due to a change of GitHub organisation.
* Option in the 'attach' function to force molecules to be added only to single binding sites.


# v0.3.0
* Reworked generation methodology for easier shape combination and pattern usage
* Added new alpha-cristobalite pattern
* Generalization of connectivity search
* New Lennard Jones parameter (Coasne, B.; Fourkas, J. T. from 10.1021/jp203831q)

# v0.2.5
* Add properties output in yaml format

# v0.2.4
* Add PoreAmorphCylinder class for amorphous generation
* Rewrite exterior surface site determination - now works independet of exterior surface
* Distance search now takes in a range instead of a distance and error
* Bugfixes and quality of life improvements

# v0.2.3
* Bugfixes and quality of life improvements

# v0.2.2
* System Class: Move duplicate methods of child classes to parent class
* System Class: Remove other output tables except full table
* Slit-Pore: Possibility to add reservoir
* Slit-Pore: Possibility to add siloxan bridges
* Cylinder-Pore: Possibility to remove reservoir
* Capsule-Pore: Possibility to remove reservoir
* Capsule-Pore: Possibility to add siloxan bridges
* Store Class: Add LAMMPS structure storage function

# v0.2.1
* Bugfixes and Improvements
