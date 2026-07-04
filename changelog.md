# v1.0.0

### Performance

Benchmarked on a 16-core macOS machine, v1.0.0 vs PyPI v0.3.0
(`python tests/bench_compare.py`):

| Benchmark | v1.0.0 | v0.3.0 | Speedup |
|---|---|---|---|
| BetaCristobalit pattern 4³ | 0.033 s | 0.046 s | **1.4×** |
| BetaCristobalit pattern 8³ | 0.137 s | 0.160 s | **1.2×** |
| PoreCylinder 4 nm | 0.213 s | 1.131 s | **5.3×** |
| PoreCylinder 6 nm | 0.686 s | 1.441 s | **2.1×** |
| PoreCylinder 8 nm | 1.398 s | 1.967 s | **1.4×** |
| PoreSlit 6 nm | 0.676 s | 1.458 s | **2.2×** |
| **Geometric mean** | | | **~2×** |

Key changes driving the speedup:

* `geometry.py`: scalar functions (`cross_product`, `length`, `dot_product`, `unit`, `angle`, `rotate`) rewritten with `math` module and manual formulas — eliminates NumPy Python-layer overhead (`normalize_axis_tuple`, `moveaxis`) for 3D vectors called ~16k times per pore; `rotate()` uses a fast list path for single vectors and NumPy `einsum` for arrays
* `Dice.find_parallel()`: replaced bound-method `apply_async(self.find_bond, ...)` with module-level `_find_bond_task()` and `starmap()` to avoid pickling the full `Molecule` object per worker; uses `fork` context on single-threaded Unix processes (zero-copy via copy-on-write) with `spawn` fallback on Windows and in test environments
* `Dice.find_bond()`: vectorized pairwise distance check with NumPy broadcasting; pre-computed `mol_box` array moved to constructor
* `Dice.neighbor()`: replaced fragile `pop(13)` magic index with explicit `n != cube_id` filter
* `Molecule`: positions stored as a NumPy array; bulk operations avoid per-atom Python loops
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
* `PoreKit.table()`: new `fmt` parameter — `fmt="plain"` returns a Unicode-box styled string (μmol/m², nm³ units; section headers; aligned columns), `fmt="latex"` returns a `longtable` environment (booktabs rules, `\quad`-indented sub-rows, LaTeX math units) suitable for direct inclusion in papers

### CI / tooling
* GitHub Actions: added ruff linting workflow (`lint.yml`)
* GitHub Actions: added pip-audit security scan workflow (`security.yml`)
* GitHub Actions: publishing handled via existing `python-publish.yml` (manually triggered)
* Dependabot configuration added for automated dependency updates
* Python version matrix updated: 3.12–3.13; `python_requires` bumped to `>=3.12`
* Migrated from `setup.py` + `MANIFEST.in` + `pytest.ini` to a single `pyproject.toml` (PEP 517/621); `[tool.ruff]` config added
* CI migrated from `pip` to `uv` (`astral-sh/setup-uv@v5`); `requirements.txt` removed (deps resolved via `pyproject.toml`)

### Documentation
* RST source files migrated to MyST Markdown
* Sphinx theme updated to furo; API docs via sphinx-autoapi
* Copyright year updated to 2026

### Administrative
* `pyproject.toml`: version 1.0.0, `requires-python = ">=3.12"`, author email updated; `[project.optional-dependencies]` dev group added (`pytest`, `pytest-cov`)
* README: fixed `docsrc/` → `docs/` image paths, updated Python version, added testing and benchmarking commands
* `tests/bench_compare.py`: standalone script comparing generation speed of the local build vs any PyPI release
* `tests/test_bench.py`: benchmark time limits tightened to reflect real performance; `POREMS_BENCH_SCALE` env var allows scaling limits on slow machines


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
