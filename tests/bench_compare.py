"""
Generation speed benchmark: local optimized build vs PyPI porems.

Creates a temporary virtualenv with the latest PyPI release, runs identical
pore-generation benchmarks in both environments, then prints a side-by-side
comparison table.

Usage
-----
    python tests/bench_compare.py              # compare local vs latest PyPI
    python tests/bench_compare.py --pypi 0.3.0 # pin a specific PyPI version
    python tests/bench_compare.py --quick      # fewer / smaller cases
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
import venv

# ---------------------------------------------------------------------------
# Benchmark suite
# This code is sent verbatim to both Python interpreters via -c.
# Keep it dependency-free beyond porems itself.
# ---------------------------------------------------------------------------

_BENCH_CASES_FULL = [
    ("pattern_4³",    "pms.BetaCristobalit().generate([4,4,4], 'z')"),
    ("pattern_8³",    "pms.BetaCristobalit().generate([8,8,8], 'z')"),
    ("cylinder_4nm",  "pms.PoreCylinder([4,4,4], 2.0, res=0).finalize()"),
    ("cylinder_6nm",  "pms.PoreCylinder([6,6,6], 3.0, res=0).finalize()"),
    ("cylinder_8nm",  "pms.PoreCylinder([8,8,8], 4.0, res=0).finalize()"),
    ("slit_6nm",      "pms.PoreSlit([6,6,6], 3.0, res=0).finalize()"),
    ("capsule",       "pms.PoreCapsule([6,6,8], 3.0, 4.0, res=0).finalize()"),
]

_BENCH_CASES_QUICK = [
    ("pattern_4³",   "pms.BetaCristobalit().generate([4,4,4], 'z')"),
    ("cylinder_4nm", "pms.PoreCylinder([4,4,4], 2.0, res=0).finalize()"),
    ("cylinder_6nm", "pms.PoreCylinder([6,6,6], 3.0, res=0).finalize()"),
    ("slit_6nm",     "pms.PoreSlit([6,6,6], 3.0, res=0).finalize()"),
]


def _make_bench_script(cases):
    """Return a Python -c script that runs cases and prints JSON results."""
    lines = [
        "import json, time, sys, warnings",
        "warnings.filterwarnings('ignore')",
        "import porems as pms",
        "print('version:' + getattr(pms, '__version__', '?'), file=sys.stderr)",
        # Warm up: trigger any lazy imports (scipy KDTree, etc.) before timing.
        # Run a small overlap check so compiled extension modules are fully loaded.
        "try:",
        "    _m = pms.BetaCristobalit().generate([2,2,2], 'z')",
        "    _m.overlap([0,0,0], 0.1)",
        "except Exception:",
        "    pass",
        "results = {}",
    ]
    for label, expr in cases:
        lines += [
            f"t0 = time.perf_counter()",
            f"try:",
            f"    {expr}",
            f"    results[{label!r}] = round(time.perf_counter() - t0, 3)",
            f"except Exception as e:",
            f"    results[{label!r}] = f'ERROR: {{e}}'",
        ]
    lines.append("print(json.dumps(results))")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Venv helpers
# ---------------------------------------------------------------------------

def _create_pypi_venv(version, tmpdir):
    """Create a venv and pip-install porems==version from PyPI."""
    venv_dir = os.path.join(tmpdir, "pypi_venv")
    print(f"  Creating venv at {venv_dir} ...", flush=True)
    venv.create(venv_dir, with_pip=True, clear=True)

    pip = os.path.join(venv_dir, "bin", "pip")
    spec = f"porems=={version}"
    print(f"  Installing {spec} from PyPI ...", flush=True)
    # Run pip from tmpdir so it cannot find the local porems source tree.
    # install scipy separately — older porems releases don't pin it but use it
    subprocess.run(
        [pip, "install", "--quiet", "--no-cache-dir", spec, "scipy"],
        check=True,
        cwd=tmpdir,
    )
    return os.path.join(venv_dir, "bin", "python")


def _create_local_venv(tmpdir):
    """Create a venv and install the local porems source tree into it."""
    venv_dir = os.path.join(tmpdir, "local_venv")
    # Find the project root (parent of this tests/ directory)
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print(f"  Creating venv at {venv_dir} ...", flush=True)
    venv.create(venv_dir, with_pip=True, clear=True)

    pip = os.path.join(venv_dir, "bin", "pip")
    print(f"  Installing local build from {project_root} ...", flush=True)
    subprocess.run(
        [pip, "install", "--quiet", project_root],
        check=True,
    )
    return os.path.join(venv_dir, "bin", "python")


# ---------------------------------------------------------------------------
# Run benchmarks
# ---------------------------------------------------------------------------

def run_benchmarks(python_bin, script, cwd=None):
    """Run bench script in the given interpreter; return dict of results.

    cwd should be a neutral directory (not the project root) so that Python's
    ``-c`` mode does not add the project source tree to sys.path[0].
    """
    result = subprocess.run(
        [python_bin, "-c", script],
        capture_output=True,
        text=True,
        cwd=cwd,
    )
    if result.returncode != 0:
        # Print stderr so failures are visible
        print(f"    stderr: {result.stderr.strip()[:400]}", file=sys.stderr)
        raise RuntimeError(f"Benchmark subprocess failed (rc={result.returncode})")
    # Last line of stdout is the JSON payload
    for line in reversed(result.stdout.strip().splitlines()):
        if line.startswith("{"):
            return json.loads(line)
    raise ValueError(f"No JSON found in output:\n{result.stdout}")


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def _fmt(val):
    if isinstance(val, float):
        return f"{val:.3f}s"
    return str(val)


def print_table(cases, local_res, pypi_res, pypi_version):
    labels = [c[0] for c in cases]
    col_w = max(len(lb) for lb in labels) + 2
    ver_local = "local (optimized)"
    ver_pypi  = f"PyPI {pypi_version}"

    header = f"{'Benchmark':<{col_w}}  {ver_local:>18}  {ver_pypi:>14}  {'speedup':>8}"
    sep = "-" * len(header)

    print()
    print(sep)
    print(header)
    print(sep)

    for label in labels:
        t_local = local_res.get(label)
        t_pypi  = pypi_res.get(label)

        if isinstance(t_local, float) and isinstance(t_pypi, float) and t_pypi > 0:
            speedup = f"{t_pypi / t_local:.2f}×"
            faster  = "  ✓" if t_pypi > t_local else "  ="
        else:
            speedup = "n/a"
            faster  = ""

        print(
            f"{label:<{col_w}}  {_fmt(t_local):>18}  {_fmt(t_pypi):>14}  {speedup:>8}{faster}"
        )

    print(sep)

    # Summary: geometric mean speedup over numeric results
    ratios = [
        pypi_res[lb] / local_res[lb]
        for lb in labels
        if isinstance(local_res.get(lb), float)
        and isinstance(pypi_res.get(lb), float)
        and local_res[lb] > 0
    ]
    if ratios:
        import math
        gmean = math.exp(sum(math.log(r) for r in ratios) / len(ratios))
        print(f"  Geometric mean speedup: {gmean:.2f}×")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pypi", metavar="VERSION", default="0.3.0",
                        help="PyPI version to compare against (default: 0.3.0)")
    parser.add_argument("--quick", action="store_true",
                        help="Run a smaller subset of benchmarks")
    args = parser.parse_args()

    cases = _BENCH_CASES_QUICK if args.quick else _BENCH_CASES_FULL
    script = _make_bench_script(cases)

    with tempfile.TemporaryDirectory(prefix="porems_bench_") as tmpdir:
        # ---- PyPI version -------------------------------------------------
        pypi_version = args.pypi
        print(f"\n[1/2] Installing PyPI porems ({pypi_version}) ...", flush=True)
        pypi_python = _create_pypi_venv(args.pypi, tmpdir)

        print(f"\n[1/2] Running benchmarks with PyPI version ...", flush=True)
        pypi_res = run_benchmarks(pypi_python, script, cwd=tmpdir)

        # Detect the actual installed version from the venv via pip show
        pip_pypi = os.path.join(os.path.dirname(pypi_python), "pip")
        ver_out = subprocess.run(
            [pip_pypi, "show", "porems"],
            capture_output=True, text=True, cwd=tmpdir,
        )
        detected = pypi_version
        for line in ver_out.stdout.splitlines():
            if line.startswith("Version:"):
                detected = line.split(":", 1)[1].strip()
                break

        # ---- Local version ------------------------------------------------
        print(f"\n[2/2] Installing local build ...", flush=True)
        local_python = _create_local_venv(tmpdir)
        print(f"\n[2/2] Running benchmarks with local build ...", flush=True)
        local_res = run_benchmarks(local_python, script, cwd=tmpdir)

        # ---- Report -------------------------------------------------------
        print_table(cases, local_res, pypi_res, detected)


if __name__ == "__main__":
    main()
