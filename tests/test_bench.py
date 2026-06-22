"""Benchmark tests for pore generation speed.

Run with:  pytest tests/test_bench.py -v --tb=short
These are marked slow so they are excluded from the fast unit-test run.
"""
import time
import pytest
import porems as pms


pytestmark = pytest.mark.slow


def _time(fn):
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


@pytest.mark.parametrize("size,diam,label", [
    ([4, 4, 4],  2.0, "small"),
    ([6, 6, 6],  3.0, "medium"),
    ([8, 8, 8],  4.0, "large"),
])
def test_bench_cylinder_generation(size, diam, label):
    """PoreCylinder end-to-end generation time must be under 120 s."""
    elapsed = _time(lambda: pms.PoreCylinder(size, diam, res=0).finalize())
    print(f"\nPoreCylinder({size}, d={diam}) [{label}]: {elapsed:.2f}s")
    assert elapsed < 120, f"Generation took {elapsed:.1f}s — expected < 120s"


def test_bench_slit_generation():
    elapsed = _time(lambda: pms.PoreSlit([6, 6, 6], 3.0, res=0).finalize())
    print(f"\nPoreSlit([6,6,6], h=3.0): {elapsed:.2f}s")
    assert elapsed < 120


def test_bench_multi_channel_generation():
    elapsed = _time(lambda: pms.PoreMultiChannel([8, 4, 6], 2, 2.0, res=0).finalize())
    print(f"\nPoreMultiChannel([8,4,6], 2ch, d=2.0): {elapsed:.2f}s")
    assert elapsed < 180


def test_bench_cone_generation():
    elapsed = _time(lambda: pms.PoreCone([6, 6, 8], 2.0, 4.0, res=0).finalize())
    print(f"\nPoreCone([6,6,8], d=2→4): {elapsed:.2f}s")
    assert elapsed < 120


def test_bench_pattern_generation():
    """BetaCristobalit pattern generation for various block sizes."""
    for size, label in [([4,4,4], "4³"), ([8,8,8], "8³"), ([11,11,10], "11×11×10")]:
        elapsed = _time(lambda s=size: pms.BetaCristobalit().generate(s, "z"))
        print(f"\nBetaCristobalit.generate({label}): {elapsed:.2f}s")
        assert elapsed < 60, f"Pattern generation {label} took {elapsed:.1f}s"
