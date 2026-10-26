"""Tests of the .out and .pttrn readers and the heat-map classes."""

import io
import math

import pytest

from cladecomparelib import core, report


def test_parse_out_reads_what_write_pvalues_writes(tmp_path):
    hits = [
        ("A", "C", 0.5),
        ("D", "-", 0.003),
        ("-", "E", 1.0),
        ("F", "G", 2.6e-74),
        ("H", "I", 0.0),
        ("K", "L", 0.001),
        ("M", "N", 0.005),
    ]
    buf = io.StringIO()
    core.write_pvalues(hits, buf, 0.005)
    path = tmp_path / "test.out"
    path.write_text(buf.getvalue())
    # One star per order of magnitude below 1, only for 0 < p < alpha.
    assert list(report.parse_out(str(path))) == [
        ("A", "C", 1, 0.5, 0),
        ("D", "-", 2, 0.003, 2),
        ("-", "E", 3, 1.0, 0),
        ("F", "G", 4, 2.6e-74, 73),
        ("H", "I", 5, 0.0, 0),
        ("K", "L", 6, 0.001, 3),
        ("M", "N", 7, 0.005, 0),
    ]


def test_parse_out_reads_nan(tmp_path):
    path = tmp_path / "test.out"
    path.write_text("A (C) 1 : prob=nan\t\n")
    [(_fg, _bg, _posn, prob, nstars)] = report.parse_out(str(path))
    assert math.isnan(prob) and nstars == 0


@pytest.mark.parametrize(
    "pval, css",
    [
        (1.0, "b5"),
        (0.9999, "b4"),
        (0.99, "b3"),
        (0.9, "b2"),
        (0.6, "b1"),
        (0.1, "wt"),
        (0.05, "r1"),
        (1e-3, "r2"),
        (1e-6, "r3"),
        (1e-8, "r4"),
        (1e-9, "r9"),
        (2.6e-74, "r9"),
        (0.0, "r9"),
    ],
)
def test_p2class_scale(pval, css):
    assert report.p2class("A", "C", pval) == css


def test_p2class_marks_gaps_before_p_values():
    assert report.p2class("-", "C", 0.0) == "del"
    assert report.p2class("A", "-", 0.0) == "ins"


def test_parse_pattern_returns_sorted_positions(tmp_path):
    path = tmp_path / "test.pttrn"
    path.write_text("1:D34,G144,T149,V1")
    assert report.parse_pattern(str(path)) == [1, 34, 144, 149]


def test_parse_pattern_accepts_an_empty_pattern(tmp_path):
    path = tmp_path / "test.pttrn"
    path.write_text("1:")
    assert report.parse_pattern(str(path)) == []
