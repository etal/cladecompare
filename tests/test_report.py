"""Tests of the .out and .pttrn readers and the heat-map classes."""

import io

import pytest

from cladecomparelib import InputError, core, report


def test_parse_out_reads_what_write_pvalues_writes(tmp_path):
    # The .out file names the consensus residues, whichever residue names the
    # site in the pattern file.
    hits = [
        ("A", "C", 0.5, "A"),
        ("D", "-", 0.003, "S"),
        ("-", "E", 1.0, "-"),
        ("F", "G", 2.6e-74, "F"),
        ("H", "I", 0.0, "H"),
        ("K", "L", 0.001, "K"),
        ("M", "N", 0.005, "M"),
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


@pytest.mark.parametrize(
    "text, message",
    [
        ("", "holds no lines of cladecompare output"),
        ("\n\n", "holds no lines of cladecompare output"),
        ("A (C) 1 : prob=0.5\t\nA (C) 2 : prob=nan\t\n", "line 2: the value 'nan'"),
        ("A (C) 1 : prob=inf\t\n", "line 1: the value 'inf'"),
        ("A (C) 1 : prob=0.5\t\ngarbage line\n", "line 2: not a line of"),
        ("A (C) x : prob=0.5\t\n", "line 1: not a line of"),
        ("A (C) 1 : prob=0.003**\n", "line 1: not a line of"),
    ],
)
def test_parse_out_rejects_a_file_it_cannot_draw(tmp_path, text, message):
    path = tmp_path / "test.out"
    path.write_text(text)
    with pytest.raises(InputError, match=message):
        report.parse_out(str(path))


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


@pytest.mark.parametrize(
    "text, message",
    [
        ("", "is not a pattern file"),
        ("D34,G144", "is not a pattern file"),
        ("1:D34,G", "the site 'G' has no column number"),
    ],
)
def test_parse_pattern_rejects_what_write_mcbpps_does_not_write(
    tmp_path, text, message
):
    path = tmp_path / "test.pttrn"
    path.write_text(text)
    with pytest.raises(InputError, match=message):
        report.parse_pattern(str(path))


def test_pattern_marks_columns_by_their_numbers(tmp_path):
    path = tmp_path / "test.pttrn"
    path.write_text("1:D12")
    assert report.pattern_marks(str(path), (11, 12, 13), "test.out") == ["", "*", ""]


@pytest.mark.parametrize("site", ["D0", "D4"])
def test_pattern_marks_rejects_a_site_outside_the_columns(tmp_path, site):
    path = tmp_path / "test.pttrn"
    path.write_text("1:" + site)
    with pytest.raises(
        InputError,
        match="which test.out does not have; its columns are numbered 1 to 3",
    ):
        report.pattern_marks(str(path), (1, 2, 3), "test.out")


def test_tr_class_escapes_the_label():
    """tr_class escapes the label, which do_multi takes from a file name."""
    row = report.tr_class("R&D <1>", ["A"], [0.5], ["wt"])
    assert '<td class="label">R&amp;D &lt;1&gt;</td>' in row


def test_do_multi_draws_each_file_as_its_own_row(tmp_path, monkeypatch):
    """Row k of the table holds the residues, classes and p-values of file k.

    The one pattern row marks its columns above all the rows.
    """
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.out").write_text("A (C) 1 : prob=0.5\t\nD (-) 2 : prob=1e-12\t\n")
    (tmp_path / "b.out").write_text("- (E) 1 : prob=1\t\nF (G) 2 : prob=0.001\t\n")
    (tmp_path / "a.pttrn").write_text("1:D2")
    title, contents = report.do_multi(["a.out", "b.out"], "a.pttrn")
    assert title == "a.out, b.out"
    pattern = '<tr>\n<td class="label"></td><td></td><td>*</td>\n</tr>'
    row_a = (
        '<tr>\n<td class="label">a</td>'
        '<td class="wt" title="p=0.5">A</td>'
        '<td class="ins" title="p=1e-12">D</td>\n</tr>'
    )
    row_b = (
        '<tr>\n<td class="label">b</td>'
        '<td class="del" title="p=1.0">-</td>'
        '<td class="r2" title="p=0.001">F</td>\n</tr>'
    )
    assert contents.count("<table>") == 1
    assert contents.index(pattern) < contents.index(row_a) < contents.index(row_b)


@pytest.mark.parametrize(
    "other",
    [
        "A (C) 1 : prob=0.5\t\n",
        "A (C) 2 : prob=0.5\t\nA (C) 3 : prob=0.5\t\n",
    ],
)
def test_do_multi_rejects_files_numbering_other_columns(tmp_path, other):
    first = tmp_path / "a.out"
    first.write_text("A (C) 1 : prob=0.5\t\nA (C) 2 : prob=0.5\t\n")
    second = tmp_path / "b.out"
    second.write_text(other)
    with pytest.raises(InputError, match="do not number the same columns"):
        report.do_multi([str(first), str(second)], None)
