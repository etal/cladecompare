"""Tests of the input, cleaning, site-selection and output functions in core."""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pytest
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import false_discovery_control

from cladecomparelib import InputError, core

INPUTS = Path(__file__).resolve().parent / "data" / "inputs"


# --- CMA input ---


def write_cma(tmp_path, text):
    path = tmp_path / "test.cma"
    path.write_text(text)
    return str(path)


@pytest.mark.parametrize(
    "seq_line",
    [
        "{()ACDE()}*",
        "MENINNI{()ACDE()}ASATT*",
        "{(HY)ACDE(AA)}*",
        "MK{(HY)ACDE(AA)}LV*",
    ],
)
def test_cma_flanking_residues_are_not_aligned(tmp_path, seq_line):
    fname = write_cma(
        tmp_path, f">rec1 description\n{seq_line}\n>rec2\n{{()FGHI()}}*\n"
    )
    assert core.cma_records(fname) == [("rec1", "ACDE"), ("rec2", "FGHI")]
    aln = core.cma_blocks(fname)
    assert [(rec.id, str(rec.seq)) for rec in aln] == [
        ("rec1", "ACDE"),
        ("rec2", "FGHI"),
    ]


def test_cma_inserts_kept_by_records_and_dropped_by_blocks(tmp_path):
    fname = write_cma(tmp_path, ">rec1\n{()AC.DEfgh-()}*\n")
    assert core.cma_records(fname) == [("rec1", "AC.DEfgh-")]
    assert [str(rec.seq) for rec in core.cma_blocks(fname)] == ["AC.DE-"]


@pytest.mark.parametrize(
    "text, record",
    [
        (">rec1\n{()ACDE\nFGHI()}*\n", "rec1"),  # sequence wrapped over lines
        (">rec1\n{()ACDE\n", "rec1"),  # no closing '()}*'
        (">rec1\n{()ACDE()}*\n>rec2\n", "rec2"),  # header with no sequence line
        (">rec1\n{()ACDE()}* trailing\n", "rec1"),  # text after the block
    ],
)
def test_cma_malformed_sequence_line_names_record(tmp_path, text, record):
    fname = write_cma(tmp_path, text)
    with pytest.raises(InputError, match=rf"record {record}\b"):
        core.cma_records(fname)
    with pytest.raises(InputError, match=rf"record {record}\b"):
        core.cma_blocks(fname)


@pytest.mark.parametrize("name", ["scttlre-domain", "cdc2-domain"])
def test_cma_blocks_match_consensus_column_fasta(name):
    # tests/data/make_inputs.py derives the .aln.fasta files from the CMA
    # files independently of cma_blocks.
    blocks = core.cma_blocks(str(INPUTS / (name + ".fasta_aln.cma")))
    fasta = AlignIO.read(str(INPUTS / (name + ".aln.fasta")), "fasta")
    assert [(rec.id, str(rec.seq)) for rec in blocks] == [
        (rec.id, str(rec.seq)) for rec in fasta
    ]


# --- Unreadable input ---


@pytest.mark.parametrize(
    "name, content, fmt, message",
    [
        ("aln.fasta", None, "fasta", ": No such file or directory"),
        ("aln.cma", None, "fasta", ": No such file or directory"),
        ("aln.fasta", b"\xca\xfe\n", "fasta", ": it is not a text file"),
        ("aln.cma", b">a\n\xca\xfe\n", "fasta", ": it is not a text file"),
        ("aln.fasta", b"", "fasta", " as a fasta alignment: No records found"),
        (
            "aln.fasta",
            b">a\nACDE\n>b\nACD\n",
            "fasta",
            " as a fasta alignment: Sequences must all be the same length",
        ),
        ("aln.fasta", b">a\nACDE\n", "stockholm", " as a stockholm alignment: "),
        ("aln.fasta", b"hello world\n", "fasta", " as a fasta alignment: "),
    ],
)
def test_read_aln_names_a_file_it_cannot_read(tmp_path, name, content, fmt, message):
    """The message names the file, and keeps a parser's own message to one
    line."""
    path = tmp_path / name
    if content is not None:
        path.write_bytes(content)
    with pytest.raises(InputError) as excinfo:
        core.read_aln(str(path), fmt)
    assert str(excinfo.value).startswith(f"Cannot read {path}{message}")
    assert "\n" not in str(excinfo.value)


# Biopython warns about each record of these files that it does not recognize.
@pytest.mark.filterwarnings("ignore::Bio.BiopythonWarning")
@pytest.mark.parametrize("text", ["", "HEADER    KINASE\n", "not a structure\n"])
def test_read_pdb_seq_rejects_a_file_without_atoms(tmp_path, text):
    path = tmp_path / "1abc.pdb"
    path.write_text(text)
    with pytest.raises(InputError, match="as a PDB file"):
        with core.read_pdb_seq(str(path)):
            pass


def test_hmmalign_names_a_missing_profile_before_running(tmp_path, monkeypatch):
    """hmmalign is not run on a file that cannot be read."""

    def run(*args):
        raise AssertionError(f"ran {args}")

    monkeypatch.setattr(core, "call_quiet", run)
    seqs = tmp_path / "seqs.fasta"
    seqs.write_text(">a\nACDE\n")
    with pytest.raises(InputError, match="Cannot read missing.hmm: No such file"):
        core.hmmalign_a2m("missing.hmm", str(seqs))


# --- clean_alignments ---


def aln(*rows):
    return MultipleSeqAlignment(
        [SeqRecord(Seq(seq), id=id_, description="") for id_, seq in rows]
    )


def rows(alignment):
    return [(rec.id, str(rec.seq)) for rec in alignment]


def test_clean_removes_foreground_ids_from_background():
    fg, bg = core.clean_alignments(
        aln(("a", "ACDE"), ("b", "ACDF")),
        aln(("a", "GCDE"), ("c", "ACDH"), ("d", "KCDH")),
    )
    assert rows(fg) == [("a", "ACDE"), ("b", "ACDF")]
    assert rows(bg) == [("c", "ACDH"), ("d", "KCDH")]


def test_clean_purges_identical_sequences_across_both_sets():
    # A background sequence identical to a foreground one goes even under a
    # different ID, as does a repeat within either set.
    fg, bg = core.clean_alignments(
        aln(("a", "ACDE"), ("b", "ACDE"), ("c", "ACDF")),
        aln(("x", "ACDF"), ("y", "ACDH"), ("z", "ACDH"), ("w", "KCDH")),
    )
    assert rows(fg) == [("a", "ACDE"), ("c", "ACDF")]
    assert rows(bg) == [("y", "ACDH"), ("w", "KCDH")]


def test_clean_forgets_sequences_between_calls():
    # Each call purges duplicates only among its own inputs; the multi-
    # alignment mode calls it once per foreground with the same sequences.
    first = core.clean_alignments(aln(("a", "ACDE")), aln(("b", "ACDF")))
    second = core.clean_alignments(aln(("a", "ACDE")), aln(("b", "ACDF")))
    assert [rows(part) for part in second] == [rows(part) for part in first]
    assert rows(second[0]) == [("a", "ACDE")]


@pytest.mark.parametrize("gap", ["-", "."])
def test_clean_drops_columns_that_are_gaps_in_both_sets(gap):
    # Column 2 is a gap everywhere; column 4 only in the foreground. '.', the
    # gap of an insert column, counts as a gap as '-' does.
    def with_gap(*rows_):
        return [(id_, seq.replace("-", gap)) for id_, seq in rows_]

    fg, bg = core.clean_alignments(
        aln(*with_gap(("a", "A-C-E"), ("b", "A-D-E"))),
        aln(*with_gap(("c", "A-CGE"), ("d", "K-D-E"))),
    )
    assert rows(fg) == with_gap(("a", "AC-E"), ("b", "AD-E"))
    assert rows(bg) == with_gap(("c", "ACGE"), ("d", "KD-E"))


def test_clean_keeps_profile_columns_that_are_gaps_in_both_sets():
    # Under a profile each column is a match state, occupied or not.
    fg, bg = core.clean_alignments(
        aln(("a", "A-C-E"), ("b", "A-D-E")),
        aln(("c", "A-CGE"), ("d", "K-D-E")),
        profile_aligned=True,
    )
    assert rows(fg) == [("a", "A-C-E"), ("b", "A-D-E")]
    assert rows(bg) == [("c", "A-CGE"), ("d", "K-D-E")]


@pytest.mark.parametrize(
    "fg, bg, message",
    [
        ([], [("b", "ACDE")], "Foreground set is empty"),
        ([("a", "ACDE")], [], "Background set is empty"),
        ([("a", "ACDE"), ("b", "ACDF")], [("a", "ACDE")], "subset"),
        # Identical under other IDs: the background is empty only after the
        # duplicate sequences are purged.
        ([("a", "ACDE"), ("b", "ACDF")], [("c", "ACDF"), ("d", "ACDE")], "identical"),
    ],
)
def test_clean_rejects_empty_sets(fg, bg, message):
    with pytest.raises(InputError, match=message):
        core.clean_alignments(aln(*fg), aln(*bg))


# --- Letter case ---


def write_fasta(path, rows_):
    path.write_text("".join(f">{id_}\n{seq}\n" for id_, seq in rows_))
    return str(path)


def lowercased(src, dest):
    """Copy a FASTA alignment with its sequence lines in lowercase."""
    lines = Path(src).read_text().splitlines()
    dest.write_text(
        "".join((ln if ln.startswith(">") else ln.lower()) + "\n" for ln in lines)
    )
    return str(dest)


@pytest.mark.parametrize("strategy", sorted(core.STRATEGIES))
def test_alignment_without_uppercase_gives_the_results_of_uppercase(tmp_path, strategy):
    """A wholly lowercase file marks no inserts, so it is compared as aligned,
    alone and against an uppercase background."""
    module = core.STRATEGIES[strategy][1]
    fg = str(INPUTS / "scttlre-domain.aln.fasta")
    bg = str(INPUTS / "cdc2-domain.aln.fasta")
    fg_lower = lowercased(fg, tmp_path / "fg.fasta")

    if hasattr(module, "compare_one"):  # hypg compares two alignments only
        _, expected = core.process_one(core.read_aln(fg, "fasta"), module, False)
        _, hits = core.process_one(core.read_aln(fg_lower, "fasta"), module, False)
        assert hits == expected

    *_, expected = core.process_pair(
        core.read_aln(fg, "fasta"), core.read_aln(bg, "fasta"), module, False
    )
    *_, hits = core.process_pair(
        core.read_aln(fg_lower, "fasta"), core.read_aln(bg, "fasta"), module, False
    )
    assert hits == expected


def test_lowercase_copy_of_a_foreground_sequence_is_purged(tmp_path):
    fg = write_fasta(tmp_path / "fg.fasta", [("a", "ACDE"), ("b", "ACDF")])
    bg = write_fasta(tmp_path / "bg.fasta", [("x", "acdf"), ("y", "kcdh")])
    _, bg_clean = core.clean_alignments(
        core.read_aln(fg, "fasta"), core.read_aln(bg, "fasta")
    )
    assert rows(bg_clean) == [("y", "KCDH")]


@pytest.mark.parametrize("strategy", sorted(core.STRATEGIES))
def test_set_left_without_uppercase_is_compared_as_aligned(strategy):
    """The background's only uppercase row duplicates a foreground ID and is
    removed, so its other rows no longer mark anything by case."""
    module = core.STRATEGIES[strategy][1]
    fg = [("a", "SCDEF"), ("b", "SCDEW"), ("c", "TCDEF")]
    bg = [("x", "scdkf"), ("y", "gcdew"), ("z", "acdkf")]
    *_, expected = core.process_pair(
        aln(*fg), aln(*[(id_, seq.upper()) for id_, seq in bg]), module, False
    )
    *_, hits = core.process_pair(aln(*fg), aln(("a", "SCDEF"), *bg), module, False)
    assert hits == expected


def test_lowercase_beside_uppercase_marks_an_insert_column(tmp_path):
    """A file with uppercase residues keeps its case, and a column of
    lowercase residues alone is an insert column, reported with p = 1."""
    fname = write_fasta(
        tmp_path / "aln.fasta", [("a", "SCdEF"), ("b", "SCeEW"), ("c", "TC-EF")]
    )
    alignment = core.read_aln(fname, "fasta")
    assert rows(alignment)[0] == ("a", "SCdEF")
    _, hits = core.process_one(alignment, core.STRATEGIES["urn"][1], False)
    assert hits[2] == ("-", "_", 1.0, "-")


@pytest.mark.parametrize("weight", [False, True])
@pytest.mark.parametrize("strategy", sorted(core.STRATEGIES))
def test_mixed_case_column_gives_the_results_of_uppercase(strategy, weight):
    """A column that mixes cases is counted as uppercase, as its consensus is
    read. K occurs nowhere else in the foreground in uppercase; the c in
    column 2 moves the sequence weights unless it is uppercased first; and in
    column 6 the stop, *, is the aligned character beside lowercase s."""
    module = core.STRATEGIES[strategy][1]
    fg = [("a", "kCDEFs"), ("b", "kcDEWs"), ("c", "ACsEF*"), ("d", "SCDEFs")]
    bg = [("x", "GCDKFS"), ("y", "ACDEWA"), ("z", "TCDKFS")]
    fg_upper = [(id_, seq.upper()) for id_, seq in fg]

    if hasattr(module, "compare_one"):  # hypg compares two alignments only
        _, expected = core.process_one(aln(*fg_upper), module, weight)
        _, hits = core.process_one(aln(*fg), module, weight)
        assert hits == expected

    *_, expected = core.process_pair(aln(*fg_upper), aln(*bg), module, weight)
    *_, hits = core.process_pair(aln(*fg), aln(*bg), module, weight)
    assert hits == expected


# --- Untested columns ---


@pytest.mark.parametrize("strategy", sorted(core.STRATEGIES))
def test_consensus_of_unknown_residues_is_not_tested(strategy):
    """The counts skip X, so a column whose consensus is X, in either set, has
    no consensus residue to count; it is reported with p = 1, not with a value
    from the residues beside the X."""
    module = core.STRATEGIES[strategy][1]
    fg = [("a", "XDAKS"), ("b", "XDARS"), ("c", "SEAKT")]
    bg = [
        ("p", "GDAKX"),
        ("q", "ADARX"),
        ("r", "SEAKX"),
        ("s", "TDCKA"),
        ("u", "ADAWA"),
    ]
    if hasattr(module, "compare_one"):  # hypg compares two alignments only
        _, hits = core.process_one(aln(*fg), module, False)
        assert hits[0] == ("X", "_", 1.0, "X")
    *_, hits = core.process_pair(aln(*fg), aln(*bg), module, False)
    assert (hits[0][0], hits[0][2]) == ("X", 1.0)
    assert (hits[4][1], hits[4][2]) == ("X", 1.0)


def test_process_one_refuses_a_strategy_for_two_alignments():
    alignment = aln(("a", "ACDE"), ("b", "ACDF"))
    with pytest.raises(InputError, match="hypg compares two alignments"):
        core.process_one(alignment, core.STRATEGIES["hypg"][1], False)


# --- top_hits ---


def hits_from(pvalues):
    return [("A", "C", p, "A") for p in pvalues]


# The most sites that top_hits returns.
MAX_SITES = 50


def bh_expected(pvalues, alpha):
    """Positions of the first min(k*, MAX_SITES) sites by the BH procedure."""
    pvalues = np.asarray(pvalues)
    k_star = int((false_discovery_control(pvalues, method="bh") <= alpha).sum())
    order = np.argsort(pvalues, kind="stable")
    return [int(idx) + 1 for idx in order[: min(k_star, MAX_SITES)]]


def positions(tophits):
    return [posn for posn, _residue in tophits]


@pytest.mark.parametrize("seed", range(20))
@pytest.mark.parametrize("alpha", [0.005, 0.05, 0.2])
@pytest.mark.parametrize("n_signals", [35, 150])
def test_top_hits_agrees_with_bh(seed, alpha, n_signals):
    # Up to about n_signals strong signals among uniform nulls and gap columns
    # (p = 1). k* exceeds the cap of 50 in every draw with n_signals = 150 and
    # m = 250, and in no other.
    rng = np.random.default_rng(seed)
    m = int(rng.choice([10, 60, 250]))
    pvalues = rng.random(m)
    signals = rng.random(m) < n_signals / m
    pvalues[signals] *= 10.0 ** -rng.uniform(2, 8, signals.sum())
    pvalues[rng.random(m) < 0.3] = 1.0
    expected = bh_expected(pvalues, alpha)
    assert positions(core.top_hits(hits_from(pvalues), alpha)) == expected


def test_top_hits_breaks_ties_by_position_and_returns_tested_residues():
    """Each site is named by the residue that names it in hits, which in column
    2 is not the foreground consensus."""
    hits = [
        ("A", "C", 0.5, "A"),
        ("D", "E", 0.001, "S"),
        ("F", "G", 0.001, "F"),
        ("H", "I", 0.0, "H"),
    ]
    assert core.top_hits(hits, 0.05) == [(4, "H"), (2, "S"), (3, "F")]


def test_top_hits_counts_every_column_as_a_hypothesis():
    # Columns skipped as gaps enter with p = 1 and raise m, so that each
    # rank's threshold alpha * k / m falls.
    assert positions(core.top_hits(hits_from([0.001, 0.02]), 0.05)) == [1, 2]
    assert positions(core.top_hits(hits_from([0.001, 0.02] + [1.0] * 8), 0.05)) == [1]


def test_top_hits_accepts_p_equal_to_the_threshold():
    # p_(2) = alpha * 2 / m and p_(1) = alpha * 1 / m exactly.
    assert positions(core.top_hits(hits_from([0.025, 0.05]), 0.05)) == [1, 2]


@pytest.mark.parametrize(
    "pvalues",
    [
        # k* = 100, and p > alpha * k / m at every rank k <= 50, so a search
        # confined to those ranks accepts nothing.
        pytest.param([0.004] * 100 + [0.5] * 900, id="none-in-the-top-50"),
        # k* = 100 again, and a search confined to ranks k <= 50 stops at
        # rank 10, accepting 10 of the 50 sites that pass.
        pytest.param([1e-5] * 10 + [0.004] * 90 + [0.5] * 900, id="ten-in-the-top-50"),
    ],
)
def test_top_hits_returns_50_sites_when_more_pass(pvalues):
    expected = bh_expected(pvalues, 0.05)
    assert len(expected) == MAX_SITES
    assert positions(core.top_hits(hits_from(pvalues), 0.05)) == expected


# --- Structure residue numbering ---


def test_aln_offset_matches_the_whole_segment_including_x():
    # The chain also begins with AGS, so the residues before the X alone
    # have no unique place.
    assert core.get_aln_offset("AGSAGSXKLMN", "AGSXKL", "chain A of 1xyz.pdb") == 3


def test_aln_offset_ignores_case_deletions_and_insert_padding():
    assert core.get_aln_offset("QQACDEFGHQQ", "ac-D.EfGH", "chain A of 1xyz.pdb") == 2


@pytest.mark.parametrize(
    "full, aln, message",
    [
        ("ACDACD", "ACD", "more than one place"),
        # Overlapping repeats, which str.count does not count twice.
        ("KAKAK", "KAK", "more than one place"),
        ("MXXXK", "XX", "more than one place"),
        ("ACDEFG", "WWWW", "do not occur"),
        ("ACDEFG", "--", "No residues"),
    ],
)
def test_aln_offset_rejects_a_segment_it_cannot_place(full, aln, message):
    with pytest.raises(InputError, match=message) as excinfo:
        core.get_aln_offset(full, aln, "chain A of 1xyz.pdb")
    assert "chain A of 1xyz.pdb " in str(excinfo.value)


def pdb_chain(chain_id, seq, start=1):
    """A chain as read_pdb_seq yields it: the FASTA name and the record."""
    record = SeqRecord(
        Seq(seq),
        id="1XYZ:%s" % chain_id,
        annotations={"chain": chain_id, "start": start},
    )
    return record.id, record


def test_pdb_ref_resnums_numbers_residues_from_the_chain_start():
    # Q101 Q102 A103 C104 D105 E106 F107 G108 H109 Q110 Q111
    chains = dict([pdb_chain("A", "QQACDEFGHQQ", start=101)])
    ref, resnums, inserts = core.pdb_ref_resnums(
        {"1XYZ:A": "ac-DEfGH"}, chains, "1xyz.pdb"
    )
    assert ref.id == "1XYZ:A"
    assert resnums == [("-", None), ("D", 105), ("E", 106), ("G", 108), ("H", 109)]
    assert inserts == [(103, 104), (107, 107)]


def test_pdb_ref_resnums_names_the_chain_it_cannot_place():
    chains = dict([pdb_chain("A", "ACDACD")])
    with pytest.raises(InputError, match=r"chain A of 1xyz.pdb .*more than one place"):
        core.pdb_ref_resnums({"1XYZ:A": "ACD"}, chains, "1xyz.pdb")


# A kinase, A, aligned at three match states, beside a partner, B, aligned
# at one, and an identical copy of the kinase, C; X is not counted.
CHAINS = dict([pdb_chain("A", "KLM"), pdb_chain("B", "E"), pdb_chain("C", "KLMX")])
ROWS = {"1XYZ:A": "KLM-", "1XYZ:B": "-E--", "1XYZ:C": "KLMX"}


@pytest.mark.parametrize("order", [["A", "B", "C"], ["C", "B", "A"], ["B", "C", "A"]])
def test_choose_best_aligned_takes_the_most_matches_then_the_first_id(order):
    """Ties go to the chain ID that sorts first, in any order of the rows."""
    rows = {"1XYZ:" + c: ROWS["1XYZ:" + c] for c in order}
    record, row = core.choose_best_aligned(rows, CHAINS, "1xyz.pdb")
    assert (record.annotations["chain"], row) == ("A", "KLM-")
    del rows["1XYZ:A"]
    record, row = core.choose_best_aligned(rows, CHAINS, "1xyz.pdb")
    assert (record.annotations["chain"], row) == ("C", "KLMX")


@pytest.mark.parametrize("rows", [{}, {"1XYZ:A": "---", "1XYZ:B": "xX-"}])
def test_choose_best_aligned_rejects_a_structure_without_matches(rows):
    """--mapgaps can return no chain; hmmalign can return only gaps and X."""
    with pytest.raises(InputError, match="No chain of 1xyz.pdb has a residue"):
        core.choose_best_aligned(rows, CHAINS, "1xyz.pdb")


# --- Profile-aligned runs ---

# Five match states of a profile, aligned as by --hmm or --mapgaps. State 2
# is a gap in every sequence of every set, and states 3 and 4 hold W in the
# first set and G in the others.
PROFILE_ALIGNED = {
    "fg.fasta": [("a", "A-WWK"), ("b", "C-WWK"), ("c", "D-WWK")],
    "bg.fasta": [("x", "E-GGK"), ("y", "F-GGK"), ("z", "H-GGK")],
    "bg2.fasta": [("p", "I-GGK"), ("q", "L-GGK"), ("r", "M-GGK")],
}
# The structure's residue at each of the five match states.
STRUCTURE = [("M", 10), ("N", 11), ("P", 12), ("Q", 13), ("R", 14)]


def profile_run(tmp_path, monkeypatch, option, background):
    """Run cladecompare with a profile and a structure on the sets above."""
    chain = SeqRecord(
        Seq("MNPQR"), id="1abc:A", annotations={"chain": "A", "idcode": "1abc"}
    )
    monkeypatch.setattr(
        core,
        "%s_align_and_read" % option,
        lambda _profile, name: aln(*PROFILE_ALIGNED[name]),
    )
    monkeypatch.setattr(
        core, "pdb_%s" % option, lambda _profile, _pdb: (chain, STRUCTURE, [])
    )
    monkeypatch.chdir(tmp_path)
    core.process_args(
        argparse.Namespace(
            foreground="fg.fasta",
            background=background,
            format="fasta",
            hmm="x.hmm" if option == "hmm" else None,
            mapgaps="x" if option == "mapgaps" else None,
            pdb=["1abc.pdb"],
            strategy="gtest",
            weight=False,
            alpha=0.05,
            output="x.out",
            pattern="x.pttrn",
            pmlout=None,
        )
    )


def column_labels(out_file):
    # Column 1 holds three residue types in each set, so its consensus
    # depends on how biofrills breaks the tie; only its number is checked.
    columns = out_file.read_text().splitlines()
    assert columns[0].split(" : ")[0].endswith(" 1")
    return [line.split(" : ")[0] for line in columns[1:]]


@pytest.mark.parametrize("option", ["hmm", "mapgaps"])
def test_profile_run_numbers_columns_by_match_state(tmp_path, monkeypatch, option):
    """A match state that no sequence occupies keeps its column, so each site
    draws the structure's residue at its own match state."""
    profile_run(tmp_path, monkeypatch, option, ["bg.fasta"])
    assert column_labels(tmp_path / "x.out") == [
        "- (-) 2",
        "W (G) 3",
        "W (G) 4",
        "K (K) 5",
    ]
    sites = (tmp_path / "x.pttrn").read_text().strip().split(":")[1].split(",")
    assert sorted(sites) == ["W3", "W4"]
    script = (tmp_path / "1abc.pdb.pml").read_text()
    drawn = re.findall(r"^create (\w+)-1abcA", script, re.MULTILINE)
    assert sorted(drawn) == ["P12", "Q13"]


def test_profile_run_of_three_sets_numbers_columns_by_match_state(
    tmp_path, monkeypatch
):
    """Each set compared with the other two keeps every match state."""
    background = ["bg.fasta", "bg2.fasta"]
    profile_run(tmp_path, monkeypatch, "hmm", background)
    for name in ["fg.fasta"] + background:
        labels = column_labels(tmp_path / (name + ".out"))
        assert len(labels) == 4
        assert labels[0] == "- (-) 2"


def test_process_output_rejects_a_structure_of_another_width(tmp_path):
    """Sites that cannot name match states stop the run before any output."""
    chain = SeqRecord(Seq("MNPQ"), id="1abc:A", annotations={"chain": "A"})
    hits = [("W", "G", 1e-6)] * 5
    pdb_data = [("1abc.pdb", chain, STRUCTURE[:4], [])]
    output = tmp_path / "x.out"
    with pytest.raises(InputError, match="1abc.pdb is aligned to 4 match states"):
        core.process_output(None, None, hits, 0.05, str(output), None, pdb_data)
    assert not output.exists()


# --- External programs ---


def test_call_quiet_names_a_missing_program():
    with pytest.raises(RuntimeError, match="Is no-such-program-x9 installed"):
        core.call_quiet("no-such-program-x9", "--help")


def test_call_quiet_reports_the_output_of_a_failed_command():
    with pytest.raises(RuntimeError, match=r"(?s)partial output.*went wrong"):
        core.call_quiet(
            sys.executable,
            "-c",
            "import sys; print('partial output'); sys.exit('went wrong')",
        )
