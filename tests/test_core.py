"""Tests of the input, cleaning and site-selection functions in core."""

import sys
from pathlib import Path

import numpy as np
import pytest
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import false_discovery_control

from cladecomparelib import core

INPUTS = Path(__file__).resolve().parent / "data" / "inputs"


# --- CMA input ---

def write_cma(tmp_path, text):
    path = tmp_path / "test.cma"
    path.write_text(text)
    return str(path)


@pytest.mark.parametrize("seq_line", [
    "{()ACDE()}*",
    "MENINNI{()ACDE()}ASATT*",
    "{(HY)ACDE(AA)}*",
    "MK{(HY)ACDE(AA)}LV*",
])
def test_cma_flanking_residues_are_not_aligned(tmp_path, seq_line):
    fname = write_cma(tmp_path,
                      f">rec1 description\n{seq_line}\n>rec2\n{{()FGHI()}}*\n")
    assert core.cma_records(fname) == [("rec1", "ACDE"), ("rec2", "FGHI")]
    aln = core.cma_blocks(fname)
    assert [(rec.id, str(rec.seq)) for rec in aln] == [("rec1", "ACDE"),
                                                       ("rec2", "FGHI")]


def test_cma_inserts_kept_by_records_and_dropped_by_blocks(tmp_path):
    fname = write_cma(tmp_path, ">rec1\n{()AC.DEfgh-()}*\n")
    assert core.cma_records(fname) == [("rec1", "AC.DEfgh-")]
    assert [str(rec.seq) for rec in core.cma_blocks(fname)] == ["AC.DE-"]


@pytest.mark.parametrize("text, record", [
    (">rec1\n{()ACDE\nFGHI()}*\n", "rec1"),     # sequence wrapped over lines
    (">rec1\n{()ACDE\n", "rec1"),               # no closing '()}*'
    (">rec1\n{()ACDE()}*\n>rec2\n", "rec2"),    # header with no sequence line
    (">rec1\n{()ACDE()}* trailing\n", "rec1"),  # text after the block
])
def test_cma_malformed_sequence_line_names_record(tmp_path, text, record):
    fname = write_cma(tmp_path, text)
    with pytest.raises(ValueError, match=rf"record {record}\b"):
        core.cma_records(fname)
    with pytest.raises(ValueError, match=rf"record {record}\b"):
        core.cma_blocks(fname)


@pytest.mark.parametrize("name", ["scttlre-domain", "cdc2-domain"])
def test_cma_blocks_match_consensus_column_fasta(name):
    # tests/data/make_inputs.py derives the .aln.fasta files from the CMA
    # files independently of cma_blocks.
    blocks = core.cma_blocks(str(INPUTS / (name + ".fasta_aln.cma")))
    fasta = AlignIO.read(str(INPUTS / (name + ".aln.fasta")), "fasta")
    assert ([(rec.id, str(rec.seq)) for rec in blocks]
            == [(rec.id, str(rec.seq)) for rec in fasta])


# --- clean_alignments ---

def aln(*rows):
    return MultipleSeqAlignment([SeqRecord(Seq(seq), id=id_, description="")
                                 for id_, seq in rows])


def rows(alignment):
    return [(rec.id, str(rec.seq)) for rec in alignment]


def test_clean_removes_foreground_ids_from_background():
    fg, bg = core.clean_alignments(
        aln(("a", "ACDE"), ("b", "ACDF")),
        aln(("a", "GCDE"), ("c", "ACDH"), ("d", "KCDH")))
    assert rows(fg) == [("a", "ACDE"), ("b", "ACDF")]
    assert rows(bg) == [("c", "ACDH"), ("d", "KCDH")]


def test_clean_purges_identical_sequences_across_both_sets():
    # A background sequence identical to a foreground one goes even under a
    # different ID, as does a repeat within either set.
    fg, bg = core.clean_alignments(
        aln(("a", "ACDE"), ("b", "ACDE"), ("c", "ACDF")),
        aln(("x", "ACDF"), ("y", "ACDH"), ("z", "ACDH"), ("w", "KCDH")))
    assert rows(fg) == [("a", "ACDE"), ("c", "ACDF")]
    assert rows(bg) == [("y", "ACDH"), ("w", "KCDH")]


def test_clean_forgets_sequences_between_calls():
    # Each call purges duplicates only among its own inputs; the multi-
    # alignment mode calls it once per foreground with the same sequences.
    first = core.clean_alignments(aln(("a", "ACDE")), aln(("b", "ACDF")))
    second = core.clean_alignments(aln(("a", "ACDE")), aln(("b", "ACDF")))
    assert [rows(part) for part in second] == [rows(part) for part in first]
    assert rows(second[0]) == [("a", "ACDE")]


def test_clean_drops_columns_that_are_gaps_in_both_sets():
    # Column 2 is a gap everywhere; column 4 only in the foreground.
    fg, bg = core.clean_alignments(
        aln(("a", "A-C-E"), ("b", "A-D-E")),
        aln(("c", "A-CGE"), ("d", "K-D-E")))
    assert rows(fg) == [("a", "AC-E"), ("b", "AD-E")]
    assert rows(bg) == [("c", "ACGE"), ("d", "KD-E")]


@pytest.mark.parametrize("fg, bg, message", [
    ([], [("b", "ACDE")], "Foreground set is empty"),
    ([("a", "ACDE")], [], "Background set is empty"),
    ([("a", "ACDE"), ("b", "ACDF")], [("a", "ACDE")], "subset"),
])
def test_clean_rejects_empty_sets(fg, bg, message):
    with pytest.raises(ValueError, match=message):
        core.clean_alignments(aln(*fg), aln(*bg))


# --- top_hits ---

def hits_from(pvalues):
    return [("A", "C", p) for p in pvalues]


# The most sites that top_hits returns.
MAX_SITES = 50


def bh_expected(pvalues, alpha):
    """Positions of the first min(k*, MAX_SITES) sites by the BH procedure."""
    pvalues = np.asarray(pvalues)
    k_star = int((false_discovery_control(pvalues, method="bh")
                  <= alpha).sum())
    order = np.argsort(pvalues, kind="stable")
    return [int(idx) + 1 for idx in order[:min(k_star, MAX_SITES)]]


def positions(tophits):
    return [posn for posn, _faa, _baa in tophits]


@pytest.mark.parametrize("seed", range(20))
@pytest.mark.parametrize("alpha", [0.005, 0.05, 0.2])
def test_top_hits_agrees_with_bh_when_at_most_50_pass(seed, alpha):
    # Up to about 35 strong signals among uniform nulls and gap columns
    # (p = 1), so that k* stays below 50; the k* > 50 case is tested below.
    rng = np.random.default_rng(seed)
    m = int(rng.choice([10, 60, 250]))
    pvalues = rng.random(m)
    signals = rng.random(m) < 35 / m
    pvalues[signals] *= 10.0 ** -rng.uniform(2, 8, signals.sum())
    pvalues[rng.random(m) < 0.3] = 1.0
    expected = bh_expected(pvalues, alpha)
    assert len(expected) < MAX_SITES
    assert positions(core.top_hits(hits_from(pvalues), alpha)) == expected


def test_top_hits_breaks_ties_by_position_and_returns_residues():
    hits = [("A", "C", 0.5), ("D", "E", 0.001), ("F", "G", 0.001),
            ("H", "I", 0.0)]
    assert core.top_hits(hits, 0.05) == [(4, "H", "I"), (2, "D", "E"),
                                         (3, "F", "G")]


def test_top_hits_counts_every_column_as_a_hypothesis():
    # Columns skipped as gaps enter with p = 1 and raise m, so that each
    # rank's threshold alpha * k / m falls.
    assert positions(core.top_hits(hits_from([0.001, 0.02]), 0.05)) == [1, 2]
    assert positions(core.top_hits(hits_from([0.001, 0.02] + [1.0] * 8),
                                   0.05)) == [1]


def test_top_hits_accepts_p_equal_to_the_threshold():
    # p_(2) = alpha * 2 / m and p_(1) = alpha * 1 / m exactly.
    assert positions(core.top_hits(hits_from([0.025, 0.05]), 0.05)) == [1, 2]


@pytest.mark.xfail(strict=True, reason="the step-up search examines only the "
                   "50 most significant ranks, so it misses k* > 50")
def test_top_hits_returns_50_sites_when_more_pass():
    # k* = 100, but p > alpha * k / m at every rank k <= 50.
    pvalues = [0.004] * 100 + [0.5] * 900
    expected = bh_expected(pvalues, 0.05)
    assert len(expected) == MAX_SITES
    assert positions(core.top_hits(hits_from(pvalues), 0.05)) == expected


# --- External programs ---

def test_call_quiet_names_a_missing_program():
    with pytest.raises(RuntimeError, match="Is no-such-program-x9 installed"):
        core.call_quiet("no-such-program-x9", "--help")


def test_call_quiet_reports_the_output_of_a_failed_command():
    with pytest.raises(RuntimeError, match=r"(?s)partial output.*went wrong"):
        core.call_quiet(sys.executable, "-c",
                        "import sys; print('partial output'); "
                        "sys.exit('went wrong')")
