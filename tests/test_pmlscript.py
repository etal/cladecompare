"""Tests of pmlscript: PDB file name parsing and the residues drawn."""

import logging

import pytest

from cladecomparelib import InputError, pmlscript


@pytest.mark.parametrize(
    "fname, expected",
    [
        ("1hck.pdb", ("1hck", None)),
        ("dir/00_1hck.pdb", ("1hck", None)),
        ("1hckA.pdb", ("1hck", "A")),
        ("1hck_a.pdb", ("1hck", "A")),
        ("1hck_chain_A.pdb", ("1hck", None)),
        ("pdb1hck.ent", ("1hck", None)),
    ],
)
def test_pdbid_from_fname(fname, expected):
    assert pmlscript.pdbid_from_fname(fname) == expected


@pytest.mark.parametrize(
    "fname",
    ["tmpq355bkts", "x.pdb", "1HCKA_.pdb", "1hc_A.pdb", "1FIN-AB-CA.pdb", "1hck.ent"],
)
def test_pdbid_from_fname_rejects_other_names(fname):
    """Names that cannot be read raise InputError, even under python -O.

    An assertion would not run under python -O, and each name would then
    give a wrong PDB ID.
    """
    with pytest.raises(InputError, match="Cannot read a PDB ID"):
        pmlscript.pdbid_from_fname(fname)


# The five columns of a chain that begins at residue 5 and has a deletion at
# column 3, as core.aln_resnums_inserts returns them for the row "QK-GP".
MAPPING = [("Q", 5), ("K", 6), ("-", None), ("G", 7), ("P", 8)]


def residues(sites):
    return "".join(pmlscript.make_residues(MAPPING, sites, "1hckA", "A"))


def test_make_residues_draws_the_first_column():
    """Site 1 is the residue at column 1, not at column 2."""
    script = residues([1])
    assert "create Q5-1hckA, 1hckA and resi 5 and chain A" in script
    assert "K6" not in script


def test_make_residues_draws_the_last_column():
    """The last column is a valid site."""
    assert "create P8-1hckA, 1hckA and resi 8 " in residues([5])


def test_make_residues_skips_a_deletion_column(caplog):
    """A column the chain deletes is logged and the next site still drawn."""
    caplog.set_level(logging.WARNING, logger="cladecomparelib.pmlscript")
    script = residues([4, 3, 1])
    assert script.count("create ") == 2
    assert "create G7-1hckA" in script
    assert "create Q5-1hckA" in script
    assert "None" not in script
    assert caplog.messages == ["Column 3 is a deletion in chain A; no residue to show"]
