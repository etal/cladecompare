"""Tests of the PDB file name parsing in pmlscript."""

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
