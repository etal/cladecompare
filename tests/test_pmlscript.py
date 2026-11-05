"""Tests of pmlscript: PyMOL object names and the residues drawn."""

import logging

import pytest

from cladecompare import pmlscript


@pytest.mark.parametrize(
    "fname, idcode, expected",
    [
        # The PDB ID of the HEADER record wins over the file name.
        ("1FIN-AB-CA.pdb", "1FIN", "1finA"),
        ("pdb1hck.ent", "1HCK", "1hckA"),
        # Without one, the file name gives the name.
        ("dir/AF-P24941-F1-model_v6.pdb", "", "af_p24941_f1_model_v6A"),
        ("1FIN-AB-CA.pdb", "", "1fin_ab_caA"),
        ("model.v2.pdb", "", "model_v2A"),
    ],
)
def test_get_pml_name(fname, idcode, expected):
    assert pmlscript.get_pml_name(fname, "a", idcode) == expected


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
    caplog.set_level(logging.WARNING, logger="cladecompare.pmlscript")
    script = residues([4, 3, 1])
    assert script.count("create ") == 2
    assert "create G7-1hckA" in script
    assert "create Q5-1hckA" in script
    assert "None" not in script
    assert caplog.messages == ["Column 3 is a deletion in chain A; no residue to show"]
