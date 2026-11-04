"""PyMOL script generator.

The script syntax is the command language used by PyMOL (based on Python's cmd
module), not using the pymol API directly. You generate the script, then you
load the script in PyMOL and the pretty happens.
"""

import logging
import os.path
import re

logger = logging.getLogger(__name__)

# ---------------------------------------------------------
# Mapping pattern to PDB sequence alignment

# The builders below take a mapping of alignment columns to PDB residue
# numbers and a list of insert ranges, both computed by
# core.aln_resnums_inserts. build_multi receives a tuple of each, one per
# structure, and passes the tuples on unchanged, so it fails with more than
# one structure.


# ---------------------------------------------------------
# Script builders


def build_single(mapping, inserts, patterns, pdbfname, chain, idcode):
    """Build a PyMOL script showing a single structure with selected residues."""
    pml_name = get_pml_name(pdbfname, chain, idcode)
    outs = [
        mk_intro_one(pdbfname, pml_name, chain.strip()),
        HR,
        mk_struct(
            pml_name,
            "MainStruct_" + pml_name,
            chain=chain,
            color="smudge",
            transparency=0.3,
        ),
        HR,
    ]
    if inserts:
        outs.extend(make_inserts(inserts, "MainStruct_" + pml_name, chain, "gray50"))
    if patterns:
        outs.extend(make_residues(mapping, patterns, pml_name, chain))
    outs.append(mk_outro())
    return "\n".join(outs)


def build_multi(mapping, inserts, key_residues, pdbfnames, chains, idcodes):
    """Superimpose multiple structures onto a reference, showing equivalent
    selected residues in each.

    To reduce clutter, only show residues deviating from the reference side
    chain by at least `threshold` Angstroms RMS.
    """
    # TODO - ensure Pymol's automatic struct colors aren't clobbered
    pml_names = [
        get_pml_name(pfn, cid, idc)
        for pfn, cid, idc in zip(pdbfnames, chains, idcodes, strict=True)
    ]
    ref_pml_name, ref_chn = (pml_names[0], chains[0])

    outs = [
        mk_intro_multi(pdbfnames, pml_names, chains),
        HR,
        mk_struct(
            ref_pml_name,
            "RefStruct_" + ref_pml_name,
            chain=ref_chn,
            color="smudge",
            transparency=0.7,
        ),
    ]
    if inserts:
        outs.extend(
            make_inserts(inserts, "RefStruct_" + ref_pml_name, ref_chn, "gray70")
        )
    if key_residues:
        # Side chains for the reference PDB
        outs.extend(make_residues(mapping, key_residues, ref_pml_name, ref_chn))
    for eqv_pml_name, eqv_chn in zip(pml_names[1:], chains[1:], strict=True):
        outs.append(
            mk_struct(
                eqv_pml_name,
                "EqvStruct_" + eqv_pml_name,
                chain=eqv_chn,
                color="slate",
                transparency=0.7,
            )
        )
        # Side chains for the other PDBs
        if inserts:
            outs.extend(
                make_inserts(inserts, "EqvStruct_" + eqv_pml_name, eqv_chn, "marine")
            )
        if key_residues:
            # Generate PyMOL script lines
            outs.extend(make_residues(mapping, key_residues, eqv_pml_name, eqv_chn))

    outs.extend([HR, mk_outro()])

    return "\n".join(outs)  # just the script


def make_inserts(inserts, pml_name, chain, color):
    """Highlight the insert regions of a structure in PyMOL.

    Inserts are according to the CMA file -- lowercase letters.
    """
    for start, end in inserts:
        if end - start > 1:
            yield mk_range(pml_name, start, end, chain, color)


def make_residues(mapping, patterns, pml_name, chain):
    """Draw the structure's residue at each significant column.

    `mapping` is a list with one (residue type, residue number) pair per
    match state of the profile, in order, as core.aln_resnums_inserts
    returns it. Each site in `patterns` is an alignment column number counted
    from 1, so site N is entry N - 1. A site whose match state is a deletion
    in the structure has no residue to draw; it is logged and left out.
    """
    for posn in patterns:
        restype, resid = mapping[posn - 1]
        if resid is None:
            logger.warning(
                "Column %d is a deletion in chain %s; no residue to show",
                posn,
                chain_label(chain),
            )
            continue
        yield mk_residue(pml_name, resid, restype, chain, mol=pml_name, color="magenta")


# ---------------------------------------------------------
# Helper functions


def get_pml_name(pdb_fname, chain_id, idcode):
    """Name the PyMOL object of a structure.

    The name is the PDB ID that the file's HEADER record gives, or, where it
    gives none, as in a model from AlphaFold, the file name without its
    extension, followed by the chain ID. Characters other than ASCII letters,
    digits and '_' become '_', since PyMOL takes those three unchanged in an
    object name and in a selection.
    """
    base = idcode or os.path.splitext(os.path.basename(pdb_fname))[0]
    base = re.sub(r"\W", "_", base, flags=re.ASCII).lower()
    return base + chain_id.strip().upper()


def chain_label(chain):
    """A chain ID as a selection and --chain give it.

    PyMOL selects a chain whose ID is blank as chain "".
    """
    return chain.strip() or '""'


def chain_selection(chain):
    """The part of a PyMOL selection that restricts it to a chain."""
    if chain is None:
        return ""
    return "and chain %s" % chain_label(chain)


# ---------------------------------------------------------
# Script templates

HR = "\n# %s\n" % ("-" * 70)


def mk_intro_one(pdb, name, description=None):
    """Generate the initializing code for a PyMol script."""
    return """\
# PyMOL script -- %(description)s

# Load the pdb file
load %(pdb)s, %(name)s
hide all
""" % {
        "pdb": pdb,
        "name": name,
        "description": description or name,
    }


def mk_intro_multi(pdbs, names, chains=None, description=None):
    """Generate the initializing PyMol code for multiple aligned proteins."""
    if isinstance(pdbs, str):
        pdbs = [pdbs]
    if isinstance(names, str):
        names = [names]
    if not len(pdbs) == len(names):
        raise ValueError(
            "Number of PDB filenames (%s) and object names (%s) don't match!"
            % (len(pdbs), len(names))
        )
    return """\
# PyMOL script -- %(description)s

%(load)s
%(align)s
hide all
""" % {
        "description": description or names[0],
        "load": "\n".join(
            "load %s, %s" % (p, n) for p, n in zip(pdbs, names, strict=True)
        ),
        "align": (
            "\n".join(
                "align (%s %s), (%s %s)"
                % (oname, chain_selection(ochain), names[0], chain_selection(chains[0]))
                for oname, ochain in zip(names[1:], chains[1:], strict=True)
            )
            if chains and all(c is not None for c in chains)
            else "\n".join("cealign %s, %s" % (names[0], other) for other in names[1:])
        ),
    }


def mk_struct(
    name,
    label,
    resn=None,
    chain=None,
    color="white",
    style="cartoon",
    transparency=None,
):
    """Generate the PyMol code to display a structure."""
    return """
create %(label)s, %(name)s %(chain)s %(resn)s
show %(style)s, %(label)s
color %(color)s, %(label)s 
%(transparency)s
""" % {
        "label": label,  # ENH: make optional
        "name": name,
        "resn": ("and resn %s" % resn) if resn else "",
        "chain": chain_selection(chain),
        "color": color,
        "style": style,
        "transparency": (
            "set %s_transparency=%s, %s\n" % (style.rstrip("s"), transparency, label)
            if transparency is not None
            else ""
        ),
    }


def mk_range(name, startnum, endnum, chain=None, color="gray70"):
    """Generate the PyMol code to display a colorized region."""
    return """
color %(color)s, %(name)s and resi %(snum)s-%(enum)s %(chain)s 
""" % {
        "name": name,
        "snum": startnum,
        "enum": endnum,
        "chain": chain_selection(chain),
        "color": color,
    }


def mk_residue(
    name, num, restype="X", chain=None, label=None, mol="", color="magenta", style=None
):
    """Generate the PyMol code to display a residue."""
    return """
create %(label)s, %(name)s and resi %(num)s %(chain)s and not (name %(bbone)s)
show %(style)s, %(label)s
color %(color)s, %(label)s
util.cnc %(label)s
""" % {
        "label": label or ("%s%s-%s" % (restype, num, mol)),
        "name": name,
        "chain": chain_selection(chain),
        "bbone": "c,n,o" if restype != "P" else "c,o",
        "num": num,
        "color": color,
        "style": style or ("sticks" if restype != "G" else "spheres"),
    }


def mk_outro(background="black"):
    """Generate the finishing code for a PyMol script."""
    return """
# Viewing options
bg_color black
reset
"""
