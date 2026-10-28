"""PyMOL script generator.

The script syntax is the command language used by PyMOL (based on Python's cmd
module), not using the pymol API directly. You generate the script, then you
load the script in PyMOL and the pretty happens.
"""

import os.path

from . import InputError

# ---------------------------------------------------------
# Mapping pattern to PDB sequence alignment

# The builders below take a mapping of alignment columns to PDB residue
# numbers and a list of insert ranges, both computed by
# core.aln_resnums_inserts.


# ---------------------------------------------------------
# Script builders


def build_single(mapping, inserts, patterns, pdbfname, chain):
    """Build a PyMOL script showing a single structure with selected residues."""
    pml_name = get_pml_name(pdbfname, chain)
    outs = [
        mk_intro_one(pdbfname, pml_name, chain),
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


def build_multi(mapping, inserts, key_residues, pdbfnames, chains):
    """Superimpose multiple structures onto a reference, showing equivalent
    selected residues in each.

    To reduce clutter, only show residues deviating from the reference side
    chain by at least `threshold` Angstroms RMS.
    """
    # TODO - ensure Pymol's automatic struct colors aren't clobbered
    pml_names = [
        get_pml_name(pfn, cid) for pfn, cid in zip(pdbfnames, chains, strict=True)
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
    """Build residues in PyMOL for the simple cases."""
    for posn in patterns:
        restype, resid = mapping[posn]
        yield mk_residue(pml_name, resid, restype, chain, mol=pml_name, color="magenta")


# ---------------------------------------------------------
# Helper functions


def get_pml_name(pdb_fname, chain_id):
    pdb_id, chn = pdbid_from_fname(pdb_fname)
    if chain_id is None:
        chain_id = chn
    if chain_id is None:
        return "%s" % pdb_id.lower()
    return "%s%s" % (pdb_id.lower(), chain_id.upper())


def pdbid_from_fname(fname):
    """Extract a PDB ID and chain ID from a filename.

    Ex:
        1abc.pdb, 00_1abc.pdb, 01_1abc.pdb, 1abc_chain_A.pdb -> ('1abc', None)
        1abcA.pdb, 1abc_a.pdb -> ('1abc', 'A')
        pdb1abc.ent -> ('1abc', '1ABC')

    The PDB ID keeps the case of the file name, and the chain ID is
    uppercase. Any other name raises InputError.
    """
    fn = os.path.basename(fname)
    if fn.endswith(".ent"):
        recognized = len(fn) == 11 and fn.startswith("pdb")
        pdbid = fn[3:-4]
        chnid = fn[3:-4].upper()
    elif fn.endswith(".pdb"):
        if fn.startswith(("01_", "00_")):
            fn = fn[3:]
        pdbid = fn[:4]
        if len(fn) == 8:
            recognized = "_" not in fn
            chnid = None
        elif len(fn) == 9:
            recognized = "_" not in fn
            chnid = fn[4].upper()
        elif len(fn) == 10:
            recognized = fn[4] == "_"
            chnid = fn[5].upper()
        else:
            recognized = "chain" in fn
            # TODO - where's the chain ID?
            chnid = None
    else:
        recognized = False
    if not recognized:
        raise InputError(
            "Cannot read a PDB ID from the file name %s; name the file as, "
            "for example, 1abc.pdb, 1abcA.pdb, 1abc_A.pdb or pdb1abc.ent" % fname
        )
    return pdbid, chnid


# ---------------------------------------------------------
# Script templates

HR = "\n# %s\n" % ("-" * 70)


def mk_intro_one(pdb, name=None, description=None):
    """Generate the initializing code for a PyMol script."""
    return """\
# PyMOL script -- %(description)s

# Load the pdb file
load %(pdb)s, %(name)s
hide all
""" % {
        "pdb": pdb,
        "name": name or pdbid_from_fname(pdb)[0],
        "description": description or name or pdbid_from_fname(pdb)[0],
    }


def mk_intro_multi(pdbs, names=None, chains=None, description=None):
    """Generate the initializing PyMol code for multiple aligned proteins."""
    if isinstance(pdbs, str):
        pdbs = [pdbs]
    if not names:
        names = [pdbid_from_fname(p)[0] for p in pdbs]
    elif isinstance(names, str):
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
                "align (%s and chain %s), (%s and chain %s)"
                % (oname, ochain, names[0], chains[0])
                for oname, ochain in zip(names[1:], chains[1:], strict=True)
            )
            if chains and len(chains) == len([c for c in chains if c])
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
        "chain": ("and chain %s" % chain) if chain else "",
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
        "chain": ("and chain %s" % chain) if chain else "",
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
        "chain": ("and chain %s" % chain) if chain else "",
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
