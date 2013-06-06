"""PyMOL script generator.

The script syntax is the command language used by PyMOL (based on Python's cmd
module), not using the pymol API directly. You generate the script, then you
load the script in PyMOL and the pretty happens.
"""

import os.path


# ---------------------------------------------------------
# Mapping pattern to PDB sequence alignment 

# New shit:
#   - for each PDB, build a mapping of overall alignment column numbers to
#     PDB resnums
#       - like biocma.utils.get_equivs; currently in cladecompare.py
#   - Also record locations of inserts
#       - biocma.utils.get_inserts
#       - cladecompare.py still botches this


# ---------------------------------------------------------
# Script builders

def build_single(mapping, inserts, patterns, pdbfname, chain):
    """Build a PyMOL script showing a single structure with selected residues."""
    pml_name = get_pml_name(pdbfname, chain)
    outs = [mk_intro_one(pdbfname, pml_name, chain),
            HR,
            mk_struct(pml_name, 'MainStruct_'+pml_name, chain=chain,
                      color='smudge', transparency=0.3),
            HR]
    if inserts:
        outs.extend(make_inserts(inserts, 'MainStruct_'+pml_name,
                                 chain, 'gray50'))
    if patterns:
        outs.extend(make_residues(mapping, patterns, pml_name, chain))
    outs.append(mk_outro())
    return '\n'.join(outs)


def build_multi(mapping, inserts, key_residues, pdbfnames, chains):
    """Superimpose multiple structures onto a reference, showing equivalent
    selected residues in each.

    To reduce clutter, only show residues deviating from the reference side
    chain by at least `threshold` Angstroms RMS.
    """
    # TODO - ensure Pymol's automatic struct colors aren't clobbered
    pml_names = [get_pml_name(pfn, cid) for pfn, cid in zip(pdbfnames, chains)]
    ref_pml_name, ref_chn = (pml_names[0], chains[0])

    outs = [mk_intro_multi(pdbfnames, pml_names, chains),
            HR,
            mk_struct(ref_pml_name, 'RefStruct_'+ref_pml_name,
                      chain=ref_chn, color='smudge', transparency=0.7)]
    if inserts:
        outs.extend(make_inserts(inserts, 
                                 'RefStruct_'+ref_pml_name, ref_chn,
                                 'gray70'))
    if key_residues:
        # Side chains for the reference PDB
        outs.extend(make_residues(mapping, key_residues, ref_pml_name, ref_chn))
    for eqv_pml_name, eqv_chn in zip(pml_names[1:], chains[1:]):
        outs.append(mk_struct(eqv_pml_name, 'EqvStruct_'+eqv_pml_name,
                              chain=eqv_chn, color='slate',
                              transparency=0.7))
        # Side chains for the other PDBs
        if inserts:
            outs.extend(make_inserts(inserts, 
                                     'EqvStruct_'+eqv_pml_name, eqv_chn,
                                     'marine'))
        if key_residues:
            # Generate PyMOL script lines
            outs.extend(make_residues(mapping, key_residues,
                                      eqv_pml_name,
                                      eqv_chn))

    outs.extend([HR, mk_outro()])

    return '\n'.join(outs)      # just the script


def make_inserts(inserts, pml_name, chain, color):
    """Highlight the insert regions of a structure in PyMOL.

    Inserts are according to the CMA file -- lowercase letters.
    """
    for start, end in inserts:
        if end - start > 1:
            yield mk_range(pml_name, start, end, chain, color)


def make_residues(mapping, patterns, pml_name, chain, pale=False):
    """Build residues in PyMOL for the simple cases."""
    # for color, posns in colorize(patterns, pale):
        # for posn in posns:
    for posn in patterns:
        restype, resid = mapping[posn]
        yield mk_residue(pml_name, resid, restype, chain, mol=pml_name,
                         # color=color)
                         color='magenta')



# ---------------------------------------------------------
# Helper functions

def colorize(patterns, pale):
    """Map patterns to color schemes.

    Order of appearance of colors:
        magenta 2
        orange  3
        yellow  4
        green   5
        cyan    1
        blue    6
    """
    dark_colors = ['lightmagenta', 'brightorange', 'tv_yellow',
                   'tv_green', 'cyan', 'marine']
    pale_colors = ['lightpink', 'lightorange', 'paleyellow',
                   'palegreen', 'palecyan', 'lightblue']
    colorset = pale_colors if pale else dark_colors
    # assert patterns, 'No pattern files given (*.pttrn, *.ptn)'
    if not patterns:
        return
    if len(patterns) == 1:
        # Usual scBPPS result
        return [(colorset[4], patterns[0])]
    if len(patterns) == len(colorset):
        # Maxed-out mcBPPS result
        return zip(colorset, patterns)
    # Follow the spectrum but make sure the last category is shown in cyan
    return zip(colorset, patterns[:-1]) + [(colorset[-2], patterns[-1])]


def get_pml_name(pdb_fname, chain_id):
    pdb_id, chn = pdbid_from_fname(pdb_fname)
    if chain_id is None:
        chain_id = chn
    if chain_id is None:
        return "%s" % pdb_id.lower()
    return "%s%s" % (pdb_id.lower(), chain_id.upper())


def pdbid_from_fname(fname, upper=False, default_chain=None):
    """Extract an uppercase PDB ID and chain ID from a filename.

    Ex:
        1ABC.pdb
        pdb1ABC.ent
        00_1abc.pdb, 01_1abc.pdb
        1abcA.pdb
        1abc_A.pdb
        1abc_chain_A.pdb (?)

    will all return ('1ABC', None) or ('1ABC', 'A').
    """
    fn = os.path.basename(fname)
    if fn.endswith('.ent'):
        assert len(fn) == 11 and fn.startswith('pdb')
        pdbid = fn[3:-4]
        chnid = fn[3:-4].upper()
    else:
        assert fn.endswith('.pdb')
        if fn.startswith('01_') or fn.startswith('00_'):
            fn = fn[3:]
        pdbid = fn[:4]
        if len(fn) == 8:
            assert '_' not in fn
            chnid = default_chain
        elif len(fn) == 9:
            assert '_' not in fn
            chnid = fn[4].upper()
        elif len(fn) == 10:
            assert fn[4] == '_'
            chnid = fn[5].upper()
        else:
            assert 'chain' in fn
            # TODO - where's the chain ID?
            chnid = default_chain
    return (pdbid.upper() if upper else pdbid, chnid)



# ---------------------------------------------------------
# Script templates 

HR = '\n# %s\n' % ('-' * 70)


def mk_intro_one(pdb, name=None, description=None):
    """Generate the initializing code for a PyMol script."""
    return """\
# PyMOL script -- %(description)s

# Load the pdb file
load %(pdb)s, %(name)s
hide all
""" % {
        'pdb': pdb,
        'name': name or pdbid_from_fname(pdb)[0],
        'description': description or name or pdbid_from_fname(pdb)[0],
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
            % (len(pdbs), len(names)))
    return """\
# PyMOL script -- %(description)s

%(load)s
%(align)s
hide all
""" % {
    'description': description or names[0],
    'load':     '\n'.join("load %s, %s" % (p, n)
                          for p, n in zip(pdbs, names)),
    'align':    ('\n'.join("align (%s and chain %s), (%s and chain %s)"
                           % (oname, ochain, names[0], chains[0])
                           for oname, ochain in zip(names[1:], chains[1:]))
                 if chains and len(chains) == len(filter(bool, chains))
                 else '\n'.join("cealign %s, %s" % (names[0], other)
                                for other in names[1:]))
} 


def mk_struct(name, label, resn=None, chain=None,
              color='white', style='cartoon',
              transparency=None):
    """Generate the PyMol code to display a structure."""
    return """
create %(label)s, %(name)s %(chain)s %(resn)s
show %(style)s, %(label)s
color %(color)s, %(label)s 
%(transparency)s
""" % {
        'label':    label,  # ENH: make optional
        'name':     name,
        'resn':     ("and resn %s" % resn) if resn else '',
        'chain':    ("and chain %s" % chain) if chain else '',
        'color':    color,
        'style':    style,
        'transparency': ("set %s_transparency=%s, %s\n" %
                    (style.rstrip('s'), transparency, label) 
                    if transparency is not None else ''),
        }


def mk_range(name, startnum,  endnum, chain=None, color='gray70'):
    """Generate the PyMol code to display a colorized region."""
    return """
color %(color)s, %(name)s and resi %(snum)s-%(enum)s %(chain)s 
""" % {
        'name':     name,
        'snum':     startnum,
        'enum':     endnum,
        'chain':    ('and chain %s' % chain) if chain else '',
        'color':    color,
        }


def mk_residue(name, num,  restype='X', chain=None, label=None, mol='',
               color='magenta', style=None):
    """Generate the PyMol code to display a residue."""
    return """
create %(label)s, %(name)s and resi %(num)s %(chain)s and not (name %(bbone)s)
show %(style)s, %(label)s
color %(color)s, %(label)s
util.cnc %(label)s
""" % {
        'label':    label or ("%s%s-%s" % (restype, num, mol)),
        'name':     name,
        'chain':    ('and chain %s' % chain) if chain else '',
        'bbone':    "c,n,o" if restype != 'P' else "c,o",
        'num':      num,
        'color':    color,
        'style':    style or ('sticks' if restype != 'G' else 'spheres'),
        }


def mk_outro(background='black'):
    """Generate the finishing code for a PyMol script."""
    return """
# Viewing options
bg_color black
reset
"""

