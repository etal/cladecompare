#!/usr/bin/env python3
"""Write the inputs for the Python 2.7 reference capture into inputs/.

Run from any directory with Python 3.6 or later; it reads the test files
under test/ in this repository and downloads two structures from the RCSB
Protein Data Bank.
Rerunning it reproduces the committed files byte for byte, as long as the
RCSB has not remediated those two entries since.

Files written, all under tests/data/inputs/:

* scttlre-domain.fasta_aln.cma, cdc2-domain.fasta_aln.cma, scttlre-domain.fasta,
  cdc2-domain.fasta, CDK.hmm -- copied unchanged from test/.
* scttlre-domain.aln.fasta, cdc2-domain.aln.fasta -- the two CMA alignments
  reduced to their consensus columns as cladecomparelib.core.cma_blocks does:
  keep the aligned part of each sequence line and drop its lowercase insert
  residues.  Every sequence line in these two files has the plain form
  {()SEQ()}*, so the fixed slice below is equivalent.  cladeweb
  reads uploads as FASTA, and these are the same alignments in that form.
* cdc2-a.aln.fasta, cdc2-b.aln.fasta -- cdc2-domain.aln.fasta split into its
  first CDC2_SPLIT records and the rest, for the three-alignment mode.
* cdc2-trimmed.aln.fasta -- cdc2-domain.aln.fasta with columns TRIM_START to
  TRIM_END (1-based, inclusive) removed, so that it is narrower than the
  scttlre alignment and cladecompare must align the two with MUSCLE.
* 1HCK.pdb -- CDK2 bound to ATP, one chain, unchanged from the RCSB.
* 1FIN-AB-CA.pdb -- the CDK2/cyclin A complex 1FIN reduced to its HEADER
  line and the C-alpha ATOM records of chains A (CDK2) and B (cyclin A), with
  TER and END records.  It exercises a structure with more than one chain.
"""

import shutil
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
SRC = REPO / "test"
OUT = HERE / "inputs"

COPIED = [
    "scttlre-domain.fasta_aln.cma",
    "cdc2-domain.fasta_aln.cma",
    "scttlre-domain.fasta",
    "cdc2-domain.fasta",
    "CDK.hmm",
]

CDC2_SPLIT = 18
TRIM_START, TRIM_END = 121, 130
RCSB = "https://files.rcsb.org/download/{}.pdb"


def cma_blocks(path):
    """Read a CMA file of plain {()SEQ()}* lines as (id, sequence) pairs.

    On such lines this gives the same result as core.cma_blocks.
    """
    records = []
    lines = iter(path.read_text().splitlines())
    for line in lines:
        if line.startswith(">"):
            acc = line.split(None, 1)[0][1:]
            seq = next(lines).strip()
            records.append((acc, "".join(c for c in seq[3:-4] if not c.islower())))
    return records


def write_text(path, text):
    with open(path, "w", newline="\n") as handle:
        handle.write(text)


def write_fasta(path, records):
    write_text(path, "".join(f">{acc}\n{seq}\n" for acc, seq in records))


def fetch(pdb_id):
    with urllib.request.urlopen(RCSB.format(pdb_id)) as response:
        return response.read().decode("ascii")


def reduce_pdb(text, chains):
    keep = []
    n_atoms = 0
    for line in text.splitlines():
        if line.startswith("HEADER"):
            keep.append(line)
        elif (
            line.startswith("ATOM")
            and line[12:16].strip() == "CA"
            and line[21] in chains
        ):
            keep.append(line)
            n_atoms += 1
        elif line.startswith("TER") and len(line) > 21 and line[21] in chains:
            keep.append(line)
    assert n_atoms
    keep.append("END")
    return "\n".join(keep) + "\n"


def main():
    OUT.mkdir(exist_ok=True)
    for name in COPIED:
        shutil.copyfile(SRC / name, OUT / name)

    scttlre = cma_blocks(SRC / "scttlre-domain.fasta_aln.cma")
    cdc2 = cma_blocks(SRC / "cdc2-domain.fasta_aln.cma")
    write_fasta(OUT / "scttlre-domain.aln.fasta", scttlre)
    write_fasta(OUT / "cdc2-domain.aln.fasta", cdc2)
    write_fasta(OUT / "cdc2-a.aln.fasta", cdc2[:CDC2_SPLIT])
    write_fasta(OUT / "cdc2-b.aln.fasta", cdc2[CDC2_SPLIT:])
    write_fasta(
        OUT / "cdc2-trimmed.aln.fasta",
        [(acc, seq[: TRIM_START - 1] + seq[TRIM_END:]) for acc, seq in cdc2],
    )

    write_text(OUT / "1HCK.pdb", fetch("1HCK"))
    write_text(OUT / "1FIN-AB-CA.pdb", reduce_pdb(fetch("1FIN"), "AB"))


if __name__ == "__main__":
    main()
