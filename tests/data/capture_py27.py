#!/usr/bin/env python2.7
"""Capture reference outputs of the Python 2 CladeCompare code base.

This script must run under Python 2.7 with Biopython 1.76, SciPy, ReportLab,
bottle and WebLogo 3.6, with HMMER 3.0 and MUSCLE 3.8.1551 installed beside
that interpreter; see README.md beside it for the environment.  It does not
run the working tree.  It extracts the sources it runs from git:

* cladecompare at CLADECOMPARE_COMMIT, the last commit to change code before
  the reference was first captured,
* cladecompare at LEGACY_COMMIT, the last commit at which -s phospho ran
  (unweighted only), for the legacy cases,
* biofrills at BIOFRILLS_COMMIT, its last Python 2 commit, from the biofrills
  repository beside this one (or from $BIOFRILLS_REPO).

It deletes and rewrites tests/data/py27/:

* cli/<case>/ -- one command-line case: run.json holds the source tree it
  ran (a key of versions.json's "commits"), each command's arguments, return
  code, standard output and standard error, and the names of the files it
  wrote; the text files it wrote sit beside run.json.
* web/<case>/ -- one form submission to cladeweb through its WSGI
  application: run.json holds the request, the response status, the
  exception the application logged, if any, and its log output and standard
  output; response.html holds the body.
* versions.json -- the producing interpreter, packages, tools, commits and
  hash seed.

Failures of the original code are part of the reference and are recorded,
not raised.  Paths are normalized as described in README.md.  PYTHONHASHSEED
is passed through to every command and recorded; the committed reference
was captured with it unset.

Usage, from any directory:

    <python2.7> <repository>/tests/data/capture_py27.py
"""
from __future__ import print_function

import io
import json
import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import wsgiref.util
from collections import OrderedDict
from StringIO import StringIO

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
INPUTS = os.path.join(HERE, 'inputs')
OUT = os.path.join(HERE, 'py27')

CLADECOMPARE_COMMIT = 'e265b1f5321c0bcdd6ed05b9cce5ebb26483ef54'
LEGACY_COMMIT = '44c81cce57d13ef0cf4caaa8a0cb6fa371f3ced7'
BIOFRILLS_COMMIT = '36684bb6c7632f96215e8b2b4ebc86640f331bcd'
HMMER_VERSION = 'HMMER 3.0 '
MUSCLE_VERSION = 'MUSCLE v3.8.1551'

STRATEGIES = ['gtest', 'urn', 'jsd', 'hypg', 'phospho']
FG_CMA = 'scttlre-domain.fasta_aln.cma'
BG_CMA = 'cdc2-domain.fasta_aln.cma'

# The eight test/Makefile targets, keyed by (mode, strategy): the target's
# name, whether it passes -w, and its -a value (None for the default).
MAKEFILE = {
    ('pair', 'gtest'): ('scttlre-v-cdc2-gtest', True, None),
    ('pair', 'urn'): ('scttlre-v-cdc2-urn', False, None),
    ('pair', 'jsd'): ('scttlre-v-cdc2-jsd', False, '.2'),
    ('pair', 'hypg'): ('scttlre-v-cdc2-hypg', False, None),
    ('single', 'gtest'): ('cdc2-gtest', True, '.001'),
    ('single', 'urn'): ('cdc2-urn', False, None),
    ('single', 'jsd'): ('cdc2-jsd', False, '.2'),
    ('single', 'hypg'): ('cdc2-hypg', False, '.2'),
}


def git_extract(repo, commit, dest):
    """Extract the tree of `commit` in `repo` into `dest`; return full SHA."""
    sha = subprocess.check_output(
        ['git', '-C', repo, 'rev-parse', commit + '^{commit}']).strip()
    tar_bytes = subprocess.check_output(['git', '-C', repo, 'archive', sha])
    with tarfile.open(fileobj=io.BytesIO(tar_bytes)) as tar:
        tar.extractall(dest)
    return sha


class Normalizer(object):
    """Replace machine-specific paths in recorded text."""

    def __init__(self):
        self.pairs = []

    def add(self, path, label):
        for variant in set([path, os.path.realpath(path)]):
            self.pairs.append((variant, label))
        # Longest first, so that a directory inside another is replaced first.
        self.pairs.sort(key=lambda pair: -len(pair[0]))

    def __call__(self, text):
        for path, label in self.pairs:
            text = text.replace(path, label)
        # Random parts of mkstemp/mkdtemp names
        return re.sub(r'<tmp>/tmp[A-Za-z0-9_]{6,}', '<tmp>/tmpXXXXXX', text)


def check_tools():
    hmmalign = subprocess.check_output(['hmmalign', '-h']).splitlines()[1]
    muscle = subprocess.check_output(['muscle', '-version'],
                                     stderr=subprocess.STDOUT).strip()
    assert HMMER_VERSION in hmmalign, hmmalign
    assert MUSCLE_VERSION in muscle, muscle
    return hmmalign.lstrip('# '), muscle


def versions(shas, tools):
    import pkg_resources
    dists = OrderedDict()
    for name in ['biopython', 'numpy', 'scipy', 'reportlab', 'bottle',
                 'weblogo', 'Pillow']:
        dists[name] = pkg_resources.get_distribution(name).version
    info = OrderedDict()
    info['python'] = sys.version.split()[0]
    info['implementation'] = platform.python_implementation()
    info['platform'] = OrderedDict([('system', platform.system()),
                                    ('release', platform.release()),
                                    ('machine', platform.machine())])
    info['packages'] = dists
    info['tools'] = OrderedDict([('hmmalign', tools[0]), ('muscle', tools[1])])
    info['commits'] = shas
    info['hash_seed'] = os.environ.get('PYTHONHASHSEED')
    return info


# --- Command-line cases ---

def make_case(case, inputs, extra=(), source='cladecompare', report=True,
              target=None):
    """Return (case id, source key, list of argv, Makefile target or None).

    The commands are cladecompare writing <case>.out and <case>.pttrn and,
    unless `report` is false, cladereport writing <case>.html.
    """
    commands = [['cladecompare.py'] + inputs + list(extra)
                + ['-p', case + '.pttrn', '-o', case + '.out']]
    if report:
        commands.append(['cladereport.py', case + '.out',
                         '-p', case + '.pttrn', '-o', case + '.html'])
    return (case, source, commands, target)


def cli_cases():
    """Yield (case id, source key, list of argv, Makefile target or None)."""
    for mode in ['pair', 'single']:
        inputs = [FG_CMA, BG_CMA] if mode == 'pair' else [BG_CMA]
        for strategy in STRATEGIES:
            target, mk_weight, alpha = MAKEFILE.get(
                (mode, strategy), (None, False, None))
            for weight in [False, True]:
                case = '%s.%s%s%s' % (mode, strategy, '.w' if weight else '',
                                      '.a' + alpha if alpha else '')
                extra = ['-s', strategy]
                if weight:
                    extra.append('-w')
                if alpha:
                    extra += ['-a', alpha]
                yield make_case(
                    case, inputs, extra,
                    target=target if target and weight == mk_weight else None)

    # Three alignments: each in turn is the foreground against the other two.
    # The output names follow the input names, so no -o or -p is given.
    multi_in = [FG_CMA, 'cdc2-a.aln.fasta', 'cdc2-b.aln.fasta']
    yield ('multi.gtest', 'cladecompare',
           [['cladecompare.py'] + multi_in,
            ['cladereport.py'] + [name + '.out' for name in multi_in]
            + ['-o', 'multi.gtest.html']], None)

    unaligned = ['scttlre-domain.fasta', 'cdc2-domain.fasta']
    # Alignments of unequal width, aligned to each other by MUSCLE.
    yield make_case('unequal.gtest', ['scttlre-domain.aln.fasta',
                                      'cdc2-trimmed.aln.fasta'])
    # Unaligned sequences aligned by hmmalign.
    yield make_case('hmm.pair.gtest', unaligned, ['--hmm', 'CDK.hmm'])
    yield make_case('hmm.single.gtest', ['cdc2-domain.fasta'],
                    ['--hmm', 'CDK.hmm'])
    # Structures aligned by hmmalign; the PyMOL script is the output.
    yield make_case('pdb.1HCK.gtest', unaligned,
                    ['--hmm', 'CDK.hmm', '--pdb', '1HCK.pdb'])
    yield make_case('pdb.1FIN-AB-CA.gtest', unaligned,
                    ['--hmm', 'CDK.hmm', '--pdb', '1FIN-AB-CA.pdb'])

    # phospho at the last commit where it ran: unweighted only (-nw), with
    # the current default alpha given explicitly.  cladereport at that commit
    # has no -o option; the report is not the point of these cases.
    for mode in ['pair', 'single']:
        inputs = [FG_CMA, BG_CMA] if mode == 'pair' else [BG_CMA]
        yield make_case('legacy.%s.phospho.nw' % mode, inputs,
                        ['-s', 'phospho', '-nw', '-a', '0.005'],
                        source='cladecompare_legacy', report=False)


def run_cli_case(case, source, source_dir, commands, target, env, scratch,
                 norm):
    work = os.path.join(scratch, 'work', case)
    tmp = os.path.join(scratch, 'tmp', case)
    shutil.copytree(INPUTS, work)
    os.makedirs(tmp)
    inputs = set(os.listdir(work))
    env = dict(env, TMPDIR=tmp)
    norm.add(work, '<workdir>')
    norm.add(tmp, '<tmp>')

    runs = []
    for argv in commands:
        proc = subprocess.Popen(
            [sys.executable, os.path.join(source_dir, argv[0])] + argv[1:],
            cwd=work, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        runs.append(OrderedDict([('argv', argv),
                                 ('returncode', proc.returncode),
                                 ('stdout', norm(stdout)),
                                 ('stderr', norm(stderr))]))
        if proc.returncode:
            break

    dest = os.path.join(OUT, 'cli', case)
    os.makedirs(dest)
    written = sorted(set(os.listdir(work)) - inputs)
    pdfs = [name for name in written if name.endswith('.pdf')]
    texts = [name for name in written if not name.endswith('.pdf')]
    for name in texts:
        shutil.copyfile(os.path.join(work, name), os.path.join(dest, name))
    record = OrderedDict([
        ('source', source),
        ('makefile_target', target),
        ('commands', runs),
        ('written_files', texts),
        # Pair-logo PDFs embed creation times, so only their names are kept.
        ('pdfs', pdfs),
        ('leftover_temp_files', len(os.listdir(tmp))),
    ])
    write_json(os.path.join(dest, 'run.json'), record)


# --- Web cases ---

BOUNDARY = 'capture-py27-boundary'


def multipart(fields, files):
    parts = []
    for name, value in fields:
        parts.append('--%s\r\nContent-Disposition: form-data; name="%s"'
                     '\r\n\r\n%s\r\n' % (BOUNDARY, name, value))
    for name, fname in files:
        with open(os.path.join(INPUTS, fname), 'rb') as handle:
            data = handle.read()
        parts.append('--%s\r\nContent-Disposition: form-data; name="%s"; '
                     'filename="%s"\r\nContent-Type: application/octet-stream'
                     '\r\n\r\n%s\r\n' % (BOUNDARY, name, fname, data))
    return ''.join(parts) + '--%s--\r\n' % BOUNDARY


def web_cases():
    """Yield (case id, form fields, uploaded files)."""
    pair = [('seqfile1', 'scttlre-domain.aln.fasta'),
            ('seqfile2', 'cdc2-domain.aln.fasta')]
    single = [('seqfile1', 'cdc2-domain.aln.fasta')]
    # hypg is accepted by the command line but not offered by the form.
    for strategy in STRATEGIES:
        yield ('web.pair.' + strategy,
               [('name', 'web-pair-' + strategy), ('strategy', strategy),
                ('alpha', '0.005')], pair)
    for strategy in ['gtest', 'urn', 'jsd', 'phospho']:
        yield ('web.single.' + strategy,
               [('name', 'web-single-' + strategy), ('strategy', strategy),
                ('alpha', '0.005')], single)
    yield ('web.pair-cma.gtest',
           [('name', 'web-pair-cma'), ('strategy', 'gtest'),
            ('alpha', '0.005')],
           [('seqfile1', FG_CMA), ('seqfile2', BG_CMA)])
    yield ('web.hmm.pair.gtest',
           [('name', 'web-hmm'), ('strategy', 'gtest'), ('alpha', '0.005')],
           [('seqfile1', 'scttlre-domain.fasta'),
            ('seqfile2', 'cdc2-domain.fasta'), ('profile', 'CDK.hmm')])
    yield ('web.pdb.1HCK.gtest',
           [('name', 'web-pdb'), ('strategy', 'gtest'), ('alpha', '0.005')],
           [('seqfile1', 'scttlre-domain.fasta'),
            ('seqfile2', 'cdc2-domain.fasta'), ('profile', 'CDK.hmm'),
            ('pdbfile', '1HCK.pdb')])


def run_web_case(app, case, fields, files, scratch, norm):
    work = os.path.join(scratch, 'work', case)
    tmp = os.path.join(scratch, 'tmp', case)
    os.makedirs(work)
    os.makedirs(tmp)
    norm.add(work, '<workdir>')
    norm.add(tmp, '<tmp>')

    body = multipart(fields, files)
    environ = {}
    wsgiref.util.setup_testing_defaults(environ)
    errors = StringIO()
    environ.update({
        'REQUEST_METHOD': 'POST',
        'PATH_INFO': '/cladecompare',
        'CONTENT_TYPE': 'multipart/form-data; boundary=' + BOUNDARY,
        'CONTENT_LENGTH': str(len(body)),
        'wsgi.input': io.BytesIO(body),
        'wsgi.errors': errors,
    })
    status = []

    def start_response(stat, headers, exc_info=None):
        status.append(stat)

    # Log as cladeweb.py configures logging when run as a server.
    log = StringIO()
    handler = logging.StreamHandler(log)
    handler.setFormatter(
        logging.Formatter('%(module)s [@%(lineno)s]: %(message)s'))
    root = logging.getLogger()
    old_level = root.level
    root.addHandler(handler)
    root.setLevel(logging.INFO)
    old_cwd, old_tempdir, old_stdout = os.getcwd(), tempfile.tempdir, sys.stdout
    os.chdir(work)
    tempfile.tempdir = tmp
    # cladeweb prints one line per temporary file it removes.
    stdout = sys.stdout = StringIO()
    try:
        response = ''.join(app(environ, start_response))
    finally:
        sys.stdout = old_stdout
        tempfile.tempdir = old_tempdir
        os.chdir(old_cwd)
        root.removeHandler(handler)
        root.setLevel(old_level)

    traceback = errors.getvalue().strip()
    dest = os.path.join(OUT, 'web', case)
    os.makedirs(dest)
    with open(os.path.join(dest, 'response.html'), 'wb') as handle:
        handle.write(norm(response))
    record = OrderedDict([
        ('fields', fields),
        ('uploads', files),
        ('status', status[0]),
        ('error', norm(traceback) if traceback else None),
        ('log', norm(log.getvalue())),
        ('stdout', norm(stdout.getvalue())),
        ('written_files', sorted(os.listdir(work))),
        ('leftover_temp_files', len(os.listdir(tmp))),
    ])
    write_json(os.path.join(dest, 'run.json'), record)


def write_json(path, obj):
    with open(path, 'wb') as handle:
        json.dump(obj, handle, indent=2, separators=(',', ': '))
        handle.write('\n')


def main():
    assert sys.version_info[:2] == (2, 7), sys.version
    biofrills_repo = os.environ.get(
        'BIOFRILLS_REPO', os.path.join(os.path.dirname(REPO), 'biofrills'))

    os.environ['PATH'] = (os.path.dirname(sys.executable) + os.pathsep
                          + os.environ['PATH'])
    tools = check_tools()

    scratch = tempfile.mkdtemp(prefix='cc-capture-')
    try:
        sources = {
            'cladecompare': os.path.join(scratch, 'cladecompare'),
            'cladecompare_legacy': os.path.join(scratch, 'cladecompare-legacy'),
        }
        biofrills = os.path.join(scratch, 'biofrills')
        shas = OrderedDict()
        shas['cladecompare'] = git_extract(REPO, CLADECOMPARE_COMMIT,
                                           sources['cladecompare'])
        shas['cladecompare_legacy'] = git_extract(
            REPO, LEGACY_COMMIT, sources['cladecompare_legacy'])
        shas['biofrills'] = git_extract(biofrills_repo, BIOFRILLS_COMMIT,
                                        biofrills)

        norm = Normalizer()
        norm.add(sources['cladecompare'], '<cladecompare>')
        norm.add(sources['cladecompare_legacy'], '<cladecompare-legacy>')
        norm.add(biofrills, '<biofrills>')
        norm.add(sys.prefix, '<python-prefix>')

        if os.path.isdir(OUT):
            shutil.rmtree(OUT)
        os.makedirs(OUT)
        write_json(os.path.join(OUT, 'versions.json'), versions(shas, tools))

        for case, source, commands, target in cli_cases():
            src = sources[source]
            env = dict(os.environ,
                       PYTHONPATH=os.pathsep.join([src, biofrills]))
            run_cli_case(case, source, src, commands, target, env, scratch,
                         norm)
            print('cli', case)

        sys.path[:0] = [sources['cladecompare'], biofrills]
        import bottle
        import cladeweb  # registers its routes on bottle's default app
        assert cladeweb.__file__.startswith(sources['cladecompare'])
        app = bottle.default_app()
        for case, fields, files in web_cases():
            run_web_case(app, case, fields, files, scratch, norm)
            print('web', case)
    finally:
        shutil.rmtree(scratch)


if __name__ == '__main__':
    main()
