#!/usr/bin/python

import argparse
import contextlib
import os.path
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO

BOOTSTRAP_MODES = 'a',

# Some utilities
@contextlib.contextmanager
def sequences_in_format(sequences, fmt='fasta', **kwargs):
    with tempfile.NamedTemporaryFile(**kwargs) as tf:
        SeqIO.write(sequences, tf, fmt)
        tf.flush()
        yield tf.name

@contextlib.contextmanager
def temp_dir(**kwargs):
    """Maintains a temporary directory for the life of the context manager."""
    temp_dir = tempfile.mkdtemp(**kwargs)
    try:
        yield temp_dir
    finally:
        # Cleanup
        # ermm... this is breaking something (maybe bootstrapping replicates?), so leaving out for now
        #shutil.rmtree(temp_dir)
        pass

def stripext(f, basename=False):
    if basename:
        return stripext(os.path.basename(f))
    return os.path.splitext(f)[0]

def nonextant_file(path):
    if os.path.exists(path):
        raise ValueError("Exists: " + path)

    return path

def joiner(base_path):
    def p(*args):
        return os.path.join(base_path, *args)
    return p

def raxml(sequences, output_tree, stats_path=None, log_path=None, quiet=False,
        executable='raxmlHPC-SSE3', model='GTRGAMMA', threads=None,
        rapid_bootstrap=None, bootstrap_seed=None, tmp_prefix=None, outgroup=None):
    name = os.path.basename(os.path.splitext(output_tree)[0])

    def f(n):
        "Gets the RAxML file name associated with a key"
        return 'RAxML_{1}.{0}'.format(name, n)

    with temp_dir(prefix='raxml-') as td:
        with sequences_in_format(sequences, fmt='phylip-relaxed',
                prefix=tmp_prefix, dir=td) as seq_file:
            p = joiner(td)

            # note: -p is needed for some reason now but didn't use to be?
            cmd = [executable, '-n', name, '-m', model, '-s', seq_file, '-p', '9988']
            if threads and threads > 1:
                cmd.extend(('-T', str(threads)))

            if rapid_bootstrap:
                cmd.extend(('-f', 'a', '-x', bootstrap_seed,
                    '-N', rapid_bootstrap))

            if outgroup:
                cmd.extend(('-o', outgroup))

            stdout = stderr = None
            if quiet:
                stdout = stderr = open(os.path.devnull)

            cmd = map(str, cmd)

            print >> sys.stderr, "Running:", ' '.join(cmd)
            try:
                subprocess.check_call(cmd, stdout=stdout, stderr=stderr, cwd=td)
            except subprocess.CalledProcessError, e:
                raise SystemExit(e.returncode)

            # Get the result - either bootstrap-annotated tree or result
            key = 'bipartitions' if rapid_bootstrap else 'result'
            shutil.move(p(f(key)), output_tree)

            if stats_path:
                shutil.move(p(f('info')), stats_path)
            if log_path:
                shutil.move(p(f('log')), log_path)

def main():
    parser = argparse.ArgumentParser(description="""Simple wrapper around
    RAxML.  Abstracts executable selection and sequence formatting; only keeps
    desired files; name specification. Most arguments are *not* supported""")
    parser.add_argument('alignment_file', type=argparse.FileType('r'),
            help="""Input alignment""")
    parser.add_argument('--input-format', default='fasta',
            help="""Format of input file [default: %(default)s]""")
    parser.add_argument('output_tree', type=nonextant_file, help="""Destination
            for output tree""")
    parser.add_argument('--stats', type=nonextant_file, metavar="<stats file>",
            help="""Save RAxML stats to <stats file>""")
    parser.add_argument('--log', type=nonextant_file, metavar="<log file>",
            help="""Write RAxML log file to <log file>""")
    parser.add_argument('-q', '--quiet', action='store_true',
            help="""Suppress output""")

    bs_group = parser.add_argument_group("Bootstrap Options")
    bs_group.add_argument('--rapid-bootstrap', metavar='N',
            help="""Run rapid bootstrap analysis with N replicates""",
            type=int)
    bs_group.add_argument('-x', '--bootstrap-seed', help="""Bootstrap seed""",
            dest='bootstrap_seed', type=int, default=1)

    rax_group = parser.add_argument_group(title="""RAxML options""")
    rax_group.add_argument('-T', '--threads', help="""Number of
            threads to use [default: 1]""", type=int)
    rax_group.add_argument('--executable', help="""RAxML executable to use.
            [default: raxmlHPC-PTHREADS-SSE3 if threads > 1, raxmlHPC-SSE3
            otherwise]""")
    rax_group.add_argument('-m', '--model', default='GTRGAMMA', help="""RAxML
            model to use [default: %(default)s]""")
    parser.add_argument('-o', '--outgroup',
            help="""Fix output for tree""")

    args = parser.parse_args()
    if not args.executable:
        args.executable = ('raxmlHPC-PTHREADS-SSE3' if args.threads else
                'raxmlHPC-SSE3')

    with args.alignment_file as fp:
        sequences = SeqIO.parse(fp, args.input_format)
        raxml(sequences, args.output_tree, executable=args.executable,
              stats_path=args.stats, quiet=args.quiet,
              threads=args.threads, model=args.model, log_path=args.log,
              rapid_bootstrap=args.rapid_bootstrap,
              bootstrap_seed=args.bootstrap_seed,
              outgroup=args.outgroup,
              tmp_prefix=stripext(fp.name, True))

if __name__ == '__main__':
    main()
