
import SCons.Script as Script


# Set up command line arguments/options

Script.AddOption('--infiles',
        dest='infiles',
        metavar='FILE_LIST',
        default="test.yaml",
        help="""Specify ':' separated list of partis output directories to process on; if full path not specified,
        assumed to be in --base-datapath. Dataset names will be assined in relation to this --base-datapath
        if present. Note: symlinks may not work properly here unless they point to things in base-datapath also.""")

Script.AddOption('--asr-progs',
        dest='asr_progs',
        metavar='LIST',
        #default='dnaml:dnapars:raxml',
        default='dnaml',
        help="""Specify ':' separated list of ancestral state reconstruction programs to run. Options are `dnaml` and
        `dnapars`. Defaults to running dnaml.""")

Script.AddOption('--prune-strategies',
        dest='prune_strategies',
        metavar='LIST',
        default='min_adcl:seed_lineage',
        help="""Specify ':' separated list of pruning strategies. Options are 'min_adcl' and 'seed_lineage'.
        Defaults to both.""")

Script.AddOption('--test',
        dest='test_run',
        action='store_true',
        default=False,
        help="Setting this flag does a test run for just a couple of seeds")

Script.AddOption('--dataset-tag',
        dest='dataset_tag',
        metavar='TAG',
        help="Adds a tag to the beginning of the automatically generated dataset ids")

Script.AddOption('--outdir',
        dest='outdir',
        metavar='DIR',
        default='output',
        help="Directory in which to output results; defaults to `output`")


def get_options(env):
    # prefer realpath so that running latest vs explicit vN doesn't require rerun; also need for defaults below
    test_run = env.GetOption("test_run")
    return dict(
        infiles = env.GetOption('infiles').split(':'),
        asr_progs = env.GetOption('asr_progs').split(':'),
        prune_strategies = env.GetOption('prune_strategies').split(':'),
        dataset_tag = env.GetOption('dataset_tag') or ('test' if test_run else None),
        outdir_base = env.GetOption('outdir'))




