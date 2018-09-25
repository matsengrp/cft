
import subprocess


# Recording software versions
# ---------------------------

# First import some libs we'll need versions for
import ete3, Bio, dendropy #, pandas
# tripl version is a little messy because sometimes we load from local checkout
def tripl_version():
    try:
        import tripl
        return tripl.__version__
    except:
        from tripl import tripl
        return tripl.__version__

# The contract here is that a string val mapped to here is a command string to get a software version. A
# function value is called to get a version. And as long as it's not a function value, then it's assumed it's
# assumed the program name is something that `which` can be called on. If function value, assumed to be a lib
# and which is not called. For now... this is a little arbitrary and specific to our use case here.
software_versions = {
    'dnaml': None,
    'muscle': 'muscle -version',
    'seqmagick': 'seqmagick --version',
    'FastTree': None,
    #'prank': 'prank -v',
    'tripl': tripl_version,
    #'nestly': lambda: nestly.__version__,
    'ete3': lambda: ete3.__version__,
    'biopython': lambda: Bio.__version__,
    'scons': 'scons -v',
    'dendropy': lambda: dendropy.__version__,
    # For minadcl
    'rppr': 'rppr --version'
    }


def call_check_command(version_command):
    try:
        subprocess.check_output(version_command.split())
    except OSError as e:
        raise Exception('version command \'%s\' failed with \'%s\'' % (version_command, e))

def software_info(prog):
    version_command = software_versions[prog]
    return {'cft.software:name': prog,
            'cft.software:version': version_command() if callable(version_command) else (
                 call_check_command(version_command) if version_command else None),
            'cft.software:which': subprocess.check_output(['which', prog]) if not callable(version_command) else None}



def add_software_versions(w):
    @w.add_target('cft.build:software')
    def software(outdir, c):
        return [software_info(prog) for prog in software_versions]


