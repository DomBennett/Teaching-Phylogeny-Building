#! /bin/usr/env python
# Dom Bennett
# 19/02/2015
'''
Ensures dependencies are installed
'''

# PACKAGES
import re
import subprocess


# FUNCTIONS
def runCommand(args):
    """Run and return command"""
    # run, read and kill
    process = subprocess.Popen(args, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    info = process.stdout.read()
    process.kill()
    return info.strip()


def getVersion(args):
    """Return version number of program"""
    # check if program exists
    try:
        info = runCommand(args)
        # find version number, extract and strip of digits
        pattern = '(v|version)?\s?[0-9]\.[0-9]+'
        res = re.search(pattern, info)
        version = info[res.span()[0]:res.span()[1]]
        non_decimal = re.compile(r'[^\d.]+')
        version = non_decimal.sub('', version)
        return float(version)
    except OSError:
        return False

if __name__ == '__main__':
    print ('\nChecking packages and dependencies ....')
    # find missing dependencies
    missing_deps = []
    if getVersion(['mafft', '-u']) < 7.0:
        missing_deps.append('MAFFT v7+')
    raxml = getVersion(['raxml', '-version']) < 7.0
    raxmlHPC = getVersion(['raxmlHPC', '-version']) < 7.0
    if raxml and raxmlHPC:
        missing_deps.append('RAxML v7+')
    # find missing packages
    missing_packages = []
    try:
        import Bio
        del Bio
    except:
        missing_packages.append('Bio')
    try:
        import dendropy
        del dendropy
    except:
        missing_packages.append('dendropy')
    try:
        import taxon_names_resolver
        del taxon_names_resolver
    except:
        missing_packages.append('taxon-names-resolver')
    if missing_deps:
        print('Missing programs: {0}'.format(missing_deps))
    if missing_packages:
        print('Missing packages: {0}'.format(missing_packages))
    if not missing_packages and not missing_deps:
        print('Everything installed.')
    if not raxml:
        print('Your raxml cmd is: ["raxml"]')
    if not raxmlHPC:
        print('Your raxml cmd is: ["raxmlHPC"]')
    print('Done!\n')
