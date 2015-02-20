#! /bin/usr/env python
# [Your name]
# [Date]
'''
Run your pipeline
'''

# PACKAGES
import sys
import os
import subprocess

# PARAMETERS
description = '''
----------------------------------------------------------------------
Your Fabulous Pipeline (C) 2015
----------------------------------------------------------------------
Brief description and terms of use.
----------------------------------------------------------------------
'''


# FUNCTIONS
def runstage(stage):
    '''Run python stage script'''
    if os.path.isdir(stage):
        sys.exit('[{0}] already run!'.format(stage))
    cmd = ['python', '-u', '{0}.py'.format(stage)]
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        while True:
            line = proc.stdout.readline().strip()
            if line == '':
                break
            else:
                print(line)
    except Exception as errmsg:
        sys.exit(errmsg)

if __name__ == '__main__':
    # TODO -- run stages here
