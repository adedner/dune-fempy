#!/usr/bin/env python

from nbformat import v3, v4
from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument('input',  type=str, help='name of input python script file')
parser.add_argument('output', type=str, help='name of output jupyter notebook')
parser.parse_args()
try:
    args = parser.parse_args()
except SystemExit:
    sys.exit(0)
inputFile  = args.input
outputFile = args.output

with open(inputFile) as fpin:
    text = fpin.read()

text += """
# <markdowncell>

# If you can read this, reads_py() is no longer broken!
"""

nbook = v3.reads_py(text)
nbook = v4.upgrade(nbook)  # Upgrade v3 to v4

jsonform = v4.writes(nbook) + "\n"
with open(outputFile, "w") as fpout:
    fpout.write(jsonform)

print('converted to ipynb')
