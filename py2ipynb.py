#!/usr/bin/env python

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from nbformat import v3, v4
from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument('input',   type=str,  help='name of input python script file')
parser.add_argument('output',  type=str,  help='name of output jupyter notebook')
parser.add_argument('--image',   type=str,  help='name of image format',default="png")
parser.add_argument('--execute', type=bool, help='execute notebook?',default=True)
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
nb = v4.upgrade(nbook)  # Upgrade v3 to v4
if args.execute:
    ep = ExecutePreprocessor(timeout=-1,\
            extra_arguments=["--InlineBackend.figure_format="+args.image])
    ep.preprocess(nb,{})
with open(outputFile, 'w', encoding='utf-8') as f:
    nbformat.write(nb, f)

