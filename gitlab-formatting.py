import fileinput
from argparse import ArgumentParser
import os, re, itertools

parser = ArgumentParser()
parser.add_argument('filename', type=str, help='enter tex file to add gitlab math')
args = parser.parse_args()
filename = args.filename

mathEnvs = ["gather","eqnarray","align","equation"]

count = 0
input = open(filename, 'r')
output = open('output', 'w')
lines = input.readlines()

for i in range(len(lines)):
    # remove pre-existing code highlighting
    lines[i] = lines[i].replace('[language=Python]', '')
    # add pythonstyle highlighting if not already done
    if any("\\begin{"+me in lines[i] for me in mathEnvs):
        output.write("```math\n")
    elif any("\\end{"+me in lines[i] for me in mathEnvs):
        output.write("```\n")
    else:
        counter = itertools.count(1)
        output.write( re.sub('\$',
              lambda m: "$`" if next(counter)%2==1 else "`$",
              lines[i])
            )
os.rename("output", filename)
