import fileinput
from argparse import ArgumentParser
import os, re, itertools

parser = ArgumentParser()
parser.add_argument('filename', type=str, help='enter tex file to add gitlab math')
args = parser.parse_args()
filename = args.filename

count = 0
input = open(filename, 'r')
output = open('output', 'w')
lines = input.readlines()
for i in range(len(lines)):
    # remove pre-existing code highlighting
    lines[i] = lines[i].replace('[language=Python]', '')
    # add pythonstyle highlighting if not already done
    if lines[i].startswith('\begin{equation*}'):
        output.write("```math")
    elif lines[i].startswith('\begin{equation}'):
        output.write("```math")
    elif lines[i].startswith('\begin{align*}'):
        output.write("```math")
    elif lines[i].startswith('\begin{align}'):
        output.write("```math")
    elif lines[i].startswith('\end{equation*}'):
        output.write("```")
    elif lines[i].startswith('\end{equation}'):
        output.write("```")
    elif lines[i].startswith('\end{align*}'):
        output.write("```")
    elif lines[i].startswith('\end{align}'):
        output.write("```")
    else:
        counter = itertools.count(1)
        output.write( re.sub('\$',
              lambda m: "$'" if next(counter)%2==1 else "'$",
              lines[i])
            )
os.rename("output", filename)
