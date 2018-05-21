import fileinput
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('filename', type=str, help='enter tex file to add python highlighting')
args = parser.parse_args()
filename = args.filename

# add on style=pythonstyle to all lstlisting entries for code highlighting
with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
    for line in file:
        # remove pre-existing highlighting
        line = line.replace('[language=Python]', '')
        if not 'style' in line:
            # add pythonstyle highlighting if not already done
            line = line.replace('begin{lstlisting}', 'begin{lstlisting}[style=pythonstyle]')
        if not 'battery.png' in line and filename != 'laplace-la.tex':
            # rescale images to fit pagewidth
            line = line.replace('includegraphics{', 'includegraphics[width=\\textwidth]{')
        elif filename == 'laplace-la.tex':
            line = line.replace('includegraphics{', 'includegraphics[width=0.5\\textwidth]{')
        print(line, end='')

print('added lstlisting tags')
