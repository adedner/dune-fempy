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
        # add pythonstyle highlighting if not already done
        if not 'style' in line:
            line = line.replace('begin{lstlisting}', 'begin{lstlisting}[style=pythonstyle]')
        # rescale images to fit pagewidth
        if filename != 'laplace-la.tex':
            line = line.replace('includegraphics{', 'includegraphics[width=\\textwidth, height=7cm, keepaspectratio]{')
        elif filename == 'laplace-la.tex':
            line = line.replace('includegraphics{', 'includegraphics[width=0.5\\textwidth]{')
        # remove empty captions
        if 'caption{png}' in line:
            pass
        else:
            print (line, end='')

print('added lstlisting tags')
