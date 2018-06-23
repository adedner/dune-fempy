import fileinput
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('filename', type=str, help='enter tex file to add python highlighting')
args = parser.parse_args()
filename = args.filename

count = 0
with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
    for line in file:
        # remove pre-existing code highlighting
        line = line.replace('[language=Python]', '')
        # add pythonstyle highlighting if not already done
        if not 'style' in line:
            line = line.replace('begin{lstlisting}', 'begin{lstlisting}[style=pythonstyle]')
        # rescale images to fit pagewidth
        if filename != 'laplace-la.tex':
            line = line.replace('includegraphics{', 'includegraphics[width=\\textwidth, height=7cm, keepaspectratio]{')
        elif filename == 'laplace-la.tex':
            line = line.replace('includegraphics{', 'includegraphics[width=0.5\\textwidth]{')
        # add captions to images
        if filename == 'battery.tex' and line.startswith('\includegraphics'):
            if count == 0:
                print(line, end='')
            elif count == 1:
                print(line, end='')
                print('\caption{The initial plot of $c$ and $\phi$}')
            elif count ==2:
                print(line, end='')
                print('\caption{The plot after the final timestep}')
            count += 1
        elif filename == 'mcf.tex' and line.startswith('\includegraphics'):
            if count == 0:
                print(line, end='')
                print('\caption{The plot of the surface at three different timesteps}')
            elif count == 1:
                print(line, end='')
                print('\caption{Comparison of the error over time for varying levels of refinement}')
            count += 1
        elif filename == 'crystal.tex' and line.startswith('\includegraphics'):
            if count == 0:
                print(line, end='')
                print('\caption{The initial adapted grid and phase field}')
            elif count == 1:
                print(line, end='')
                print('\caption{The grid, phase field and temperature after the final timestep}')
            count += 1
        elif filename == 'laplace-adaptive.tex' and line.startswith('\includegraphics'):
            if count == 0:
                print(line, end='')
                print('\caption{The first three plots of the solution}')
            elif count == 1:
                print(line, end='')
                print('\caption{The second three plots of the solution}')
            elif count == 2:
                print(line, end='')
                print('\caption{The final three plots of the solution}')
            elif count == 3:
                print(line, end='')
                print('\caption{Zooming in on the re-entrant corner}')
            elif count == 4:
                print(line, end='')
                print('\caption{Plot of the level function of the grid}')
            count += 1
        # remove empty captions
        elif 'caption{png}' in line:
            pass
        else:
            print(line, end='')

print('added lstlisting tags')
