export DUNEPY_DISABLE_PLOTTING=1

make dune-fempy.py

find . -maxdepth 1 -name "*.py" -print -exec python {} \; &>> run.out
