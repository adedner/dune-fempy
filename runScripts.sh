export DUNEPY_DISABLE_PLOTTING=1

make dune-fempy.py

find . -name "*.py" -print -exec python {} \; &>> run.out
