export DUNEPY_DISABLE_PLOTTING=1

echo "Testing notebooks" > run.out
cd notebooks
find . -name "*.py" -print -exec python {} \; &>> ../run.out
cd ..

echo "Testing demo" >> run.out
cd demo
find . -name "*.py" -print -exec python {} \; &>> ../run.out
cd ..

echo "Testing test" >> run.out
cd test
find . -name "*.py" -print -exec python {} \; &>> ../run.out
cd ..
