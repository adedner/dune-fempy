file=`basename $1 .py`
if [ ! -f $file.py ]; then
    echo "Script file $file.py not found!"
    echo "Provide the script to convert as first argument - the file extension has to by '.py'"
    exit 1
fi
echo converting $file.py to $file\_nb.ipynb and $file.rst

python py2ipynb.py $file.py $file_nb.ipynb --image="png"
jupyter nbconvert --to rst $file_nb.ipynb --output $file.rst
# sed -i '/^=======/a :download:`(notebook) <'$file'_nb.ipynb>` :download:`(script) <'$file'.py>`' $file.rst
# sed -i "/^=======/a :download:`(notebook) <$file_nb.ipynb>` :download:`(script) <$file.py>`" $file.rst
# sed -i "s/raw:: latex/math::/g" $file.rst
# sed -i "s/raw-latex/math/g" $file.rst
