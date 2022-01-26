rm -fr dist build aroma-win aroma-win.zip
pyinstaller --windowed --onefile --name Aroma --hidden-import='PIL._tkinter_finder' --hidden-import='scipy.spatial.transform._rotation_groups' aroma_gui.py
cd dist
mkdir aroma-win
mv Aroma aroma-win/
cp ../*png aroma-win/
cp ../user_aroma_constants.py aroma-win/
7z a -r aroma-win.zip aroma-win/
mv aroma-win.zip ../
cd ..
rm -fr dist build
