rm -fr dist build aroma-win aroma-win.zip

# windows build
pyinstaller --windowed --onefile --name Aroma_GUI --hidden-import='PIL._tkinter_finder' --hidden-import='scipy.spatial.transform._rotation_groups' --hidden-import=_tkinter aroma_gui.py
pyinstaller --onefile --hidden-import='scipy.spatial.transform._rotation_groups' --name Aroma aroma.py
pyinstaller --onefile --hidden-import='scipy.spatial.transform._rotation_groups' --name Aroma_picmo aroma_picmo.py
cd dist
mkdir aroma-win
mv Aroma_GUI.exe aroma-win/
mv Aroma.exe aroma-win/
mv Aroma_picmo.exe aroma-win/
cp ../*png aroma-win/
cp ../user_aroma_constants.py aroma-win/
7z a -r aroma-win.zip aroma-win/
mv aroma-win.zip ../
cd ..
rm -fr dist build
