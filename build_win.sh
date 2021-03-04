rm -fr dist build aroma-win aroma-win.zip
pyinstaller --windowed --onefile --name Aroma --hidden-import='PIL._tkinter_finder' aroma_gui.py
cd dist
mkdir aroma-win
mv Aroma aroma-win/
cp ../*png aroma-win/
cp ../user_aroma_constants.py aroma-win/
cp -r ../tests aroma-win/
7z a -r aroma-win.zip aroma-win/
mv aroma-win.zip ../
cd ..
rm -fr dist build
