rm -fr dist build aroma_win aroma_win.zip
pyinstaller --windowed --onefile --name Aroma --hidden-import='PIL._tkinter_finder' aroma_gui.py
cd dist
mkdir aroma_win
mv Aroma aroma_win/
cp ../*png aroma_win/
cp ../user_aroma_constants.py aroma_win/
cp -r ../tests aroma_win/
zip -r aroma_win.zip aroma_win/
mv aroma_win.zip ../
cd ..
rm -fr dist build
