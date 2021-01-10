rm -fr dist build aroma_macos aroma_macos.tar.gz
pyinstaller --windowed --onefile --name Aroma --hidden-import='PIL._tkinter_finder' --osx-bundle-identifier=aroma.technion.edu aroma_gui.py
cd dist
mkdir aroma_macos
mv Aroma aroma_macos/
cp ../*png aroma_macos/
cp ../user_aroma_constants.py aroma_macos/
cp -r ../tests aroma_macos/
tar zcvf aroma_macos.tar.gz aroma_macos/
mv aroma_macos.tar.gz ../
cd ..
rm -fr dist build
