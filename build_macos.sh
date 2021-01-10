rm -fr dist build aroma-macos aroma-macos.tar.gz
pyinstaller --windowed --onefile --name Aroma --hidden-import='PIL._tkinter_finder' --hidden-import='tkinter' --osx-bundle-identifier=aroma.technion.edu --add-binary='/System/Library/Frameworks/Tcl.framework/Versions/Current/tclsh8.5':'tcl' aroma_gui.py
cd dist
mkdir aroma-macos
mv Aroma aroma-macos/
cp ../*png aroma-macos/
cp ../user_aroma_constants.py aroma-macos/
cp -r ../tests aroma-macos/
tar zcvf aroma-macos.tar.gz aroma-macos/
mv aroma-macos.tar.gz ../
cd ..
rm -fr dist build
