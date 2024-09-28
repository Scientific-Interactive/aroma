rm -fr dist build aroma-linux aroma-linux.tar.gz aroma-linux aroma-linux.tar.gz


# linux build
pyinstaller --windowed --onefile --name Aroma_GUI --hidden-import='PIL._tkinter_finder' --hidden-import=_tkinter aroma_gui.py
pyinstaller --onefile --name Aroma --hidden-import='PIL._tkinter_finder' aroma.py
pyinstaller --onefile --name Aroma_picmo aroma_picmo.py
cd dist
mkdir aroma-linux
mv Aroma aroma-linux/
mv Aroma_GUI aroma-linux/
mv Aroma_picmo aroma-linux/
cp ../*png aroma-linux/
cp ../user_aroma_constants.py aroma-linux/
cp -r ../tests aroma-linux/
tar zcvf aroma-linux.tar.gz aroma-linux/
mv aroma-linux.tar.gz ../
cd ..
rm -fr dist build
