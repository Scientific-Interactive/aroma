rm -fr dist build aroma-linux aroma-linux.tar.gz
pyinstaller --windowed --onefile --name Aroma --hidden-import='PIL._tkinter_finder' aroma_gui.py
cd dist
mkdir aroma-linux
mv Aroma aroma-linux/
cp ../*png aroma-linux/
cp ../user_aroma_constants.py aroma-linux/
cp -r ../tests aroma-linux/
tar zcvf aroma-linux.tar.gz aroma-linux/
mv aroma-linux.tar.gz ../
cd ..
rm -fr dist build
