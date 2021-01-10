rm -fr dist build
pyinstaller --windowed --onefile --name Aroma --hidden-import='PIL._tkinter_finder' aroma_gui.py
cd dist
mkdir aroma_linux
mv Aroma aroma_linux/
cp ../*png aroma_linux/
cp ../user_aroma_constants.py aroma_linux/
cp -r ../tests aroma_linux/
tar zcvf aroma_linux.tar.gz aroma_linux/
mv aroma_linux.tar.gz ../
cd ..
rm -fr dist build
