rm -fr dist build aroma-linux aroma-linux.tar.gz aroma-linux-cmd aroma-linux-cmd.tar.gz

# gui build
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

# cmd build
pyinstaller --windowed --onefile --name Aroma --hidden-import='PIL._tkinter_finder' aroma.py
cd dist
mkdir aroma-linux-cmd
mv Aroma aroma-linux-cmd/
cp ../*png aroma-linux-cmd/
cp ../user_aroma_constants.py aroma-linux-cmd/
cp -r ../tests aroma-linux-cmd/
tar zcvf aroma-linux-cmd.tar.gz aroma-linux-cmd/
mv aroma-linux-cmd.tar.gz ../
cd ..
rm -fr dist build
