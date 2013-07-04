# Install METE module
git clone https://github.com/weecology/METE.git
cd METE
python setup.py install

#Install macroecotools module
git clone https://github.com/weecology/macroecotools.git
cd macroecotools
git checkout white-etal-2012
python setup.py install
