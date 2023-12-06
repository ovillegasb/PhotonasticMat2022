#!/bin/bash

# Enviroment SimMOL
conda create --name SimMOL scipy jupyterlab numpy matplotlib seaborn pandas
conda activate SimMOL
conda config --add channels conda-forge
# conda install --channel conda-forge pymatgen
conda install pymatgen
conda install mdtraj
conda install mdanalysis
conda install networkx
conda install statsmodels

# python modules
python -m site --user-site
mkdir -p $(python -m site --user-site)

# create new .pth file with our path
echo "$HOME/foo/bar" > "$SITEDIR/somelib.pth"

# Extra
conda install -n SimMOL -c conda-forge nglview
conda install -c conda-forge nodejs
conda install -c conda-forge/label/gcc7 nodejs
conda install -c conda-forge/label/cf201901 nodejs
conda install -c conda-forge/label/cf202003 nodejs
jupyter labextension install nglview-js-widgets
jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-nglview-jsmol

jupyter labextension update --all
jupyter labextension list
jupyter labextension install @jupyterlab/celltags
jupyter serverextension enable --py jupyterlab --sys-prefix

# Others
sudo dnf install pandoc
conda install ipympl  # %matplotlib widget

conda install -n base -c conda-forge jupyterlab_widgets
conda install -n SimMOL -c conda-forge ipywidgets
conda install -c conda-forge mpi4py
conda install memory_profiler
conda install statsmodels
pip install pyppeteer

# conda install -c auto threadpool
# conda install -c auto multiprocessing


# Enviroment AmberTools
conda create --name AmberTools22
conda activate AmberTools22
conda install -c conda-forge ambertools=22 compilers
conda install -c conda-forge acpype
# acpype -i test/mol2/hpt.mol2 -a gaff2 -o gmx


# To convert Jupyter to slides
jupyter nbconvert --to slides --no-input --post serve presentation.slides.ipynb
# then to pdf using pandoc
pandoc presentation.slides.html -t latex -o slides.pdf

# OTHERS
pip install jupyter_contrib_nbextensions
jupyter contrib nbextension install --user
jupyter nbextension enable codefolding/main
python -m pip install -U notebook-as-pdf
pyppeteer-install
pip install 'PyPDF2<3.0'

jupyter nbconvert --to pdfviahtml --no-input --no-prompt ARTICLE_results.ipynb

jupyter nbconvert --to webpdf --no-input --no-prompt ARTICLE_results.ipynb

# wkhtmltopdf

sudo dnf install wkhtmltopdf

jupyter nbconvert --to html Chapter1.ipynb
jupyter nbconvert --to html Chapter2.ipynb

wkhtmltopdf --enable-internal-links -L 10mm -R 9.5mm -T 10mm -B 9.5mm Chapter1.html Chapter1.pdf
wkhtmltopdf --enable-internal-links -L 10mm -R 9.5mm -T 10mm -B 9.5mm Chapter2.html Chapter2.pdf
./cpdf Chapter1.pdf Chapter2.pdf -o Ebook.pdf 

# python modules
python -m site --user-site
mkdir -p $(python -m site --user-site)

# create new .pth file with our path
echo "$HOME/foo/bar" > "$SITEDIR/somelib.pth"

# Run
chmod +x ./GITPROYECTS/PhotonasticMat/stamptools/__main__.py
ln -s /home/ovillegas/GITPROYECTS/PhotonasticMat/stamptools/__main__.py /home/ovillegas/.local/bin/stamptools


# Widgets
# https://ipywidgets.readthedocs.io/en/7.x/user_install.html
pip install ipywidgets
#jupyter labextension install @jupyter-widgets/jupyterlab-manager
#jupyter nbextension enable --py widgetsnbextension --sys-prefix
jupyter nbextension enable --py widgetsnbextension --sys-prefix
jupyter labextension install jupyter-matplotlib


# JupyterLab
jupyter labextension --help
jupyter labextension list
jupyter lab build --minimize=False