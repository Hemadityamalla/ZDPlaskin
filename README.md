# ZDPlaskin
- You can do 0-D simulations using this. 
- The code is partially open source: http://www.zdplaskin.laplace.univ-tlse.fr/
- I try to maintain the code here for personal use. 
- There is a google group where people ask questions sometimes: https://groups.google.com/forum/#!forum/zdplaskin
- You can write stuff to your own custom output files for visualization, but they use something called QtPkaskin.
- QtPlaskin was written by Alejandro Luque: https://github.com/aluque/qtplaskin, but he doesnt really maintain it.
- These guys from france ported QtPlaskin to python 3 and they maintain their project: https://github.com/erwanp/qtplaskin

### Installing QtPlaskin on Fedora 30 and above:
These instructions were written for future reference and for the people in my group at CWI. They can be easily modified for installation in Ubuntu.
Similar installation instructions are given in the INSTALL.txt file in this repo: https://github.com/erwanp/qtplaskin
#### Creating a separate conda environment (yes, you need have anaconda installed, just install it to make life easier!)
- `conda create -n ZDPlaskin`
- `conda install -c conda-forge numpy scipy h5py matplotlib`
- `pip install pyqt5`
- `git clone https://github.com/erwanp/qtplaskin.git`
- `cd qtplaskin`
- `pip install -e .`
- To test if the installation was a success, type `qtplaskin` in the terminal and you will have a GUI opening. 

#### Running a simple ZDPlaskin simulation and visualizing the results using qtplaskin:
