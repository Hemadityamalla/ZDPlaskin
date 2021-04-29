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
- `conda install -c conda-forge numpy scipy h5py matplotlib`
- `pip install pyqt5`
- `git clone https://github.com/erwanp/qtplaskin.git`
- `cd qtplaskin`
- `pip install -e .` (pay attention to the . at the end.)
- To test if the installation was a success, type `qtplaskin` in the terminal and you will have a GUI opening. 

#### Running a simple ZDPlaskin simulation and visualizing the results using qtplaskin:
- The website of ZDPlaskin gives 3 examples. All of them are in this repository with folder names 'example1', 'example2', 'example3'. 
- You can read the manual of ZDPlaskin, type out the commands in there and move ahead. However, I tried to streamline this process a bit by using Makefiles.
- We want to simulate the evolution of chemical species under the influence of a certain electric field. To this end, we need four main files. 
  - First, the chemistry reaction mechanisms. This is specified in the `kinetics.inp` file. The auxilliary script `preprocessor` converts the `kinetics.inp` file into the `zdplaskin_m.f90` file.
  - Second, we need to obtain the rate constants for the reactions. A few of these, which depend on the E/N values are computed using BOLSIG-. For this we need to use either the `bolsig_x86_64.so` or the `bolsig.so`. 
  - Third, we need to specify the simulation parameters like gas temperature, initial charged specie densities, electric field, etc. This is done in the `main.f90` file. 
  - Finally, we need an ODE solver that solves our chemistry equations. This is done by using the `dvode_f90_m.f90` file. 
- For a simple simulation, you only need to edit the `kinetics.inp` and the `main.f90` files. Once these 4 files are setup/present you can run your simulation by typing `make run`. 
- Go to the 'Air_pulsed_simpleChemistry' folder for a simple example.  
