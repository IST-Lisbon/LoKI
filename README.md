<p align="center">----------------------------- README FILE FOR THE LoKI-B SIMULATION TOOL -----------------------------<br>
<align="center">-----------------------------  10 steps to get acquainted with the tool  -----------------------------<br>
<align="center">(updated for version LoKI-B_v2.1.0)</p>

1. What's LoKI-B distribution license ?   
   The LisbOn KInetics Boltzmann (LoKI-B) is an open-source tool, licensed under the GNU general public license.  
   LoKI-B is freely available for users to perform electron kinetics calculations, and for expert researchers who are invited to continue testing the tool and/or to contribute for its development and improvement.

2. How to contact the developers ?   
   You are much welcome to report any problem or bug you may find using the code, or to send any feedback with suggestions for further developments.

   After downloading LoKI-B, and especially if you intend to interact with us, you are invited to send a short message   
   to: loki@tecnico.ulisboa.pt   
   with subject: <b>LoKI-B</b>   
   just giving your <b>name</b> and <b>affiliation</b>.   

3. What's LoKI-B ?   
   LoKI-B solves a time and space independent form of the two-term electron Boltzmann equation (EBE), for non-magnetised non-equilibrium low-temperature plasmas excited by DC/HF electric fields from different gases or gas mixtures. LoKI-B can also solve the two-term EBE for time-dependent (non-oscillatory) electric fields.
   LoKI-B includes electron-electron collisions, it handles rotational collisions adopting either a discrete formulation or a more convenient continuous approximation, and it accounts for variations in the number of electrons due to non-conservative events (ionisation and attachment) by assuming either a space-homogeneous exponential temporal growth or a time-constant exponential spatial growth of the electron density.

4. What's the programming language of LoKI-B ?   
   LoKI-B is developed with flexible and upgradable object-oriented programming under MATLAB, to benefit from its matrix-based architecture, adopting an ontology that privileges the separation between tool and data.

5. What are the input data of LoKI-B ?   
   On input, LoKI-B defines the operating work conditions (e.g. the applied reduced electric field and frequency, the gas pressure and temperature, the electron density, and the fractions of the different gases in the case of mixtures), the distribution of populations for the electronic, vibrational and rotational levels of the atomic / molecular gases considered, and the relevant sets of electron-scattering cross sections obtained from the open-access website LXCat (http://www.lxcat.net/). 

6. What are the output results of LoKI-B ?   
   On output, LoKI-B yields the isotropic and the anisotropic parts of the electron distribution function (the former usually termed the electron energy distribution function, EEDF), the electron swarm parameters, and the electron power absorbed from the electric field and transferred to the different collisional channels.   
   The latter parameters can be calculated using either the distribution function obtained from the solution to the EBE or some other form prescribed by the user, e.g. a generalized Maxwellian EEDF.

7. How to find your way in the code ?   
   After pulling the files in the repository, the LoKI-B folder contains   
   A) Subfolder "Documentation", with important documentation files - PLEASE READ THEM BEFORE USING THE CODE !!!   
   B) Subfolder "Code" containing   
   &ensp;(a) Several '\*.m' files corresponding to the MATLAB code, of which 'lokibcl.m' is the main file.  
   &ensp;(b) A subfolder "Input", containing the input files required for the simulations, organised as follows   
   &ensp;&ensp;i. Default configuration files 'default_lokib_setup.in' and 'default_lokib_pulse_setup.in'
   &ensp;&ensp;ii. A subfolder "Databases" with '\*.txt' files, containing different properties (masses, energies of levels, atomic/molecular constants, ...) for the gases used in the simulations.   
   &ensp;&ensp;iii. Several subfolders "Helium", "Nitrogen", ... with '\*.txt' files, containing the electron-scattering cross sections for the different gases used in the simulations, usually obtained from the open-access website LXCat (http://www.lxcat.net/).   
   &ensp;(c) A subfolder "PropertyFunctions", with several '\*.m' auxiliary functions for calculating the values of universal constants, some predefined distribution of species (Boltzmann, Treanor, ...), the energy of the levels according to some models, etc.    
   &ensp;(d) A subfolder "OtherAuxFunctions", with several '\*.m' auxiliary functions that provide the calculations of other properties, rate coeffcients, working conditions, etc.       
   &ensp;(e) A subfolder "Output" (eventually), where LoKI-B will write the output files resulting from the simulations.

8. How to run LoKI-B ?   
   The minimum requirements to run the code is a PC with an installation of MATLAB (R2017b or any later version).   
   We cannot ensure that all the features of LoKI will work properly under older versions.   

   LoKI-B runs upon calling the MATLAB function 'lokibcl(setupFile)'.   
   The end-user interacts with the code by specifying a particular "setup" for the simulation.    
   This setup is sent to the 'lokibcl()' function through the required input argument 'setupFile'.   
   The setup files should be located in [repository folder]/LoKI-B/Code/Input/ with a '.in' extension   
   (this is just a recommendation in order to keep the input folder organised; the setup files are just plain text files).   

   The distribution of LoKI-B includes a couple of default configuration files to help you make a first run of the code, following the sequence of steps below:   
   A) Open MATLAB   
   B) Navigate to the "Code" folder of your local copy of the repository:   
   &ensp;&ensp;>> cd [repository folder]/LoKI-B/Code/   
   C) Execute the following command in the MATLAB console:   
   &ensp;&ensp;>> lokibcl('default_lokib_setup.in')   
   D) The graphical user interface (GUI) should appear showing the solution(s) for the default setup file.

9. How to reference the code ?   
   LoKI-B is the result of the efforts of the Portuguese group N-Plasmas Reactive: Modeling and Engineering (N-PRiME), that decided to share the outcome of its research with the members of the Low-Temperature Plasmas community.   

   When using LoKI-B in your work, please give proper credits to the main developers, by adding the following citations:     
   [] Tejero A et al "The LisbOn KInetics Boltzmann solver" 2019 Plasma Sources Sci. Technol. 28 043001 (https://doi.org/10.1088/1361-6595/ab0537)       
   [] Tejero A et al "On the quasi-stationary approach to solve the electron Boltzmann equation in pulsed plasmas" 2021 Plasma Sources Sci. Technol. 30 065008 (https://doi.org/10.1088/1361-6595/abf858)        
   [available as open-access papers]

10. Acknowledgments   
   This work was partially funded by Portuguese FCT - Fundação para a Ciência e a Tecnologia, initially under project PTDC/FISPLA/1243/2014 (KIT-PLASMEBA), and currently under projects UIDB/50010/2020 and UIDP/50010/2020.

