The FluidJ Robot Simulator serves as a comprehensive tool for analyzing both
forward and inverse kinematics in microrobots equipped with fluid joints.

Folder structure:

- /temp: saves the temporary files to create the mesh video
- /src: Code
- /data: stl file of the endplatform 
- /data/mesh: mesh of the fluid joints 
- /data/vids: output videos
- /data/figs: output figs


### How to run


Surface Evolver needs to be installed in the system.
Follow the official documentation: http://kenbrakke.com/evolver/html/install.htm

In terminal linux:

Clone this repository
'git clone https://github.com/Fran1702/FluidJ-Robot-Simulator'

Move to the folder and install the python requierements

'cd FluidJ-Robot-Simulator && python -m pip install -r requirements.txt'

Change to the src directory:

'cd src

Run the main_fk.py file as an example of the forward kinematics:

'python main_fkin.py '


