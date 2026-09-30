Code for simulating linking swarms using the Kirchhoff model for elastic fibers, using Matlab 2024a

In this repository:


-Code files

linking_swarm_simulator.m - code for simulating linking swarms. Output: simulation movie showing filaments coiling and linking.

Requires the following to be in the same directory: 
qualitative12.m
input_simulation_parameters.txt
input_initial_condition.txt
initial_conditions folder 
results folder 


IC_generator.m - code for generating initial conditions on which linking_swarm_simulator.m can be run. Output: picture of a valid initial condition of straight rods, initial condition is save to initial_conditions folder.

Requires the following to be in the same directory: 
qualitative12.m
initial_conditions folder


contactlink.m - code for plotting contact link per curve against time for a simulation of linked swarms, obtained from linking_swarm_simulator.m (above).

Requires the following to be in the same directory: 
results folder containing a file 'output_sim.mat' which will be analyzed


qualitative12.m - contains colormap used for plotting simulated swarms



-Text files

input_simulation_parameters.txt - plaintext file containing simulation input parameters and instructions on how to modify them. This input file is read by linking_swarm_simulator.m 

input_initial_condition.txt - plaintext file specifying which initial condition will be simulated, and instructions on how to choose the initial condition.



-Folders

initial_conditions - folder containing initial conditions which can be used as inputs for linking_swarm_simulator.

results_sample - sample results from 2 simulations (fast folding and slow folding from the initial condition, sticks20_1.mat given in the initial_conditions folder).

results - folder which will contain output files from running the simulation code.



-Data files

link_curvature_data.mat - file containing data for the 15 simulations runs given in results_sample. File consists of 4 arrays of size 20 x 201 x 15, which corresponds to 20 filaments x 201 timesteps x 15 runs. The arrays are labelled cLink1, cLink2, kList1, kList2. The arrays cLink* store contact link for each filament at each time across all 15 runs. The arrays kList* store the curvature for each filament at each time across all 15 runs. The arrays (cLink1, kList1) correspond to filament swarms which coil fast, and (cLink2, kList2) correspond to filament swarms which coil slowly. Fast and slow coiling correspond to Nbend = 0, 25000 respectively (see input_simulation_parameters.txt)




