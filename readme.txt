For the convenience of implementation, some instructions are provided as follows.

For Narrowband scenario,
1. First please set parameters of the mmWave MIMO system in "channel_realizaiton.m" and run it.
2. Save all the matrices as a .mat file, which will be used in further simulation, e.g., channel matrix, array response vectors, and optimal precoders.
3. Choose a hybrid precoding algorithm (different algorithms are packed into different folders), run the "main_XXX.m". Note that at the beginning of each "main_XXX.m" file, a command loading the saved channel matrices is used. Users should revise this command accordingly corresponding to the "XXX.mat" file name of the saved channel matrices.

** NOTE: For the "SDR_AltMin" and "SIC" algorithms, the authors used a cvx (Version 2.1, Build 1103 (9714d49)). Users should include a cvx package when these algorithms are implemented.

For OFDM scenario,
Just choose an algorithm folder and run the "main_SNR.m".


** Users can also carefully use the "parfor" command to implement the Matlab parallel pool to accelerate the simulation.

** The solvers for manifold optimization in the codes are from Manopt, details can be found on http://www.manopt.org/, or
[Ref] N. Boumal, B. Mishra, P.-A. Absil, and R. Sepulchre, ¡°Manopt, a Matlab toolbox for optimization on manifolds,¡± J. Mach. Learn. Research, vol. 15, pp. 1455¨C1459, Jan. 2014.