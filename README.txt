****************************** Copyright & License *****************************
CALIPER v1.1 is a software framework for the lifetime reliability evaluation of multicore architectures. Copyright (C) 2014 Politecnico di Milano.

This framework is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This framework is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details (http://www.gnu.org/licenses/).

Neither the name of Politecnico di Milano nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
********************************************************************************

** Introduction:
CALIPER v1.1 is a software framework for the lifetime reliability evaluation of N-out-of-M multicore architectures. A theoretical description of the implemented lifetime reliability model can be found in:
C. Bolchini, M. Carminati, M. Gribaudo, A. MIELE: A lightweight and open-source framework for the lifetime estimation of multicore systems. In Proc. of IEEE International Conference on Computer Design (ICCD), Seoul, South Korea, 2014, pp. 166-172. 

If you use CALIPER in your research, we would appreciate a citation to:

@inproceedings{ICCD2014,
 author = {C. Bolchini and M. Carminati and M. Gribaudo and A. Miele},
 title = {A lightweight and open-source framework for the lifetime estimation of multicore systems},
 booktitle = {Proc. IEEE Intl. Conf. Computer Design - ICCD},
 year = {2014},
 month = {Oct.},
 pages = {166--172}
}


** This folder contains:
_ this "README.txt" file
_ a "Makefile" for the project
_ a "src" folder where to find the source code
_ an example of input configuration file of a 2oo3 system ("example2oo3.txt")


** How to compile and run CALIPER (on Linux/Mac):
_ just type "make" from the root of this folder
_ run the caliper executable which has been placed in the "bin" folder ("./bin/caliper" from the root of this folder)
_ the following list of parameters is available (parameters in square brackets are to be considered as optional ones):
	-f INPUT_FILENAME: file containing the aging rate for any possible configuration
	[-m MIN_NUM_OF_CORES] : minimum number of cores for the QoS to be considered met; if not specified, it is derived from the aging rate file
	[-s RANDOM|num:num:num] : seed specification; three integer numbers separated by : need to be specified
	[-t NUMBER_OF_TESTS | -c CONFIDENCE_INTERVAL -r STOPPING_THRESHOLD] : Monte Carlo simulations to be performed: specify the exact number or the confidence interval and the stopping threshold
	[-p REMAPPING_PERIOD_IN_H] : duration (in hours) of the remapping period, in case multiple aging rates for each configuration are listed in the input file
	[-o R_OUTPUT_FILENAME] : output file where the system reliability function will be tabulated
	[-h] : print the help and specify the default values for the parameters


** Agin Rate File Generation:
When compiling the code, an executable to generate configuration files based on a dummy temperature model and a balanced workload distribution is created as well (/bin/dummythermal). By running such an executable, a valid input file for CALIPER is generated.
The command for running it is: 
	./bin/dummythermal <rows> <cols> <min_unit> <initialWorkloadXunit> <output_file>
where the parameters (all of them are mandatory) are:
	<rows> : is the number of rows in the architecture
	<cols> : is the number of cols in the architecture (for example, a 3 rows and 4 cols architecture consists of 12 cores)
	<min_unit> : the minimum number of working units to satisfy the QoS constraint
	<initialWorkloadXunit> : a number between 0 and 1 describing the initial workload of each single cores; the overall workload will be evenly redistributed among all the remaining cores after each failure
	<output_file> : the name of the file where the configuration file will be written


** Aging Rate File Format:
Each line in the aging rate input file must have the following format:
	CONFIGURATION AGING_RATES
where:
> CONFIGURATION: is a comma-separated string, listing the ids of failed cores in that configuration, sorted by time (with the oldest failure first). The initial configuration, where there is no failure, is indicated by means of the character ‘#’ (hash).
> AGING_RATES: is a blank space-separate string, listing the aging rate for the Reliability formula (based on the Weibull distribution) for each core in the architecture; the Beta parameter in the Weibull formula is specified with a macro in the source code. A number must be inserted also for failed core, even if it will be ignored by CALIPER.

The configuration file must contain all the possible failure combinations until a minimum number of core, which can be specified by the user or derived from the file itself. Moreover, in order to specify periodic remapping strategies, more than one entry for the same CONFIGURATION need to be specified on separate lines. example2oo3.txt is an example of configuration file of a 2oo3 system


** Sample Execution:
_ generate the aging rate configuration file for a 3x3 architecture working with at least 7 cores with a 0.5 initial workload by typing:
	./bin/dummythermal 3 3 7 0.5 3x3config.txt
_ then run CALIPER:
	./bin/caliper -f 3x3config.txt -s RANDOM -c 0.95 -r 0.1 -o 3x3rel.txt
Here, the above generated “3x3config.txt" configuration file is used; a RANDOM seed is selected; Monte Carlo simulations are run until a 95% confidence interval (-c 0.95) with a 0.1 stopping threshold (-r 0.1) is reached; the reliability curve is saved in the file "reliability.txt”, while other statistics are printed on the standard output.


** Updates:
1.0: 
 * Initial release
1.1:
 * Fixed various bugs (commit 4273)
 
 ** To signal bugs or for questions, please write to: antonio <dot> miele <at> polimi <dot> it
