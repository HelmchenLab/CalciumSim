CalciumSim
==========

Simulation framework for calcium imaging data

Henry Lütcke
Brain Research Institute
University of Zurich
Switzerland


To run the code add the folder modelCalcium and all sub-folders to your Matlab path.

On the command line type: S = modelCalcium

This will run the simulation and generate plots of simulated spikes and calcium traces. 

Parameters and data from the simulation are returned in the structure S.

Change any of the parameters in S and run the function again with the new parameter set: 
S = modelCalcium(S,1)

The second input argument is a plot flag, setting it to 0 suppresses all plots.

Here is an overview of the main parameters and their use:

A1 ... single AP calcium transient amplitude (in % dF/F)

A1sigma ... variability of A1 from AP to AP. This is modeled as a normal distribution with mean A1 and s.d. A1sigma*A1. Leave empty to use sam A1 for each AP.

tau1 ... decay time of calcium transient (in s)

tau1sigma ... variability of tau1 from AP to AP. This is modeled as a normal distribution with mean tau1 and s.d. tau1sigma*tau1. Leave empty to use same tau1 for each AP.

tauOn ... onset time of calcium transient (in s)

dur ... duration of the simulation (in s)

frameRate ... sampling rate of the simulated calcium signal.

snr ... Signal-to-noise ratio of the noisy calcium trace. SNR is defined as A1/SDnoise

spikeRate ... neuronal firing rate. To generate spike times, a Poisson process with this rate is used.

spikeTimes ... it is also possible to specify spike times explicitely using a cell array. All spike times in s.

recycleSpikeTimes ... flag to use the provided spike times (1) or generate new spike times at specified rate.

offset ... this is the time (in s) that is added before the beginning of the simulation to reach a steady state. important for high firing rates only.

reconAlg ... the spike reconstruction algorithm. this can be 'peeling' or 'none'.

peelOptions ... structure with options for the peeling algorithm. Possible fields are:
peelOptions.schmitt = [1.75 -1 0.3]; % schmitt trigger settings (high and low thresholds in terms of s.d. and min. duration in s)
peelOptions.peel_p.optimizeSpikeTimes = 1; % perform spike time optimization after peeling? useful for data with high temporal resolution. requires the optimization toolbox in Matlab.
peelOptions.peel_p.fitonset = 0; % fit onset of calcium transient. this is another option for improving the timing and could be used instead of the optimization approach.

In addition to modeling the calcium transient as single-exponential process, it is also possible to model a double-exponential process with a fast and a slow decay (see for example Grewe et al., Nat Meth, 2010). To do this, specify A1 / tau1 for the fast component and A2 / tau2 for the slow component.

After running the modelCalcium routine, output structure S contains additional fields:
data ... structure with simulated data:
			dff ... simulated DFF at original temporal resolution
			noisyDFF ... with noise added
			noisyDFFlowResT ... with noise and at target sampling rate
recon ... structure with output of reconstruction algorithm
			spikePredict ... predicted spike times by peeling
			peel ... the residual trace
			

