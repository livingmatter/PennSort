# PennSort
Spike Sorting Software

# Installation

1. Add /SpikeSort and all its subfolders to the MATLAB path.

2. mex /Spike Fitting/reconstructSpike.c

3. Compile /Clustering/ascii2dump.c and /Clustering/optics.cc as standalone executables.

4. When using a new MEA, it will be necessary to create a "geometry" function. For examples,
   see Arrays/2F30EL/geom_2F30EL.m and Arrays/256/geom_256.m. This should take an electrode group
   number as input and return a matrix of x-y coordinates of all the electrodes in that group.
   It should also return a cell array, ChannelGroups, listing all the electrode numbers belonging to each
   group.

   This geometry function should then be fed into the constructor for the "array" variable in SpikeSort.m.

5. Write a loading function for your data. It should take in a raw data filename and return a Matlab array
   whose rows are the time series from each channel. It should perform any filtering you deem necessary. 
   For our data, see "new_LoadAndFilterPenn.m" This function then must be assigned to the variable, "loadfun"
   in SpikeSort.m.

6. Adjust parameters as needed. A list of parameters and where to find them follows:

Parameters preceded by * likely need to be changed.

Parameters preceded by ** probably only need to be changed if using 
a sampling rate different from ours (10 kHz). They all represent some time interval, 
measured in samples, so just multiply by the appropriate factor to maintain the same
time interval. Some may need to be further changed if the actual spikes have different
time scales than ours.

Parameters preceded by ? probably don't need to be changed,
but you should think about it.

## SpikeSort.m:

- **onewaveformlength            Length, in samples, of event clips
- *numberchanns                  Number of electrodes
- **ndecimate                    Factor by which to interpolate templates
- *group                         If array has multiple electrode groups, which one to process
- *filelen                       Length, in samples, of one raw data file.
- *array                         Structure encoding array geometry. See note above.
- *threshold                     Spike threshold
- *noisethresh                   Noise threshold
- ?whiten                        Whether or not to spatially whiten

## Spike Extraction/extractSpikes.m:

*filelen                       As above
*array                         As above
*ChannelGroups                 Cell array giving electrode numbers within each group

## Statistics/getNoiseAndWhiteningFilter.m:

*ChannelGroups                 As above

## Clustering/order_points.m:

optics_minpts                  OPTICS parameter which, roughly speaking, controls "fragmentation" of clusters (bigger means less fragmentation but you miss more small clusters)      
optics_eps                     OPTICS parameter which roughly sets an upper bound on distances between elemnts of the same cluster. Nothing lost by making it really big, except maybe some speed.
interpfactor                   Interpolation factor, only used for aligning spikes prior to clustering.

## Clustering/makeTemplates.m:

?spikenbd                      Number of channels around peak to include in the template.
**ccwidth                      Maximum time shift to apply to events in aligning them for template building. Measured in samples.

## Statistics/updateClusterStats.m:

alpha                          These are all hyperparameters for the conjugate priors used in 
beta                           inferring amplitude and firing rate statistics.
nu
lambda
a
b

## Spike Fitting/preFit.m:

relevantthresh                  If a spike's amplitude is weaker than a template amplitude by more than this factor, don't bother trying to fit that template to the spike.
**correlationrange              Maximum time shift, relative to spike peak, to apply to templates in fitting. Measured in samples.
?widths                         Factors by which to stretch template widths. Maybe to get off the ground just try widths = [1] (no width fitting).

## Spike Fitting/cmultifit.m:

adjmethod                       "fast" or "slow." Use "fast" unless you run out of memory.
**timewindow                    During fitting, spikes are temporally cropped to a window of this size around the peak.
*interestingVThresh             Fitting terminates if no voltage crosses this threshold.
?interestingRThresh             Fitting terminates if the probability ratio of the best fitting template doesn't exceed this threshold.
?NMaxFitSpike                   Maximum number of templates to fit to a given event.
?spikenbd                       Number of channels around peak to include in event being fit. See makeTemplates.m.

# Instructions

## Using SpikeSort.m

The main file which controls all the different spike sorting stages is SpikeSort.m. It is organized using MATLAB's cell mode, so the idea is to advance through the cells one at a time to perform the different steps in the process.

At top there is a list of parameters, some of which will need to be changed for each spike sorting session. Many of them are described in "Installation.txt." The remaining ones, which will change from run to run, are described here: 
	datadir is the directory containing the raw data files.
	outputdir is where to put the spike times.
	clusterdir is where to put the clustered data.
	infileprefix is the prefix of the raw data files. The code assumes that they are named "infileprefix1"..."infileprefixN", and that consecutively 		      numbered files correspond to contiguously recorded data.
	clustfiles is a list of which raw data files to use in clustering. Choose ~2 min of data, spread uniformly throughout the experiment.
	noisefile is the data file used to get statistics on the noise. By default it is the first clustfile, and this probably doesn't need to be changed.
	fitfiles is the list of raw data files from which to get spike times.

The clustering code outputs to clusterdir:
	One "Cluster*.mat" file for each cluster. It contains a matrix, "C", with all the events in the cluster. "mCluster*.mat" is analogous for the merged 	clusters, and "mmCluster*.mat" for the second round of merging.

	"TemplateMatrix.mat" contains all the templates.
	"Stats.mat" contains parameters of the priors for each template.

The fitting code outputs to outputdir:
	One "ST__*" file for each raw data file. It contains "SpikeTimes," a cell array each of whose entries is the list of spike times (in samples) for the 	corresponding cluster, "amplist," a cell array whose entries are lists of amplitude scale factors for each fit of the corresponding template, and 	"stats," which gives the priors used in fitting that raw data file.
	
	"SpikeTimes" is of course the main output of the code. The other things can be useful for diagnostic purposes.


## Using the clustering GUI (plots.m)

After running OPTICS, SpikeSort.m saves the ordered events to a file prefixed by "DataToCluster_", placed in the specified cluster directory.
Load this file with File -> Open from .mat

The different events run horizontally. The waveforms on all channels are laid out vertically, with the clips from different channels concatenated together. The color gives the potential. You are looking for horizontal bands, representing a group of events that all peak at the same place with similar amplitude.

You interact with the GUI through the buttons on the toolbar. You can zoom in, zoom out, and pan horizontally with the magnifying glass and hand buttons.

To mark cluster boundaries, use the button with the empty box icon. Click once on the main screen where you want to start the cluster and once again where you want it to end. The region you marked off will become shaded. You can mark as many clusters in this way as you'd like; press the space bar when you want to stop. Note that each shaded region corresponds to a distinct cluster: you cannot mark disconnected clusters. Sometimes there are disconnected regions that obviously belong to the same cluster. Don't worry about this: mark them as separate clusters and you'll be able to merge them in the subsequent step.

To delete a cluster, click the icon with a box and an X, then click the cluster you want to remove.

The green arrow icons play movies of the events in a cluster but they're currently broken for arrays other than ours.

When done, you may use File -> Save Cluster Boundaries to save your work as a .clu file. There is also a corresponding "Open cluster boundaries." Close the GUI to advance SpikeSort to the next step.


## Merging clusters

The code will display pairs of overlapping templates in a new window. If you think they match well enough that the corresponding clusters should be merged, type 'y' at the MATLAB command line. If you don't want them merged, just hit enter to see the next one. The merged clusters are stored in files called "mCluster*.mat" in the specified cluster directory.


## Using the cluster splitting GUI (clustersplit.m)

When it first opens you'll get a file dialog box. Navigate to the cluster directory and select any one of the "mCluster" files.

Use the drop down menu at the upper-left to select which cluster you want to view. The arrow buttons (">" and "<") can also be used to move through this list. "Update plots" will display the events in the cluster. Select whether you want to see all the events or just their median with the radio buttons.

You're looking for channels on which there seem to be distinct subclusters. You'll see plenty of scatter, especially on channels far from the peak in the median. This is not important. What matters is whether the peak channel, or its near neighbors, doesn't look coherent. 

If a channel looks suspicious, click on it. The channel number should appear in the text box at the upper right. Then click "Histogram PC" to see a histogram of principal component weights on that channel. You can look at different principal components by using the arrow buttons or the drop-down menu. If a PC looks multimodal, the cluster should probably be split.

Click "Cut" to draw a threshold on the PC histogram. The cursor should change to crosshairs (there is a bug where this doesn't always happen) and you can click where on the histogram you'd like to place the threshold (that is, if the histogram is bimodal you want to click on the valley between the two peaks). You'll see a red line on the histogram where you clicked.

Click "Show L" and the main display will update to only show those events that fall to the left of the threshold. Similarly, "Show R" will show only those events that fall to the right. Inspect both of these, and if they really do look like distinct clusters, click "Save as new cluster." This will remove the currently displayed events from the cluster and place them in a new "mCluster" file at the end of the list.

If you don't want to split the cluster you can always click "Update plots" to see all the events again.


