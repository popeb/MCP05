MCP05
=====
Scripts associated with Mouse ENCODE companion paper, Pope et al. Nature 2014 "Topologically-associating domains are stable units of replication-timing regulation" (doi:10.1038/nature13986)

###DefineRDBs.R
    DefineRDBs( x, Slope_E = 2.75e-6, Slope_L = 1e-06, Ends = 1e6, Data_Points = 30
                Length_Max = 1e6, Length_Min = 2e5, RTU_Min = 0.55, Span = 35,
                Gap = 8e4, Gap_Dis = 125e3 )

DefineRDBs is designed to identify slope transitions along chromosomal profiles of DNA replication timing data but can be applied in principle to any bivariate data.

#####Input
x  -  table with columns "CHR" and "POSITION" indicating genomic coordinates of data points and at least one additional column containing replication timing data per each individual sample

#####Arguments
Slope_E  -  slope threshold for early Timing Transition Region (TTR) borders

Slope_L  -  slope threshold for late TTR borders

Ends  -  minimum distance between early TTR borders and chromosome ends

Data_Points  -  minimum number of data points within TTRs

Length_Max  -  maximum TTR size

Length_Min  -  minimum TTR size

RTU_Min  -  minimum replication timing difference across TTR

Span  -  span for loess smoothing

Gap  -  threshold to define gaps in data point spacing

Gap_Dis  -  minimum distance between early TTR borders and gaps in data point spacing

### The deep autoencode analysis
See *dae_analysis* folder and follow

### RT_states.sh
The pipeline "RT_states.sh" will run chromHMM on ChIP-seq reads from mouse and human samples, learn a species-joint model, and produce whole genome segmentation in each sample. It extends the RT boundaries on both sides by the given bin size and bin number. Then it maps the segments of states onto the bins around the RT boundaries and produces a state matrix for each boundary.

Before running the pipeline, you should make the following preparation:
1. Install chromHMM (http://compbio.mit.edu/ChromHMM/) and put the directory in your path.
2. Prepare BED format sequencing tag files for all ChIP-input samples and pulldown samples, and separate input and pulldown into two directories. 
3. Create cell market tables according to the chromHMM manual.
4. Your input RT boundary file should have three columns and a header: "chr position strand".

Specify the parameters in the "Input set up" section of the pipeline, and run it by typing "sh RT_states.sh". The final outputs will have the same as "[SampleName]_states_onbinpoints".

###qnormRT.R
An R script to quantile-normalize Repli-chip and Repli-seq data sets.

#####Usage
qnormRT( RepliSeqBedGraphs , RepliChipFiles, cores="max" )

#####Arguments
RepliSeqBedGraphs - a character vector of length 1 or more containing file paths/names of Repli-seq data in bedGraph format.

RepliChipFiles - a character vector of length 1 or more containing file paths/names of Repli-chip data with the first column containing chromosome names, the second column containing chromosome coordinates, and additional columns containing replication timing scores, with each column corresponding to a sample.

cores - the number of threads to utilize. Default is "max" which sets "cores" to 1 less than the total available cores.

#####Details
qnormRT uses Repli-seq and Repli-chip data to generate a pool of replication timing scores which is subsequently used to sample scores for quantile normalization of the input replication timing profiles. The approach of pooling and subsampling scores allows for the quantile-normalization of data sets with different numbers of data points. The output of qnormRT is a new bedGraph of quantile-normalized data for each input sample. This function expects that both Repli-chip and Repli-seq profiles have similar, zero-centered distributions with similar IQRs.

#####Value
qnormRT creates a bedGraph file in the working directory for each sample corresponding to a quantile-normalized replication timing profile, and returns a character vector of file names for each bedGraph created.
