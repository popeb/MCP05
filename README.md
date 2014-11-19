MCP05
=====

Scripts associated with Mouse ENCODE companion paper - Pope et al. Nature 2014 "Topologically-associating domains are stable units of replication-timing regulation"

###Replication domain boundary calling
DefineRDBs is designed to identify slope transitions along chromosomal profiles of DNA replication timing data, but can be applied in principle to any bivariate data.

DefineRDBs(x, Slope_E = 2.75e-6, Slope_L = 1e-06, Ends = 1e6, Data_Points = 30, Length_Max = 1e6, Length_Min = 2e5, RTU_Min = 0.55, Span = 35, Gap = 8e4, Gap_Dis = 125e3)

The input is a data table with at least three columns as described below:

x - table with columns "CHR" and "POSITION" indicating genomic coordinates of data points and additional columns containing replication timing data from individual samples

The following parameters can also be adjusted:

Slope_E - slope threshold for early Timing Transition Region (TTR) borders

Slope_L - slope threshold for late TTR borders

Ends - minimum distance between early TTR borders and chromosome ends

Data_Points - minimum number of data points within TTRs

Length_Max - maximum TTR size

Length_Min - minimum TTR size

RTU_Min - minimum RT difference across TTR

Span - span for loess smoothing

Gap - threshold to define gaps in data point spacing

Gap_Dis - minimum distance between early TTR borders and gaps in data point spacing

