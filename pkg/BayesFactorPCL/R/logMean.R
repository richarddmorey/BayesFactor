logMeanExpLogs = function(v)
{
	N = length(v)
	.Call("RLogMeanExpLogs", as.numeric(v), N, package="BayesFactor")
}

logCumMeanExpLogs = function(v)
{
	N = length(v)
	.Call("RLogCumMeanExpLogs", as.numeric(v), N, package="BayesFactor")
}

