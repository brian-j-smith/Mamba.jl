function plot(c::MCMCChains, ptype::Symbol=:trace; args...)
  ptype == :trace   ? traceplot(c) :
  ptype == :density ? densityplot(c) :
  ptype == :autocor ? autocorplot(c; args...) :
  ptype == :mean    ? meanplot(c) :
    error("unsupported plot type $ptype")
end

function traceplot(c::MCMCChains)
	nrows, nvars, nchains = size(c)
	plots=Array(Plot, nvars)
	for i in 1:nvars
		plots[i] = plot(y=[[c.value[:,i,j] for j in 1:nchains]...],
						x=[[c.range for j in 1:nchains]...], 
						Geom.line, 
						color=repeat([j for j in 1:nchains], inner=[length(c.range)]),
						Scale.discrete_color(),	Guide.colorkey("Chain"), 
						Guide.xlabel("Iteration"),Guide.ylabel(c.names[i]))
	end
	return(plots)	
end

function densityplot(c::MCMCChains)
	nrows, nvars, nchains = size(c)
	plots=Array(Plot, nvars)
	for i in 1:nvars
		plots[i] = plot(x=[[c.value[:,i,j] for j in 1:nchains]...], Geom.density, 
						color=repeat([j for j in 1:nchains], inner=[length(c.range)]),
						Scale.discrete_color(),	Guide.colorkey("Chain"), 
						Guide.xlabel(c.names[i]),Guide.ylabel("Density"))
	end
	return(plots)
end

function autocorplot(c::MCMCChains; maxlag::Integer=100)
	nrows, nvars, nchains = size(c)
	plots=Array(Plot, nvars)
	ac = autocor(c,lags=[1:maxlag])
	for i in 1:nvars
		plots[i] = plot(y=[[ac.value[i,:,j] for j in 1:nchains]...],
						x=[[1:maxlag for j in 1:nchains]...], 
						Geom.line, 
						color=repeat([j for j in 1:nchains], inner=[length(c.range)]),
						Scale.discrete_color(),	Guide.colorkey("Chain"), 
						Guide.xlabel("Lag"),Guide.ylabel("Autocorrelation"),
						Guide.title(c.names[i]))
	end
	return(plots)
end

function meanplot(c::MCMCChains)
	nrows, nvars, nchains = size(c)
	plots=Array(Plot, nvars)
	for i in 1:nvars
		plots[i] = plot(y=[[cumsum(c.value[:,i,j])/nrows for j in 1:nchains]...],
						x=[[c.range for j in 1:nchains]...], 
						Geom.line, 
						color=repeat([j for j in 1:nchains], inner=[length(c.range)]),
						Scale.discrete_color(),	Guide.colorkey("Chain"), 
						Guide.xlabel("Iteration"),Guide.ylabel(c.names[i]))
	end
	return(plots)
end