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
	dcols=distinguishable_colors(nchains,LCHab(70,60,240))
	for i in 1:nvars
		pl = Array(Layer,nchains)
		for j in 1:nchains
		   pl[j] = layer(y=c.value[:,i,j],x=c.range, Geom.line, 
						 Theme(default_color=dcols[j]))
		end
		plots[i] = plot(pl...,Guide.xlabel("Iteration"),Guide.ylabel(c.names[i]))
	end
	return(plots)	
end

function densityplot(c::MCMCChains)
	nrows, nvars, nchains = size(c)
	plots=Array(Plot, nvars)
	dcols=distinguishable_colors(nchains,LCHab(70,60,240))
	for i in 1:nvars
		pl = Array(Layer,nchains)
		for j in 1:nchains
		   pl[j] = layer(x=c.value[:,i,j], Geom.density, 
						 Theme(default_color=dcols[j]))
		end
		plots[i] = plot(pl...,Guide.xlabel(c.names[i]),Guide.ylabel("Density"))
	end
	return(plots)
end

function autocorplot(c::MCMCChains; maxlag::Integer=100)
	nrows, nvars, nchains = size(c)
	plots=Array(Plot, nvars)
	dcols=distinguishable_colors(nchains,LCHab(70,60,240))
	for i in 1:nvars
		pl = Array(Layer,nchains)
		for j in 1:nchains
			ac = autocor(c,lags=[1:maxlag])
			pl[j] = layer(y=ac.value[i,:,j], x=1:maxlag, Geom.line, 
						 Theme(default_color=dcols[j]))
		end
		plots[i] = plot(pl..., Guide.xlabel("Lag"), Guide.ylabel("Autocorrelation"))
	end
	return(plots)
end

function meanplot(c::MCMCChains)
	nrows, nvars, nchains = size(c)
	plots=Array(Plot, nvars)
	dcols=distinguishable_colors(nchains,LCHab(70,60,240))
	for i in 1:nvars
		pl = Array(Layer,nchains)
		for j in 1:nchains
		   pl[j] = layer(y=cumsum(c.value[:,i,j])/nrows, x=c.range, Geom.line, 
						 Theme(default_color=dcols[j]))
		end
		plots[i] = plot(pl...,Guide.xlabel("Iteration"),Guide.ylabel(c.names[i]))
	end
	return(plots)
end