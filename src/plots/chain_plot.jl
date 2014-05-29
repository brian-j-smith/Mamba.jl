function chain_plot(c::MCMCChains, ptype::Symbol=:trace; args...)
  ptype == :trace ? plot_trace(c) :
  ptype == :density ? plot_density(c) :
  ptype == :autocorrelation ? plot_autocorrelation(c) :
    error("unsupported plot type $type")
end

function plot_trace(c::MCMCChains)
	plots=Plot[]
	dcols=distinguishable_colors(size(c)[3])
	for j in 1:size(c)[2]
		push!(plots,plot([layer(y=c.value[:,j,i],x=[c.range], Geom.line, Theme(default_color=dcols[i])) for i in 1:size(c)[3]]...))
	end
	return(plots)	
end

function plot_density(c::MCMCChains)
	plots=Plot[]
	dcols=distinguishable_colors(size(c)[3])
	for j in 1:size(c)[2]
		dens = [kde(c.value[:,j,i]) for i in 1:size(c)[3]]
		push!(plots,plot([layer(y=dens[i].density,x=dens[i].x, Geom.line, Theme(default_color=dcols[i])) for i in 1:size(c)[3]]...))
	end
	return(plots)
end

function plot_autocorrelation(c::MCMCChains)
	plots=Plot[]
	dcols=distinguishable_colors(size(c)[3])
	for j in 1:size(c)[2]
		ac = autocor(c,lags=[1:100])
		push!(plots,plot([layer(y=ac.value[j,:,i],x=[1:100], Geom.line, Theme(default_color=dcols[i])) for i in 1:size(c)[3]]...))
	end
	return(plots)
end
