function chain_plot(c::MCMCChains, ptype::Symbol=:trace; args...)
  ptype == :trace ? plot_trace(c) :
  ptype == :density ? plot_density(c) :
  ptype == :autocorrelation ? plot_autocorrelation(c) :
    error("unsupported plot type $type")
end

function plot_trace(c::MCMCChains)
	plots=Plot[]
	for j in 1:size(c)[2]
		push!(plots,plot([layer(y=c.value[:,j,i],x=[c.range], Geom.line, Theme(default_color=distinguishable_colors(size(c)[3])[i])) for i in 1:size(c)[3]]...))
	end
	return(plots)	
end

function plot_density(c::MCMCChains)

end

function plot_autocorrelation(c::MCMCChains)

end
