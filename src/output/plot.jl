function plot(c::MCMCChains, ptype::Symbol=:trace; args...)
  ptype == :trace   ? traceplot(c) :
  ptype == :density ? densityplot(c) :
  ptype == :autocor ? autocorplot(c; args...) :
  ptype == :mean    ? meanplot(c) :
    error("unsupported plot type $ptype")
end

function traceplot(c::MCMCChains)
  nrows, nvars, nchains = size(c.value)
  plots=Array(Plot, nvars)
  for i in 1:nvars
    plots[i] = plot(y=[[c.value[:,i,j] for j in 1:nchains]...],
                    x=[[c.range for j in 1:nchains]...], 
                    Geom.line, 
                    color=repeat([1:nchains], inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"), 
                    Guide.xlabel("Iteration"),Guide.ylabel(c.names[i]))
  end
  return(plots)  
end

function densityplot(c::MCMCChains)
  nrows, nvars, nchains = size(c.value)
  plots=Array(Plot, nvars)
  for i in 1:nvars
    plots[i] = plot(x=[[c.value[:,i,j] for j in 1:nchains]...], Geom.density, 
                    color=repeat([1:nchains], inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"), 
                    Guide.xlabel(c.names[i]),Guide.ylabel("Density"))
  end
  return(plots)
end

function autocorplot(c::MCMCChains; maxlag::Integer=100)
  nrows, nvars, nchains = size(c.value)
  plots=Array(Plot, nvars)
  lags = [1:maxlag]*step(c.range)
  ac = autocor(c,lags=[1:maxlag])
  for i in 1:nvars
    plots[i] = plot(y=[[ac.value[i,:,j]' for j in 1:nchains]...],
                    x=[[lags for j in 1:nchains]...], 
                    Geom.line, 
                    color=repeat([1:nchains], inner=[maxlag]),
                    Scale.discrete_color(), Guide.colorkey("Chain"), 
                    Guide.xlabel("Lag"),Guide.ylabel("Autocorrelation"),
                    Guide.title(c.names[i]))
  end
  return(plots)
end

function meanplot(c::MCMCChains)
  nrows, nvars, nchains = size(c.value)
  plots=Array(Plot, nvars)
  for i in 1:nvars
    plots[i] = plot(y=[[cumsum(c.value[:,i,j])./[1:nrows] for j in 1:nchains]...],
                    x=[[c.range for j in 1:nchains]...], 
                    Geom.line, 
                    color=repeat([1:nchains], inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"), 
                    Guide.xlabel("Iteration"),Guide.ylabel(c.names[i]))
  end
  return(plots)
end

function draw{T<:Plot}(p::Vector{T}; fmt::Symbol=:svg, 
                                     filename::String="chainplot."string(fmt), 
                                     width::MeasureOrNumber=12cm, 
                                     height::MeasureOrNumber=8cm,
                                     nrow::Integer=2, ncol::Integer=2, 
                                     byrow::Bool=false)
  if !(fmt in [:png, :svg, :pdf, :ps])
    error("$(fmt) is not a supported plot format")
  end
  pp = nrow*ncol       ## plots per page
  ps = length(p)       ## number of plots
  np = div(ps,pp)+1    ## number of pages
  ex = pp - (ps % pp)  ## number of blank plots

  for page in 1:np
    nrem = ps - (page-1)*pp
    mat = Canvas[]

    for j in 1:min(nrem,pp)
      push!(mat, render(p[(page-1)*pp+j]))
    end
    while length(mat) < pp
      push!(mat, canvas())
    end
    if page > 1
      println("Press any key to see next plot")
      key = readline(STDIN)
    end
    if fmt == :png
      draw(PNG(filename,width,height),gridstack(reshape(mat,nrow,ncol)))
    elseif fmt == :svg
      draw(SVG(filename,width,height),gridstack(reshape(mat,nrow,ncol)))
    elseif fmt == :pdf
      draw(PDF(filename,width,height),gridstack(reshape(mat,nrow,ncol)))
    elseif fmt == :ps
      draw(PS(filename,width,height),gridstack(reshape(mat,nrow,ncol)))
    end
  end

end

