function plot(c::MCMCChains, ptype::Symbol=:summary; args...)
  ptype == :summary ? [traceplot(c; args...) densityplot(c; args...)].' :
  ptype == :trace   ? traceplot(c; args...) :
  ptype == :density ? densityplot(c; args...) :
  ptype == :autocor ? autocorplot(c; args...) :
  ptype == :mean    ? meanplot(c; args...) :
    error("unsupported plot type $ptype")
end

function traceplot(c::MCMCChains; legend::Bool=false)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    plots[i] = plot(y=[[c.value[:,i,j] for j in 1:nchains]...],
                    x=[[c.range for j in 1:nchains]...],
                    Geom.line,
                    color=repeat([1:nchains], inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration"), Guide.ylabel("Value"),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function densityplot(c::MCMCChains; legend::Bool=false, 
                                    trim::(Real,Real)=(0.01,0.99))
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    qs = [quantile(c.value[:,i,j],[trim[1],trim[2]]) for j in 1:nchains]
    val = [c.value[ qs[j][1] .<= c.value[:,i,j] .<= qs[j][2],i,j] 
            for j in 1:nchains]
    plots[i] = plot(x=[val...], Geom.density,
                    color=repeat([1:nchains], inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Value"), Guide.ylabel("Density"),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function autocorplot(c::MCMCChains; maxlag::Integer=100, legend::Bool=false)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  lags = [(0:maxlag) * step(c.range)]
  ac = autocor(c, lags=[0:maxlag])
  for i in 1:nvars
    plots[i] = plot(y=[[ac.value[i,:,j]' for j in 1:nchains]...],
                    x=[[lags for j in 1:nchains]...],
                    Geom.line,
                    color=repeat([1:nchains], inner=[maxlag]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Lag"), Guide.ylabel("Autocorrelation"),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function meanplot(c::MCMCChains; legend::Bool=false)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    plots[i] = plot(y=[[cumsum(c.value[:,i,j])./[1:nrows] for j in 1:nchains]...],
                    x=[[c.range for j in 1:nchains]...],
                    Geom.line,
                    color=repeat([1:nchains], inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration"), Guide.ylabel("Mean"),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function draw(p::Array{Plot}; fmt::Symbol=:svg,
              filename::String="chainplot."string(fmt),
              width::MeasureOrNumber=8inch, height::MeasureOrNumber=8inch,
              nrow::Integer=3, ncol::Integer=2, byrow::Bool=true)

  in(fmt, [:png, :svg, :pdf, :ps]) ||
    error("$(fmt) is not a supported plot format")

  pp = nrow*ncol       ## plots per page
  ps = length(p)       ## number of plots
  np = iceil(ps/pp)    ## number of pages
  ex = pp - (ps % pp)  ## number of blank plots

  mat = Array(Canvas, pp)
  for page in 1:np
    nrem = ps - (page-1)*pp

    for j in 1:pp
      if j <= nrem
        mat[j] = render(p[(page-1)*pp+j])
      else
        mat[j] = canvas() ## pad with blank plots
      end
    end
    if page > 1
      println("Press ENTER to draw next plot")
      readline(STDIN)
    end
    if byrow
      result = reshape(mat,ncol,nrow).'
    else
      result = reshape(mat,nrow,ncol)
    end
    if fmt == :png
      draw(PNG(filename,width,height),gridstack(result))
    elseif fmt == :svg
      draw(SVG(filename,width,height),gridstack(result))
    elseif fmt == :pdf
      draw(PDF(filename,width,height),gridstack(result))
    elseif fmt == :ps
      draw(PS(filename,width,height),gridstack(result))
    end
  end

end
