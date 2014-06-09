function plot(c::MCMCChains, ptype::Symbol=:summary; args...)
  ptype == :summary ? hcat(traceplot(c), densityplot(c)).' :
  ptype == :trace   ? traceplot(c) :
  ptype == :density ? densityplot(c) :
  ptype == :autocor ? autocorplot(c; args...) :
  ptype == :mean    ? meanplot(c) :
    error("unsupported plot type $ptype")
end

function traceplot(c::MCMCChains)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  for i in 1:nvars
    plots[i] = plot(y=[[c.value[:,i,j] for j in 1:nchains]...],
                    x=[[c.range for j in 1:nchains]...],
                    Geom.line,
                    color=repeat([1:nchains], inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration"), Guide.ylabel(c.names[i]))
  end
  return plots
end

function densityplot(c::MCMCChains)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  for i in 1:nvars
    plots[i] = plot(x=[[c.value[:,i,j] for j in 1:nchains]...], Geom.density,
                    color=repeat([1:nchains], inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel(c.names[i]), Guide.ylabel("Density"))
  end
  return plots
end

function autocorplot(c::MCMCChains; maxlag::Integer=100)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  lags = [(0:maxlag) * step(c.range)]
  ac = autocor(c, lags=[0:maxlag])
  for i in 1:nvars
    plots[i] = plot(y=[[ac.value[i,:,j]' for j in 1:nchains]...],
                    x=[[lags for j in 1:nchains]...],
                    Geom.line,
                    color=repeat([1:nchains], inner=[maxlag]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Lag"), Guide.ylabel("Autocorrelation"),
                    Guide.title(c.names[i]))
  end
  return plots
end

function meanplot(c::MCMCChains)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  for i in 1:nvars
    plots[i] = plot(y=[[cumsum(c.value[:,i,j])./[1:nrows] for j in 1:nchains]...],
                    x=[[c.range for j in 1:nchains]...],
                    Geom.line,
                    color=repeat([1:nchains], inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration"), Guide.ylabel(c.names[i]))
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
