function plot(c::MCMCChains, ptype::Symbol=:summary; args...)
  ptype == :summary ? [traceplot(c) densityplot(c; args...)].' :
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
                                    trim::(Real,Real)=(0.025,0.975))
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

function autocorplot(c::MCMCChains;
                     maxlag::Integer=int(10*log10(length(c.range))),
                     legend::Bool=false)
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

function draw(p::Array{Plot}; fmt::Symbol=:svg, filename::String="",
              width::MeasureOrNumber=8inch, height::MeasureOrNumber=8inch,
              nrow::Integer=3, ncol::Integer=2, byrow::Bool=true)

  in(fmt, [:pdf, :png, :ps, :svg]) ||
    error("$(fmt) is not a supported draw format")

  if length(filename) == 0
    img = fmt == :pdf ? PDF(width, height) :
          fmt == :png ? PNG(width, height) :
          fmt == :ps  ? PS(width, height) :
                        SVG(width, height)
  else
    if search(filename, '.') == 0
      filename = string(filename, '.', fmt)
    end
    img = fmt == :pdf ? PDF(filename, width, height) :
          fmt == :png ? PNG(filename, width, height) :
          fmt == :ps  ? PS(filename, width, height) :
                        SVG(filename, width, height)
  end

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
    result = byrow ? reshape(mat, ncol, nrow).' : reshape(mat, nrow, ncol)
    draw(img, gridstack(result))
  end

end
