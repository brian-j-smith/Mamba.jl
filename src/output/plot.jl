function plot(c::Chains, ptype::Vector{Symbol}; legend::Bool=false)
  n = length(ptype)
  p = Array(Plot, n, size(c, 2))
  for i in 1:n
    showlegend = legend && i == n
    p[i,:] = plot(c, ptype[i]; legend=showlegend)
  end
  p
end

function plot(c::Chains, ptype::Symbol=:summary; args...)
  ptype == :summary ? [traceplot(c) densityplot(c; args...)].' :
  ptype == :trace   ? traceplot(c; args...) :
  ptype == :density ? densityplot(c; args...) :
  ptype == :autocor ? autocorplot(c; args...) :
  ptype == :mean    ? meanplot(c; args...) :
    error("unsupported plot type $ptype")
end

function traceplot(c::Chains; legend::Bool=false)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    plots[i] = plot(y=[[c.value[:,i,j] for j in 1:nchains]...],
                    x=[[c.range for j in 1:nchains]...],
                    Geom.line,
                    color=repeat(c.chains, inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration", orientation=:horizontal),
                    Guide.ylabel("Value", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function densityplot(c::Chains; legend::Bool=false,
                     trim::(Real,Real)=(0.025,0.975))
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    qs = [quantile(c.value[:,i,j],[trim[1],trim[2]]) for j in 1:nchains]
    val = [c.value[ qs[j][1] .<= c.value[:,i,j] .<= qs[j][2],i,j] 
            for j in 1:nchains]
    plots[i] = plot(x=[val...], Geom.density,
                    color=repeat(c.chains, inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Value", orientation=:horizontal),
                    Guide.ylabel("Density", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function autocorplot(c::Chains; maxlag::Integer=int(10*log10(length(c.range))),
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
                    color=repeat(c.chains, inner=[maxlag]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Lag", orientation=:horizontal),
                    Guide.ylabel("Autocorrelation", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function meanplot(c::Chains; legend::Bool=false)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    plots[i] = plot(y=[[cumsum(c.value[:,i,j])./[1:nrows] for j in 1:nchains]...],
                    x=[[c.range for j in 1:nchains]...],
                    Geom.line,
                    color=repeat(c.chains, inner=[length(c.range)]),
                    Scale.discrete_color(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration", orientation=:horizontal),
                    Guide.ylabel("Mean", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function draw(p::Array{Plot}; fmt::Symbol=:svg, filename::String="",
              width::MeasureOrNumber=8inch, height::MeasureOrNumber=8inch,
              nrow::Integer=3, ncol::Integer=2, byrow::Bool=true)

  in(fmt, [:pdf, :png, :ps, :svg]) ||
    error("$(fmt) is not a supported draw format")

  f(args...) = fmt == :pdf ? PDF(args...) :
               fmt == :png ? PNG(args...) :
               fmt == :ps  ? PS(args...)  : SVG(args...)

  isexternalfile = length(filename) > 0
  addextension = isexternalfile && search(filename, '.') == 0
  args = (width, height)

  pp = nrow * ncol       ## plots per page
  ps = length(p)         ## number of plots
  np = iceil(ps / pp)    ## number of pages

  mat = Array(Context, pp)
  for page in 1:np
    if page > 1 && !addextension
      println("Press ENTER to draw next plot")
      readline(STDIN)
    end

    if isexternalfile
      fname = filename
      if addextension
        fname = string(fname, '-', page, '.', fmt)
      end
      args = (fname, width, height)
    end
    img = f(args...)

    nrem = ps - (page - 1) * pp
    for j in 1:pp
      if j <= nrem
        mat[j] = render(p[(page - 1) * pp + j])
      else
        mat[j] = context() ## pad with blank plots
      end
    end
    result = byrow ? reshape(mat, ncol, nrow).' : reshape(mat, nrow, ncol)

    draw(img, gridstack(result))
  end

end
