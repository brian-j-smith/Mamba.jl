function plot(c::AbstractChains, ptype::Vector{Symbol}=[:trace, :density];
              legend::Bool=false, args...)
  n = length(ptype)
  p = Array(Plot, n, size(c, 2))
  for i in 1:n
    showlegend = legend && i == n
    p[i,:] = plot(c, ptype[i]; legend=showlegend, args...)
  end
  p
end

function plot(c::AbstractChains, ptype::Symbol; legend::Bool=false, args...)
  ptype == :trace   ? traceplot(c; legend=legend, args...) :
  ptype == :density ? densityplot(c; legend=legend, args...) :
  ptype == :autocor ? autocorplot(c; legend=legend, args...) :
  ptype == :mean    ? meanplot(c; legend=legend, args...) :
  ptype == :summary ? error("use plot type [:trace, :density] instead of :summary") :
    error("unsupported plot type $ptype")
end

function traceplot(c::AbstractChains; legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    plots[i] = plot(y=vec(c.value[:,i,:]),
                    x=repeat([c.range;], outer=[nchains]),
                    Geom.line,
                    color=repeat(c.chains, inner=[length(c.range)]),
                    Scale.color_discrete(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration", orientation=:horizontal),
                    Guide.ylabel("Value", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function densityplot(c::AbstractChains; legend::Bool=false,
                     trim::Tuple{Real,Real}=(0.025,0.975), na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    val = Array(Vector{Float64}, nchains)
    for j in 1:nchains
      qs = quantile(c.value[:,i,j], [trim[1], trim[2]])
      val[j] = c.value[qs[1] .<= c.value[:,i,j] .<= qs[2], i, j]
    end
    plots[i] = plot(x=[val...;], Geom.density,
                    color=repeat(c.chains, inner=[length(c.range)]),
                    Scale.color_discrete(), Guide.colorkey("Chain"),
                    Guide.xlabel("Value", orientation=:horizontal),
                    Guide.ylabel("Density", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function autocorplot(c::AbstractChains;
                     maxlag::Integer=round(Int, 10*log10(length(c.range))),
                     legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  lags = 0:maxlag
  ac = autocor(c, lags=[lags;])
  for i in 1:nvars
    plots[i] = plot(y=vec(ac.value[i,:,:]),
                    x=repeat([lags * step(c.range);], outer=[nchains]),
                    Geom.line,
                    color=repeat(c.chains, inner=[length(lags)]),
                    Scale.color_discrete(), Guide.colorkey("Chain"),
                    Guide.xlabel("Lag", orientation=:horizontal),
                    Guide.ylabel("Autocorrelation", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function cummean{T<:Real}(x::Array{T})
  y = similar(x, Float64)
  xs = 0.0
  for i in 1:length(x)
    xs += x[i]
    y[i] = xs / i
  end
  y
end

function meanplot(c::AbstractChains; legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  val = mapslices(cummean, c.value, [1])
  for i in 1:nvars
    plots[i] = plot(y=vec(val[:,i,:]),
                    x=repeat([c.range;], outer=[nchains]),
                    Geom.line,
                    color=repeat(c.chains, inner=[length(c.range)]),
                    Scale.color_discrete(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration", orientation=:horizontal),
                    Guide.ylabel("Mean", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function draw(p::Array{Plot}; fmt::Symbol=:svg, filename::String="",
              width::MeasureOrNumber=8inch, height::MeasureOrNumber=8inch,
              nrow::Integer=3, ncol::Integer=2, byrow::Bool=true,
              ask::Bool=true)

  in(fmt, [:pdf, :png, :ps, :svg]) ||
    error("$(fmt) is not a supported draw format")

  f(args...) = fmt == :pdf ? PDF(args...) :
               fmt == :png ? PNG(args...) :
               fmt == :ps  ? PS(args...)  : SVG(args...)

  isexternalfile = length(filename) > 0
  addextension = isexternalfile && search(filename, '.') == 0
  args = (width, height)

  pp = nrow * ncol               ## plots per page
  ps = length(p)                 ## number of plots
  np = ceil(Int, ps / pp)    ## number of pages

  mat = Array(Context, pp)
  for page in 1:np
    if ask && page > 1 && !addextension
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
