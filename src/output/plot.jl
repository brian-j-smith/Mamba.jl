#################### Posterior Plots ####################

#################### Generic Methods ####################

function draw(p::Array{Plot}; fmt::Symbol=:svg, filename::AbstractString="",
              width::MeasureOrNumber=8inch, height::MeasureOrNumber=8inch,
              nrow::Integer=3, ncol::Integer=2, byrow::Bool=true,
              ask::Bool=true)

  in(fmt, [:pdf, :pgf, :png, :ps, :svg]) ||
    error("$(fmt) is not a supported draw format")

  f(args...) = fmt == :pdf ? PDF(args...) :
               fmt == :pgf ? PGF(args...) :
               fmt == :png ? PNG(args...) :
               fmt == :ps  ? PS(args...)  : SVG(args...)

  isexternalfile = length(filename) > 0
  addextension = isexternalfile && search(filename, '.') == 0
  args = (width, height)

  pp = nrow * ncol               ## plots per page
  ps = length(p)                 ## number of plots
  np = ceil(Int, ps / pp)        ## number of pages

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
        mat[j] = context()
      end
    end
    result = byrow ? reshape(mat, ncol, nrow)' : reshape(mat, nrow, ncol)

    draw(img, gridstack(result))
  end

end

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
  ptype == :autocor      ? autocorplot(c; legend=legend, args...) :
  ptype == :bar          ? barplot(c; legend=legend, args...) :
  ptype == :density      ? densityplot(c; legend=legend, args...) :
  ptype == :mean         ? meanplot(c; legend=legend, args...) :
  ptype == :mixeddensity ? mixeddensityplot(c; legend=legend, args...) :
  ptype == :trace        ? traceplot(c; legend=legend, args...) :
                           error("unsupported plot type $ptype")
end


#################### Plot Engines ####################

function autocorplot(c::AbstractChains;
                     maxlag::Integer=round(Int, 10*log10(length(c.range))),
                     legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  lags = 0:maxlag
  ac = autocor(c, lags=collect(lags))
  for i in 1:nvars
    plots[i] = plot(y=vec(ac.value[i,:,:]),
                    x=repeat(collect(lags * step(c.range)), outer=[nchains]),
                    Geom.line,
                    color=repeat(c.chains, inner=[length(lags)]),
                    Scale.color_discrete(), Guide.colorkey("Chain"),
                    Guide.xlabel("Lag", orientation=:horizontal),
                    Guide.ylabel("Autocorrelation", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function barplot(c::AbstractChains; legend::Bool=false,
                 position::Symbol=:stack, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    S = unique(c.value[:,i,:])
    n = length(S)
    x = repmat(S, 1, nchains)
    y = zeros(n, nchains)
    for j in 1:nchains
      m = countmap(c.value[:,i,j])
      for k in 1:n
        if in(S[k], keys(m))
          y[k,j] = m[S[k]] / nrows
        end
      end
    end
    ymax = maximum(position == :stack ? mapslices(sum, y, [2]) : y)
    plots[i] = plot(x=vec(x), y=vec(y), Geom.bar(position=position),
                    color=repeat(c.chains, inner=[n]),
                    Scale.color_discrete(), Guide.colorkey("Chain"),
                    Guide.xlabel("Value", orientation=:horizontal),
                    Guide.ylabel("Density", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos),
                    Scale.y_continuous(minvalue=0.0, maxvalue=ymax))
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

function meanplot(c::AbstractChains; legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  val = mapslices(cummean, c.value, [1])
  for i in 1:nvars
    plots[i] = plot(y=vec(val[:,i,:]),
                    x=repeat(collect(c.range), outer=[nchains]),
                    Geom.line,
                    color=repeat(c.chains, inner=[length(c.range)]),
                    Scale.color_discrete(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration", orientation=:horizontal),
                    Guide.ylabel("Mean", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function mixeddensityplot(c::AbstractChains;
                          barbounds::Tuple{Real,Real}=(0,Inf), args...)
  plots = Array(Plot, size(c, 2))
  discrete = indiscretesupport(c, barbounds)
  plots[discrete] = plot(c[:,discrete,:], :bar; args...)
  plots[!discrete] = plot(c[:,!discrete,:], :density; args...)
  return plots
end

function traceplot(c::AbstractChains; legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array(Plot, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    plots[i] = plot(y=vec(c.value[:,i,:]),
                    x=repeat(collect(c.range), outer=[nchains]),
                    Geom.line,
                    color=repeat(c.chains, inner=[length(c.range)]),
                    Scale.color_discrete(), Guide.colorkey("Chain"),
                    Guide.xlabel("Iteration", orientation=:horizontal),
                    Guide.ylabel("Value", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function contourplot(c::ModelChains; bins::Int=100)
  nrows, nvars, nchains = size(c)
  combs = collect(combinations(1:nvars, 2))
  p = Array(Plot, length(combs))

  for i in 1:length(combs)
    D = vcat([c.value[:,combs[i],j] for j in 1:nchains]...)
    qx = collect(linspace(minimum(D[:,1]),maximum(D[:,1]),bins+1))
    qy = collect(linspace(minimum(D[:,2]),maximum(D[:,2]),bins+1))

    x = zeros(Float64, size(D,1))
    y = zeros(Float64, size(D,1))
    tot = zeros(Int, bins, bins)
    mx = map(j-> mean([qx[j], qx[j+1]]), 1:bins)
    my = map(j-> mean([qy[j], qy[j+1]]), 1:bins)
    for j in 1:size(D,1)
      idx1 = find(x-> qx[x] <= D[j,1] <= qx[x+1], collect(1:bins))[1]
      idx2 = find(x-> qy[x] <= D[j,2] <= qy[x+1], collect(1:bins))[1]
      tot[idx1, idx2] += 1
      x[j] = mx[idx1]
      y[j] = my[idx2]
    end
    p[i] = plot(x=mx, y=my, z=tot,
                Geom.contour, 
                Guide.xlabel(c.names[combs[i][1]], orientation=:horizontal),
                Guide.ylabel(c.names[combs[i][2]], orientation=:vertical),
                Guide.title("Contour of $(c.names[combs[i][1]]) vs $(c.names[combs[i][2]])")
               )
  end
  return p
end

#################### Auxiliary Functions ####################

function cummean{T<:Real}(x::Array{T})
  y = similar(x, Float64)
  xs = 0.0
  for i in 1:length(x)
    xs += x[i]
    y[i] = xs / i
  end
  y
end
