#################### ChainSummary Type ####################

immutable ChainSummary
  value::Array{Float64,3}
  rownames::Vector{String}
  colnames::Vector{String}
  header::String

  function ChainSummary(value::Array{Float64,3}, rownames::Vector{String},
                        colnames::Vector{String}, header::String)
    dim = size(value)
    length(rownames) == dim[1] ||
      error("length of rownames not equal to number of rows")
    length(colnames) == dim[2] ||
      error("length of colnames not equal to number of columns")
    new(value, rownames, colnames, header)
  end
end


#################### ChainSummary Constructors ####################

function ChainSummary{T<:String,U<:String}(value::Array{Float64,3},
           rownames::Vector{T}, colnames::Vector{U}, header::String)
  ChainSummary(copy(value), String[rownames...], String[colnames...], header)
end

function ChainSummary{T<:String,U<:String}(value::Matrix{Float64},
           rownames::Vector{T}, colnames::Vector{U}, header::String)
  dim = size(value)
  ChainSummary(reshape(value, dim[1], dim[2], 1), String[rownames...],
               String[colnames...], header)
end


#################### ChainSummary Base Methods ####################

function Base.showall(io::IO, s::ChainSummary)
  println(io, s.header)
  show(io, s)
end

# write n ' ' characters to io
wrtsp(io::IO, n) = while (n -= 1) >= 0 write(io, ' ') end

function Base.show(io::IO, s::ChainSummary)
  rnwid = map(length,s.rownames)   # rowname widths
  mxrnwid = maximum(rnwid)
  cnwid = map(length,s.colnames)   # column name widths
  charv = mapslices(showoff, s.value, 1) # s.value as right alignable strings
  colwid = 1 + max(cnwid,vec(maximum(map(length,charv),[1,3])))
  m,n,f = size(charv)
  for k in 1:f
    ## write the column headers centered on the column widths
    wrtsp(io, mxrnwid)
    for j in 1:n
      nspace = colwid[j] - cnwid[j] - 1 # don't count the leading space
      nright = nspace >> 1              # divide by 2 rounding down
      wrtsp(io, 1 + nspace - nright)
      print(io, s.colnames[j])
      wrtsp(io, nright)
    end
    println(io)
    for i in 1:m
      wrtsp(io, mxrnwid - rnwid[i])
      print(io, s.rownames[i])
      for j in 1:n
        wrtsp(io, colwid[j] - length(charv[i,j,k]))
        print(io, charv[i,j,k])
      end
      println(io)
    end
    println(io)
  end
end


#################### PosteriorSummaries Type ####################

immutable PosteriorSummaries
  stats::ChainSummary
  quantiles::ChainSummary
end


#################### PosteriorSummaries Base Methods ####################

function Base.show(io::IO, ps::PosteriorSummaries)
  println(io, ps.stats.header)
  print(io, "Empirical Posterior Estimates:\n")
  show(io, ps.stats)
  print(io, "Quantiles:\n")
  show(io, ps.quantiles)
end
