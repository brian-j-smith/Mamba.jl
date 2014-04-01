#################### ChainSummary Type ####################

type ChainSummary
  data::Array{Float64,3}
  rownames::Vector{String}
  colnames::Vector{String}
  header::String

  function ChainSummary(data::Array{Float64,3}, rownames::Vector{String},
                        colnames::Vector{String}, header::String)
    dim = size(data)
    length(rownames) == dim[1] ||
      error("length of rownames not equal to number of rows")
    length(colnames) == dim[2] ||
      error("length of colnames not equal to number of columns")
    new(data, rownames, colnames, header)
  end
end


#################### ChainSummary Constructors ####################

function ChainSummary{T<:String,U<:String}(data::Array{Float64,3},
           rownames::Vector{T}, colnames::Vector{U}, header::String)
  ChainSummary(deepcopy(data), String[rownames...], String[colnames...], header)
end

function ChainSummary{T<:String,U<:String}(data::Matrix{Float64},
           rownames::Vector{T}, colnames::Vector{U}, header::String)
  dim = size(data)
  ChainSummary(reshape(data, dim[1], dim[2], 1), String[rownames...],
               String[colnames...], header)
end


#################### ChainSummary Base Methods ####################

function Base.show(io::IO, s::ChainSummary)
  if size(s.data)[3] == 1
    x = annotate(s.data[:,:,1], s.rownames, s.colnames)
  else
    x = mapslices(x -> annotate(x, s.rownames, s.colnames), s.data, [1,2])
  end
  showall(io, x)
  print("\n")
end

function Base.showall(io::IO, s::ChainSummary)
  println(io, s.header)
  show(io, s)
end


