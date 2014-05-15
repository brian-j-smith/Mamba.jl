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
  ChainSummary(deepcopy(value), String[rownames...], String[colnames...],
               header)
end

function ChainSummary{T<:String,U<:String}(value::Matrix{Float64},
           rownames::Vector{T}, colnames::Vector{U}, header::String)
  dim = size(value)
  ChainSummary(reshape(value, dim[1], dim[2], 1), String[rownames...],
               String[colnames...], header)
end


#################### ChainSummary Base Methods ####################

function Base.show(io::IO, s::ChainSummary)
  value = map(x -> float16(x), s.value)
  if size(value, 3) == 1
    x = annotate(value[:,:,1], s.rownames, s.colnames)
  else
    x = mapslices(x -> annotate(x, s.rownames, s.colnames), value, [1,2])
  end
  Base.showarray(x, limit=false)
  print("\n")
end

function Base.showall(io::IO, s::ChainSummary)
  println(io, s.header)
  show(io, s)
end

function annotate(x::Matrix, rownames::Vector, colnames::Vector)
  hcat(["", rownames], vcat(colnames.', x))
end
