#################### Model Graph ####################

#################### Display ####################

function draw(m::Model; filename::AbstractString="")
  dot = graph2dot(m)
  if length(filename) == 0
    print(dot)
  else
    if search(filename, '.') == 0
      filename = string(filename, ".dot")
    end
    f = open(filename, "w")
    write(f, dot)
    close(f)
  end
end

function graph(m::Model)
  g = graph(KeyVertex{Symbol}[], Edge{KeyVertex{Symbol}}[])
  lookup = Dict{Symbol,Int}()
  for key in keys(m, :all)
    lookup[key] = length(lookup) + 1
    add_vertex!(g, KeyVertex(lookup[key], key))
  end
  V = vertices(g)
  for key in keys(m)
    node = m[key]
    if isa(node, AbstractDependent)
      for src in node.sources
        add_edge!(g, V[lookup[src]], V[lookup[key]])
      end
    end
  end
  g
end

function graph2dot(m::Model)
  g = graph(m)
  io = IOBuffer()
  write(io, "digraph MambaModel {\n")
  deps = keys(m, :dependent)
  for v in vertices(g)
    attr = Dict{AbstractString,AbstractString}()
    if v.key in deps
      node = m[v.key]
      if isa(node, AbstractLogical)
        attr["shape"] = "diamond"
      elseif isa(node, AbstractStochastic)
        attr["shape"] = "ellipse"
      end
      if isempty(node.monitor)
        attr["style"] = "filled"
        attr["fillcolor"] = "gray85"
      end
    else
      attr["shape"] = "box"
      attr["style"] = "filled"
      attr["fillcolor"] = "gray85"
    end
    write(io, "\t\"")
    write(io, v.key)
    write(io, "\" [")
    write(io, join(map(x -> "$(x[1])=\"$(x[2])\"", attr), ", "))
    write(io, "];\n")
    for e in out_edges(v, g)
      t = target(e, g)
      write(io, "\t\t\"")
      write(io, v.key)
      write(io, "\" -> \"")
      write(io, t.key)
      write(io, "\";\n")
     end
  end
  write(io, "}\n")
  bytestring(io)
end


#################### Auxiliary Functions ####################

function any_stochastic(v::KeyVertex{Symbol}, g::AbstractGraph, m::Model)
  found = false
  for v in out_neighbors(v, g)
    if isa(m[v.key], AbstractStochastic) || any_stochastic(v, g, m)
      found = true
      break
    end
  end
  found
end

function gettargets(v::KeyVertex{Symbol}, g::AbstractGraph, m::Model)
  values = Symbol[]
  for v in out_neighbors(v, g)
    push!(values, v.key)
    if !isa(m[v.key], AbstractStochastic)
      values = union(values, gettargets(v, g, m))
    end
  end
  values
end

function tsort{T}(g::AbstractGraph{KeyVertex{T}, Edge{KeyVertex{T}}})
  V = topological_sort_by_dfs(g)
  map(v -> v.key, V)
end

function tsort(m::Model)
  tsort(graph(m))
end
