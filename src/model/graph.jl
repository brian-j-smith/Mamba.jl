#################### Model Graph Methods ####################

function any_stochastic(v::KeyVertex{Symbol}, g::AbstractGraph, m::Model)
  found = false
  for v in out_neighbors(v, g)
    if isa(m[v.key], MCMCStochastic) || any_stochastic(v, g, m)
      found = true
      break
    end
  end
  found
end

function draw(m::Model; filename::String="")
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

function gettargets(v::KeyVertex{Symbol}, g::AbstractGraph, m::Model)
  values = Symbol[]
  for v in out_neighbors(v, g)
    push!(values, v.key)
    if !isa(m[v.key], MCMCStochastic)
      values = union(values, gettargets(v, g, m))
    end
  end
  values
end

function graph(m::Model)
  g = graph(KeyVertex{Symbol}[], Edge{KeyVertex{Symbol}}[])
  lookup = (Symbol => Integer)[]
  for key in keys(m, :all)
    lookup[key] = length(lookup) + 1
    add_vertex!(g, KeyVertex(lookup[key], key))
  end
  V = vertices(g)
  for dep in keys(m, :dependent)
    for src in m[dep].sources
      add_edge!(g, V[lookup[src]], V[lookup[dep]])
    end
  end
  g
end

function graph2dot(m::Model)
  g = graph(m)
  str = "digraph Mamba.Model {\n"
  deps = keys(m, :dependent)
  for v in vertices(g)
    attr = (String => String)[]
    if in(v.key, deps)
      node = m[v.key]
      if isa(node, Logical)
        attr["shape"] = "diamond"
      elseif isa(node, MCMCStochastic)
        attr["shape"] = "ellipse"
      end
      if length(node.monitor) == 0
        attr["style"] = "filled"
        attr["fillcolor"] = "gray85"
      end
    else
      attr["shape"] = "box"
      attr["style"] = "filled"
      attr["fillcolor"] = "gray85"
    end
    str *= string(
      "\t\"", v.key, "\" [",
      join(map(x -> "$(x[1])=\"$(x[2])\"", attr), ", "),
      "];\n"
    )
    for e in out_edges(v, g)
      t = target(e, g)
      str *= "\t\t\"$(v.key)\" -> \"$(t.key)\";\n"
     end
  end
  str * "}\n"
end

function tsort{T}(g::AbstractGraph{KeyVertex{T}, Edge{KeyVertex{T}}})
  V = topological_sort_by_dfs(g)
  map(v -> v.key, V)
end

function tsort(m::Model)
  tsort(graph(m))
end
