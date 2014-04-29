#################### MCMCModel Graph Methods ####################

function any_stochastic(v::ExVertex, g::AbstractGraph, m::MCMCModel)
  found = false
  for v in out_neighbors(v, g)
    if isa(m[v.label], MCMCStochastic) || any_stochastic(v, g, m)
      found = true
      break
    end
  end
  found
end

function getlinks(v::ExVertex, g::AbstractGraph, m::MCMCModel)
  values = String[]
  for v in out_neighbors(v, g)
    push!(values, v.label)
    if !isa(m[v.label], MCMCStochastic)
      values = union(values, getlinks(v, g, m))
    end
  end
  values
end

function graph(m::MCMCModel)
  g = graph(ExVertex[], Edge{ExVertex}[])
  lookup = (String=>Integer)[]
  for key in keys(m, :all)
    lookup[key] = length(lookup) + 1
    add_vertex!(g, ExVertex(lookup[key], key))
  end
  for key in keys(m, :input)
    v = vertices(g)[lookup[key]]
    v.attributes["shape"] = "box"
    v.attributes["style"] = "filled"
    v.attributes["fillcolor"] = "gray85"
  end
  for key in keys(m, :dep)
    v = vertices(g)[lookup[key]]
    node = m[key]
    if isa(node, MCMCLogical)
      v.attributes["shape"] = "diamond"
    elseif isa(node, MCMCStochastic)
      v.attributes["shape"] = "ellipse"
    end
    if !any(node.monitor)
      v.attributes["style"] = "filled"
      v.attributes["fillcolor"] = "gray85"
    end
    for key in node.sources
      add_edge!(g, vertices(g)[lookup[key]], v)
    end
  end
  g
end

function graph2dot(m::MCMCModel)
  g = graph(m)
  str = "digraph MCMCModel {\n"
  for v in vertices(g)
    str = str * string(
      "\t\"", v.label, "\" [",
      join(map(x -> "$(x[1])=\"$(x[2])\"", v.attributes), ", "),
      "];\n"
    )
    for e in out_edges(v, g)
      t = target(e, g)
      str = str * "\t\t\"$(v.label)\" -> \"$(t.label)\";\n"
     end
  end
  str * "}\n"
end

function graph2dot(m::MCMCModel, filename::String)
  f = open(filename, "w")
  write(f, graph2dot(m))
  close(f)
end

function plot(m::MCMCModel)
  stream, process = writesto(`dot -Tx11`)
  write(stream, graph2dot(m))
  close(stream)
end

function tsort(g::AbstractGraph{ExVertex, Edge{ExVertex}})
  V = topological_sort_by_dfs(g)
  map(v -> v.label, V)
end

function tsort(m::MCMCModel)
  tsort(graph(m))
end
