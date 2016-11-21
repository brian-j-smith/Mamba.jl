#################### Model Graph ####################

function ModelGraph(m::Model)
  allkeys = keys(m, :all)
  g = DiGraph(length(allkeys))
  lookup = Dict(allkeys[i] => i for i in 1:length(allkeys))
  for key in keys(m)
    node = m[key]
    if isa(node, AbstractDependent)
      for src in node.sources
        add_edge!(g, lookup[src], lookup[key])
      end
    end
  end
  ModelGraph(g, allkeys)
end


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

graph(m::Model) = ModelGraph(m)

function graph2dot(m::Model)
  dag = ModelGraph(m)
  io = IOBuffer()
  write(io, "digraph MambaModel {\n")
  deps = keys(m, :dependent)
  for v in vertices(dag.graph)
    attr = Tuple{AbstractString, AbstractString}[]
    vkey = dag.keys[v]
    if vkey in deps
      node = m[vkey]
      if isa(node, AbstractLogical)
        push!(attr, ("shape", "diamond"))
      elseif isa(node, AbstractStochastic)
        push!(attr, ("shape", "ellipse"))
      end
      if isempty(node.monitor)
        push!(attr, ("style", "filled"),
                    ("fillcolor", "gray85"))
      end
    else
      push!(attr, ("shape", "box"),
                  ("style", "filled"),
                  ("fillcolor", "gray85"))
    end
    write(io, "\t\"")
    write(io, vkey)
    write(io, "\" [")
    write(io, join(map(x -> "$(x[1])=\"$(x[2])\"", attr), ", "))
    write(io, "];\n")
    for t in out_neighbors(dag.graph, v)
      write(io, "\t\t\"")
      write(io, vkey)
      write(io, "\" -> \"")
      write(io, dag.keys[t])
      write(io, "\";\n")
     end
  end
  write(io, "}\n")
  String(io)
end


#################### Auxiliary Functions ####################

function any_stochastic(dag::ModelGraph, v::Int, m::Model)
  found = false
  for t in out_neighbors(dag.graph, v)
    tkey = dag.keys[t]
    if isa(m[tkey], AbstractStochastic) || any_stochastic(dag, t, m)
      found = true
      break
    end
  end
  found
end

function gettargets(dag::ModelGraph, v::Int, terminalkeys::Vector{Symbol})
  values = Symbol[]
  for t in out_neighbors(dag.graph, v)
    tkey = dag.keys[t]
    push!(values, tkey)
    if !(tkey in terminalkeys)
      values = union(values, gettargets(dag, t, terminalkeys))
    end
  end
  values
end

function tsort(m::Model)
  dag = ModelGraph(m)
  dag.keys[topological_sort_by_dfs(dag.graph)]
end
