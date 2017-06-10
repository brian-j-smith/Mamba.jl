macro runtest(dir, prefix)
  quote
    fname = $(esc(dir)) * $(esc(prefix)) * ".jl"
    print("\n>>> Testing $fname\n\n")
    include(fname)
  end
end
