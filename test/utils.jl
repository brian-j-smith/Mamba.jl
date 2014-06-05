macro runtest(dir, prefix)
  quote
    fname = $dir * $prefix * ".jl"
    print("\n>>> Testing $fname\n\n")
    include(fname)
  end
end
