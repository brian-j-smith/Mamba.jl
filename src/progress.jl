#################### ChainProgress Type ####################

type ChainProgress
  chain::Integer
  iters::Integer
  counter::Integer

  ChainProgress(chain::Integer, iters::Integer) = new(chain, iters, 0)
end


#################### ChainProgress Methods ####################

function next!(p::ChainProgress)
  counter = p.counter
  p.counter += 1

  dt = 10.0 / p.iters
  if floor(p.counter * dt) > counter * dt || counter == 0
    print(STDOUT, p)
  end

  p
end

function Base.print(io::IO, p::ChainProgress)
  pct = iround(100.0 * p.counter / p.iters)
  str = string(
    "Chain ", p.chain, ": ", lpad(p.counter, length(string(p.iters)), ' '),
    "/$(p.iters) [", lpad(pct, 3, ' '), "%] @ $(strftime(time()))"
  )
  println(io, str)
end
