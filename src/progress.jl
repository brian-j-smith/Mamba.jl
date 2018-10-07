#################### ChainProgress Meter ####################

#################### Types and Constructors ####################

struct ChainProgressFrame
  verbose::Bool

  function ChainProgressFrame(title::AbstractString, verbose::Bool)
    verbose && print(title * "...\n\n")
    new(verbose)
  end
end

mutable struct ChainProgress
  frame::ChainProgressFrame
  chain::Int
  iters::Int
  counter::Int
  runin::Int
  threshold::Float64
  t0::Float64

  function ChainProgress(frame::ChainProgressFrame, chain::Integer,
                         iters::Integer)
    new(frame, chain, iters, 0, max(1, min(10, round(Int, 0.01 * iters))),
        0.0, time())
  end
end


#################### Base Methods ####################

function reset!(p::ChainProgress)
  p.counter = 0
  p.threshold = 0.0
  p.t0 = time()
  p
end

function next!(p::ChainProgress)
  p.counter += 1
  if p.counter / p.iters >= p.threshold && p.counter >= p.runin
    p.threshold += 0.10
    p.frame.verbose && print(stdout, p)
  end
  p
end

function Base.print(io::IO, p::ChainProgress)
  elapsed = time() - p.t0
  remaining = elapsed * (p.iters / p.counter - 1.0)
  str = @sprintf("Chain %u: %3u%% [%s of %s remaining]%s",
                 p.chain,
                 100.0 * p.counter / p.iters,
                 strfsec(remaining),
                 strfsec(elapsed + remaining),
                 "\n"^(p.counter == p.iters))
  println(io, str)
end

function strfsec(sec::Real)
  hr, sec = divrem(sec, 3600)
  min, sec = divrem(sec, 60)
  @sprintf "%u:%02u:%02u" hr min sec
end
