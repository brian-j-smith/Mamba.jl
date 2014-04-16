#################### Variate Numerical Operators ####################

const ArithMethods2 = [
  :(Base.(:(+))),
  :(Base.(:(.+))),
  :(Base.(:(-))),
  :(Base.(:(.-))),
  :(Base.(:(*))),
  :(Base.(:(.*))),
  :(Base.(:(/))),
  :(Base.(:(./))),
  :(Base.(:(\))),
  :(Base.(:(.\))),
  :(Base.(:(^))),
  :(Base.(:(.^))),
  :(Base.(:(%))),
  :(Base.(:(.%))),
  :(Base.dot)
]

Base.(:^)(v::Variate, n::Integer) = Base.(:^)(v.data, n)

const ArrayMethods = [
  :(Base.cummax),
  :(Base.cummin),
  :(Base.cumprod),
  :(Base.cumsum),
  :(Base.cumsum_kbn),
  :(Base.diff),
  :(Base.maximum),
  :(Base.minimum),
  :(Base.norm),
  :(Base.prod),
  :(Base.sum),
  :(Base.sum_kbn)
]

const CompareMethods2 = [
  :(Base.(:(==))),
  :(Base.(:(.==))),
  :(Base.(:(!=))),
  :(Base.(:(.!=))),
  :(Base.(:(<))),
  :(Base.(:(.<))),
  :(Base.(:(<=))),
  :(Base.(:(.<=))),
  :(Base.(:(>))),
  :(Base.(:(.>))),
  :(Base.(:(>=))),
  :(Base.(:(.>=)))
]

const DivideMethods2 = [
  :(Base.div),
  :(Base.divrem),
  :(Base.fld),
  :(Base.mod),
  :(Base.rem)
]

const MathMethods = [
  :(Base.airy),
  :(Base.airyprime),
  :(Base.airyai),
  :(Base.airyaiprime),
  :(Base.airybi),
  :(Base.airybiprime),
  :(Base.besselh),
  :(Base.besselj0),
  :(Base.besselj1),
  :(Base.bessely0),
  :(Base.bessely1),
  :(Base.beta),
  :(Base.dawson),
  :(Base.digamma),
  :(Base.erf),
  :(Base.erfc),
  :(Base.erfcinv),
  :(Base.erfcx),
  :(Base.erfi),
  :(Base.erfinv),
  :(Base.erf),
  :(Base.erfc),
  :(Base.eta),
  :(Base.gamma),
  :(Base.lbeta),
  :(Base.lfact),
  :(Base.lgamma),
  :(Base.zeta)
]

const MathMethods2 = [
  :(Base.airy),
  :(Base.besselh),
  :(Base.besseli),
  :(Base.besselj),
  :(Base.besselk),
  :(Base.bessely),
  :(Base.hankelh1),
  :(Base.hankelh2)
]

const PowerMethods = [
  :(Base.sqrt),
  :(Base.cbrt),
  :(Base.exp),
  :(Base.expm1),
  :(Base.exponent),
  :(Base.log),
  :(Base.log10),
  :(Base.log1p),
  :(Base.log2),
  :(Base.significand)
]

const PowerMethods2 = [
  :(Base.hypot)
]

Base.ldexp(v::Variate, n::Integer) = Base.ldexp(v.data, n)

const RoundMethods = [
  :(Base.ceil),
  :(Base.iceil),
  :(Base.floor),
  :(Base.ifloor),
  :(Base.round),
  :(Base.iround),
  :(Base.trunc),
  :(Base.itrunc)
]

const RoundMethods2 = [
  :(Base.ceil),
  :(Base.floor),
  :(Base.round),
  :(Base.trunc)
]

const SignMethods = [
  :(Base.abs),
  :(Base.abs2),
  :(Base.sign),
  :(Base.signbit)
]

const SignMethods2 = [
  :(Base.copysign),
  :(Base.flipsign)
]

const StatMethods = [
  :(Base.cor),
  :(Base.cov),
  :(Base.mean),
  :(Base.median),
  :(Base.std),
  :(Base.var),
  :(StatsBase.kurtosis),
  :(StatsBase.mad),
  :(StatsBase.skewness)
]

const StatMethods2 = [
  :(Base.cor),
  :(Base.cov),
  :(StatsBase.corspearman)
]

const TrigMethods = [
  :(Base.acos),
  :(Base.acosh),
  :(Base.acot),
  :(Base.acoth),
  :(Base.acsc),
  :(Base.acsch),
  :(Base.asec),
  :(Base.asech),
  :(Base.asin),
  :(Base.asinh),
  :(Base.atan),
  :(Base.atanh),
  :(Base.cos),
  :(Base.cosc),
  :(Base.cosh),
  :(Base.cospi),
  :(Base.cot),
  :(Base.coth),
  :(Base.csc),
  :(Base.csch),
  :(Base.sec),
  :(Base.sech),
  :(Base.sin),
  :(Base.sinc),
  :(Base.sinh),
  :(Base.sinpi),
  :(Base.tan),
  :(Base.tanh)
]

const TrigMethods2 = [
  :(Base.atan2)
]

const UnaryMethods = [
  :(Base.(:+)),
  :(Base.(:-)),
  :(Base.isfinite),
  :(Base.isinf),
  :(Base.isnan),
  :(Base.ctranspose),
  :(Base.transpose)
]

for op in [ArithMethods2, CompareMethods2, DivideMethods2, MathMethods2,
           PowerMethods2, SignMethods2, StatMethods2, TrigMethods2]
  @eval begin
    ($op)(x::Union(Array, Number), v::Variate) = ($op)(x, v.data)
    ($op)(v::Variate, x::Union(Array, Number)) = ($op)(v.data, x)
    ($op)(u::Variate, v::Variate) = ($op)(u.data, v.data)
  end
end

for op in [ArrayMethods, MathMethods, PowerMethods, RoundMethods,
           SignMethods, StatMethods, TrigMethods, UnaryMethods]
  @eval ($op)(v::Variate) = ($op)(v.data)
end

for op in RoundMethods2
  @eval ($op)(v::Variate, digits::Integer) = ($op)(v.data, digits)
end
