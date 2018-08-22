module UnitSystems

import Base.*, Base./, Base.^, Base.show
export Unit, *, show

struct Unit{T}
    unit::T
end

function *{T}(u::Unit{T},v::Unit{T})::Unit{T}
    return Unit(T(u.unit.base, u.unit.exps + v.unit.exps))
end

function /{T}(u::Unit{T},v::Unit{T})::Unit{T}
    return Unit(T(u.unit.base, u.unit.exps - v.unit.exps))
end

function ^{T}(u::Unit{T},n::Int)::Unit{T}
    return Unit(T(u.unit.base, n * u.unit.exps))
end


supscripts = Dict("1"=>'\U00B9',"2"=>'\U00B2',"3"=>'\U00B3',
                  "4"=>'\U2074',"5"=>'\U2075',"6"=>'\U2076',
                  "7"=>'\U2077',"8"=>'\U2078',"9"=>'\U2079',
                  "0"=>'\U2070',"-"=>'\U207B')

function supscriptInt(n::Int)::String
    sn = string(n)
    asn = [supscripts[string(x)] for x in sn]
    return string(reduce(*,asn))
end

function showExp(s::Symbol, n::Int)::String
    n==0 && return string()
    n==1 && return string(s)
    return string(s,supscriptInt(n))
end

function show{T}(io::IO, a::Unit{T})
    sp(n) = if (n!=0) string(" ") else string() end
    l = length(a.unit.base)
    x = [(a.unit.base[i],a.unit.exps[i]) for i in 1:l]
    for s in x print(io,showExp(s[1],s[2])*sp(s[2])) end
    return Void()
    #print(io, "$(a.unit.base[1])^$(a.unit.exps[1])")
    #print(io, "$(a.unit.base[2])^$(a.unit.exps[2])")
    #Unicode .... print('\U00B2'), print('m','\U00B3')
    # +- : print('m','\U00B1')
    # https://www.fileformat.info/info/unicode/char/2070/index.htm
    # https://www.fileformat.info/info/unicode/char/b2/index.htm
    # print('m','\U2070') - print('m','\U2079')
end #show

end #module

module SIunits
importall UnitSystems
export SI, si

SIbase = [:m,:kg,:s,:A,:K,:mol,:cd]
SIexps = [0 for i in 1:length(SIbase)]

struct SI
    base::Array{Symbol}
    exps::Array{Int}
end

# Base units
const u1 = Unit(SI(SIbase,SIexps))
const m   = metre    = Unit(SI(SIbase,[1,0,0,0,0,0,0]))
const kg  = kilogram = Unit(SI(SIbase,[0,1,0,0,0,0,0]))
const s   = second   = Unit(SI(SIbase,[0,0,1,0,0,0,0]))
const A   = ampere   = Unit(SI(SIbase,[0,0,0,1,0,0,0]))
const K   = kelvin   = Unit(SI(SIbase,[0,0,0,0,1,0,0]))
const mol = mole     = Unit(SI(SIbase,[0,0,0,0,0,1,0]))
const cd  = candela  = Unit(SI(SIbase,[0,0,0,0,0,0,1]))

si = Dict([(:m,m),(:kg,kg),(:s,s),(:A,A),(:K,K),(:mol,mol),(:cd,cd)])

macro m(n) return(metre^n) end
end #module

# Howto: import Siunits ...
# Required: Pkg.add("IntervalArithmetic")

module PhysQty
include("units_doc.jl")
importall UnitSystems
import Base.*, Base.show
import SIunits.SI, SIunits.m
using IntervalArithmetic
export *

@doc PQTY_doc ->
struct PQTY{T,S}
    value::Interval{S}
    unit::Unit{T}
end

function *{T}(v::Real,u::Unit{T})::PQTY{T}
    return PQTY(interval(v),u)
end

function*{T,S}(v::Interval{S},u::Unit{T})::PQTY{T}
    return PQTY(v,u)
end

#=
function show{T}(io::IO, q::PQTY{T})
    pm = '\U00B1'
    sp = " "
    print(io,q.val)
    print(io,sp*pm*sp)
    print(io,q.err)
    print(io,sp^2)
    print(io,q.unit)
    return Void()
end #show
=#

IntervalArithmetic.setformat(:midpoint)

function show{T}(io::IO, q::PQTY{T})
    print(io,"$(q.value) $(q.unit)")
    return Void()
end #show


end #module


# import SIunits.si
# IntervalArithmetic.setformat(:midpoint,sigfigs=4)
# import IntervalArithmetic.±
# Alt-0177 Windows
# https://en.wikipedia.org/wiki/Plus-minus_sign
# (6.777±0.5)*si[:kg]



# Example
using SIunits
using PhysQty
a=3.5*si[:m]
typeof(a)
import IntervalArithmetic.±
b=(3.5±0.1)*(si[:m]^2*si[:kg]^(-5))
typeof(b)
Q=IntervalArithmetic.@biginterval(10.0)*(si[:m]^2*si[:kg]^(-5))
typeof(Q)
Q.value.hi
Q.unit
QQ = (Q.value^3)*Q.unit
sin(Q.value)
IntervalArithmetic.interval(8)*Q.unit
