#=
struct UnitSystem
    symbols::Vector{String}
    dims::Vector{Int}
    UnitSystem(s,d) = length(s)!=length(d) ? error("length") : new(s,d)
end

struct SI{T<:UnitSystem}
    name::String
    unit::T
end

#struct USystem(T<:)


struct BaseQuantity
    name::String
    qsym::String
    dsym::String
end

struct QuantitySystem{T<:Array{BaseQuantity}}
    name::String
    base::T
end

struct BaseUnit
    name::String
    usym::String
end

struct UnitSystem{T<:Dict{BaseUnit,BaseQuantity}}
    name::String
    base::T
end
=#

#=
struct SI{}
    base::Array{String}
    exps::Array{Int}
end

struct UnitSystem{T}
    sys::T
end

import Base.*
function *(u::UnitSystem,v::UnitSystem)::UnitSystem

sis=["m","kg","s"]
=#

import Base.*

SIbase=["m","kg","s","A"]
SIexps=[0,0,0,0]

struct SI
    base::Array{String}
    exps::Array{Int}
end

struct Unit{T}
    unit::T
end

function *(u::Unit,v::Unit)::Unit
    C = typeof(u.unit)
    return Unit(C(u.unit.base,u.unit.exps+v.unit.exps))
end

function test{T}(x::Unit{T})::Unit{T}
    return Unit(T(x.unit.base,x.unit.exps))
end

function *{T}(u::Unit{T},v::Unit{T})::Unit{T}
    return Unit(T(u.unit.base,u.unit.exps+v.unit.exps))
end
