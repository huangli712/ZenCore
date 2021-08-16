#
# Project : Pansy
# Source  : pwscf.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/08/17
#

function pwscf_adaptor()
end

"""
    Input

An abstract type representing an input object of ab initio software.
All other input types should subtype `Input`.
"""
abstract type Input end

"""
    InputEntry

Represent any component of an `Input`. The fields of an `Input` should
all be either `InputEntry` or `Nothing` (no value provided).
"""
abstract type InputEntry end

"""
    Namelist <: InputEntry

Represent a component of an `Input`, a basic Fortran data structure.
"""
abstract type Namelist <: InputEntry end

"""
    Card <: InputEntry

Represent cards of an `Input` in Quantum ESPRESSO.
"""
abstract type Card <: InputEntry end