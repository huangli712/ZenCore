#
# Project : Pansy
# Source  : global.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Stable
#
# Last modified: 2025/05/10
#

#=
### *Global Constants* : *Numerical Types*
=#

"""
    I32

Alias of Integer type (32 bit).

See also: [`R32`](@ref), [`R64`](@ref), [`N32`](@ref), [`N64`](@ref).
"""
const I32 = Int32

"""
    I64

Alias of Integer type (64 bit).

See also: [`R32`](@ref), [`R64`](@ref), [`N32`](@ref), [`N64`](@ref).
"""
const I64 = Int64

"""
    F32

Alias of Float type (32 bit).

See also: [`R32`](@ref), [`R64`](@ref), [`N32`](@ref), [`N64`](@ref).
"""
const F32 = Float32

"""
    F64

Alias of Float type (64 bit).

See also: [`R32`](@ref), [`R64`](@ref), [`N32`](@ref), [`N64`](@ref).
"""
const F64 = Float64

"""
    C32

Alias of Complex type (32 bit).

See also: [`N32`](@ref), [`N64`](@ref).
"""
const C32 = ComplexF32

"""
    C64

Alias of Complex type (64 bit).

See also: [`N32`](@ref), [`N64`](@ref).
"""
const C64 = ComplexF64

"""
    R32

Alias of Integer and Float types (32 bit). Here `R` means Real.

See also: [`N32`](@ref), [`N64`](@ref).
"""
const R32 = Union{I32,F32}

"""
    R64

Alias of Integer and Float types (64 bit). Here `R` means Real.

See also: [`N32`](@ref), [`N64`](@ref).
"""
const R64 = Union{I64,F64}

"""
    N32

Alias of Integer, Float, and Complex types (32 bit). Here `N` means Number.

See also: [`R32`](@ref), [`R64`](@ref).
"""
const N32 = Union{I32,F32,C32}

"""
    N64

Alias of Integer, Float, and Complex types (64 bit). Here `N` means Number.

See also: [`R32`](@ref), [`R64`](@ref).
"""
const N64 = Union{I64,F64,C64}

#=
### *Global Constants* : *Literal Strings*
=#

"""
    __LIBNAME__

Name of this julia package.

See also: [`__VERSION__`](@ref).
"""
const __LIBNAME__ = "ZenCore"

"""
    __VERSION__

Version of this julia package.

See also: [`__RELEASE__`](@ref).
"""
const __VERSION__ = v"0.8.6-devel.250510"

"""
    __RELEASE__

Release date of this julia package.

See also: [`__AUTHORS__`](@ref).
"""
const __RELEASE__ = "2025/05"

#=
*Remarks* :

The elements of the Array `__AUTHORS__` should be a `NamedTuple` object,
such as:

```julia
(name = "author's name", email = "author's email")
```
=#

"""
    __AUTHORS__

Core authors of this julia package.

See also: [`__LIBNAME__`](@ref).
"""
const __AUTHORS__ = [(name = "Li Huang", email = "huangli@caep.cn")]

"""
    authors()

Print authors / contributors of the `ZenCore` package.

See also: [`__AUTHORS__`](@ref).
"""
function authors()
    println("Authors (Until $__RELEASE__):")
    for a in __AUTHORS__
        println("  $(a.name) (email: $(a.email))")
    end
end
