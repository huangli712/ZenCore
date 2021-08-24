#
# Project : Pansy
# Source  : pwscf.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/08/24
#

#=
### *Abstract Types*
=#

"""
    PWInputEntry

An abstract type representing an input component of `pwscf`. Note that
all other input types (such as `PWCard` and `PWNamelist`) should subtype
`PWInputEntry`.  It is used to build the internal type system.

See also: [`PWCard`](@ref), [`PWNamelist`](@ref).
"""
abstract type PWInputEntry end

"""
    PWCard

Represent abstract cards in the input file of `pwscf`.  It is used to
build the internal type system. The input file of `pwscf` consists of
various cards and namelists, represented by `PWCard` and `PWNamelist`,
respectively.

See also: [`PWNamelist`](@ref).
"""
abstract type PWCard <: PWInputEntry end

"""
    KPointsCard

Represent abstract `K-POINTS` card in the input file of `pwscf`.
"""
abstract type KPointsCard <: PWCard end

#=
### *Customized Structs : K-Grid*
=#

"""
    ReciprocalPoint

Represent a special point of the 3D Brillouin zone. Each of them has
a weight `w`.

### Members

* coord  -> Coordinates, i.e., ``k_x``, ``k_y``, and ``k_z``.
* weight -> Weight for the ``k``-point.

See also: [`MonkhorstPackGrid`](@ref).
"""
struct ReciprocalPoint
    coord  :: Vector{F64}
    weight :: F64

    # Inner constructor
    function ReciprocalPoint(coord, weight)
        @assert length(coord) == 3
        @assert weight > 0.0
        return new(coord, weight)
    end
end

"""
    ReciprocalPoint(x::F64, y::F64, z::F64, w::F64)

Constructor for `ReciprocalPoint`.

See also: [`MonkhorstPackGrid`](@ref).
"""
function ReciprocalPoint(x::F64, y::F64, z::F64, w::F64)
    return ReciprocalPoint([x, y, z], w)
end

"""
    MonkhorstPackGrid

Represent the Monkhorst-Pack grid.

### Members

* mesh  -> A length-three vector specifying the ``k``-point grid
           (``nk_1 × nk_2 × nk_3``) as in Monkhorst-Pack grids.
* shift -> A length-three vector specifying whether the grid is displaced
           by half a grid step in the corresponding directions.

See also: [`ReciprocalPoint`](@ref).
"""
struct MonkhorstPackGrid
    mesh  :: Vector{I64}
    shift :: Vector{Bool}

    # Inner constructor
    function MonkhorstPackGrid(mesh, shift)
        @assert length(mesh) == 3
        @assert length(shift) == 3
        @assert all(mesh .>= 1)
        if eltype(shift) != Bool
            shift = Bool.(shift)
        end
        return new(mesh, shift)
    end
end

"""
    MonkhorstPackGrid(k1::I64, k2::I64, k3::I64, s1::I64, s2::I64, s3::I64)

Constructor for `MonkhorstPackGrid`.

See also: [`ReciprocalPoint`](@ref).
"""
function MonkhorstPackGrid(k1::I64, k2::I64, k3::I64, s1::I64, s2::I64, s3::I64)
    k = [k1, k2, k3]
    s = [s1, s2, s3]
    return MonkhorstPackGrid(k, s)
end

#=
### *Customized Structs : Input Entry*
=#

"""
    AtomicSpecies

Represent each line of the `ATOMIC_SPECIES` card in the input file of
`pwscf`. The `atom` field accepts at most 3 characters.

### Members

* atom -> Label of the atom. Maximum total length cannot exceed
          3 characters.
* mass -> Mass of the atomic species in atomic unit. Used only when
          performing molecular dynamics (MD) run or structural
          optimization runs using damped MD.
* upf  -> File containing pseudopotential for this species.

### Examples

```julia-repl
julia> AtomicSpecies("C1", 12, "C.pbe-n-kjpaw_psl.1.0.0.UPF")
AtomicSpecies("C1", 12.0, "C.pbe-n-kjpaw_psl.1.0.0.UPF")
```

See also: [`AtomicSpeciesCard`](@ref).
"""
struct AtomicSpecies
    atom :: String
    mass :: F64
    upf  :: String

    # Inner constructor
    function AtomicSpecies(atom::Union{AbstractChar,AbstractString}, mass, upf)
        @assert length(atom) ≤ 3
        return new(string(atom), mass, upf)
    end
end

"""
    AtomicPosition

Represent each line of the `ATOMIC_POSITIONS` card in the input file of
`pwscf`. The `atom` field accepts at most 3 characters.

### Members

* atom   -> Label of the atom as specified in `AtomicSpecies`.
* pos    -> Atomic positions. A three-element vector of floats.
* if_pos -> Component `i` of the force for this atom is multiplied
            by `if_pos(i)`, which must be either `0` or `1`.  Used
            to keep selected atoms and/or selected components fixed
            in MD dynamics or structural optimization run.

### Examples

```julia-repl
julia> AtomicPosition('O', [0, 0, 0])
AtomicPosition("O", [0.0, 0.0, 0.0], Bool[1, 1, 1])
```

See also: [`AtomicPositionsCard`](@ref).
"""
struct AtomicPosition
    atom   :: String
    pos    :: Vector{F64}
    if_pos :: Vector{Bool}

    # Inner constructor
    function AtomicPosition(atom::Union{AbstractChar,AbstractString}, pos, if_pos)
        @assert length(atom) ≤ 3
        @assert length(pos) == 3
        @assert length(if_pos) == 3
        return new(string(atom), pos, if_pos)
    end
end

"""
    AtomicSpecies(x::AtomicPosition, mass, upf)

Constructor for `AtomicSpecies`.

### Examples

```julia-repl
julia> AtomicSpecies(
           AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764]),
           32.066,
           "S.pz-n-rrkjus_psl.0.1.UPF",
       )
AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
```

See also: [`AtomicSpeciesCard`](@ref).
"""
AtomicSpecies(x::AtomicPosition, mass, upf) = AtomicSpecies(x.atom, mass, upf)

"""
    AtomicPosition(atom, pos)

Constructors for `AtomicPosition`.

See also: [`AtomicPositionsCard`](@ref).
"""
AtomicPosition(atom, pos) = AtomicPosition(atom, pos, trues(3))

"""
    AtomicPosition(x::AtomicSpecies, pos, if_pos)

Constructors for `AtomicPosition`.

### Examples

```julia-repl
julia> AtomicPosition(
           AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
           [0.500000000, 0.288675130, 1.974192764],
       )
AtomicPosition("S", [0.5, 0.28867513, 1.974192764], Bool[1, 1, 1])
```

See also: [`AtomicPositionsCard`](@ref).
"""
AtomicPosition(x::AtomicSpecies, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)

#=
### *Customized Structs : Input Blocks*
=#

"""
    PWNamelist

Represent a namelist in the input file of `pwscf`, a basic Fortran
data structure. It is used to build the internal type system.

### Members

* name -> Name of the namelist. It should be `control`, `system`, or
          `electrons`. If you want to support more namelists, please
          make your own modifications.
* data -> A dict containing pairs of key and value.

See also: [`PWCard`](@ref).
"""
mutable struct PWNamelist <: PWInputEntry
    name :: AbstractString
    data :: Dict{AbstractString,Any}
end

"""
    Base.setindex!(pnl::PWNamelist, value, key::AbstractString)

Modify an entry (specified by `key`) in the namelist object (`pnl`).

See also: [`PWNamelist`](@ref).
"""
function Base.setindex!(pnl::PWNamelist, value, key::AbstractString)
    pnl.data[key] = value
end

"""
    Base.delete!(pnl::PWNamelist, key::AbstractString)

Remove an entry (specified by `key`) in the namelist object (`pnl`).

See also: [`PWNamelist`](@ref).
"""
function Base.delete!(pnl::PWNamelist, key::AbstractString)
    delete!(pnl.data, key)
end

"""
    AtomicSpeciesCard

Represent the `ATOMIC_SPECIES` card in the input file of `pwscf`.

### Members

* data -> A vector containing `AtomicSpecies`.

See also: [`AtomicSpecies`](@ref).
"""
struct AtomicSpeciesCard <: PWCard
    data :: Vector{AtomicSpecies}
end

"""
    AtomicPositionsCard

Represent the `ATOMIC_POSITIONS` card in the input file of `pwscf`.

### Members

* data   -> A vector containing `AtomicPosition`s.
* option -> The scheme about how to define atomic positions.

See also: [`AtomicPosition`](@ref).
"""
struct AtomicPositionsCard <: PWCard
    data   :: Vector{AtomicPosition}
    option :: String

    # Inner constructor
    function AtomicPositionsCard(data, option = "alat")
        @assert option in ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
        return new(data, option)
    end
end

"""
    AutoKmeshCard

Represent the `K_POINTS` card in the input file of `pwscf` (within
the `automatic` mode).

See also: [`KPointsCard`](@ref).
"""
struct AutoKmeshCard <: KPointsCard
    data :: MonkhorstPackGrid
end

"""
    GammaPointCard

Represent the `K_POINTS` card in the input file of `pwscf` (within
the `gamma` mode).

See also: [`KPointsCard`](@ref).
"""
struct GammaPointCard <: KPointsCard end

"""
    SpecialKPointsCard

Represent the `K_POINTS` card in the input file of `pwscf`.

### Members

* data   -> A vector containing `ReciprocalPoint`s.
* option -> The way about how to define ``k``-mesh.

See also: [`KPointsCard`](@ref).
"""
struct SpecialPointsCard <: KPointsCard
    data   :: Vector{ReciprocalPoint}
    option :: String

    # Inner constructor
    function SpecialPointsCard(data, option = "tpiba")
        @assert option in ("tpiba",
                           "crystal",
                           "tpiba_b",
                           "crystal_b",
                           "tpiba_c",
                           "crystal_c")
        return new(data, option)
    end
end

#=
### *Constants Regex*

### *Remarks* :

Note that these regular expressions are followed by various combinations
of the `i`, `m`, and `x` flags. These flags have the following meanings:

* `i` : Do case-insensitive pattern matching.
* `m` : Treat string as multiple lines.
* `s` : Treat string as single line.
* `x` : Tells the regular expression parser to ignore most whitespace
        that is neither backslashed nor within a character class.

See: https://docs.julialang.org/en/v1/manual/strings/#Regular-Expressions
=#

#=
### *Remarks* : For ATOMIC_SPECIES Card
=#

const ATOMIC_SPECIES_BLOCK = r"""
^ [ \t]* ATOMIC_SPECIES [ \t]* \R+
(?P<block>
    (?:
        ^ [ \t]* \S+ [ \t]+ (?:[-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*) [ \t]+ \S+ [ \t]* \R?
    )+
)
"""imx

const ATOMIC_SPECIES_ITEM = r"""
^ [ \t]* (?P<name>\S+) [ \t]+ (?P<mass>[-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*) [ \t]+ (?P<pseudo>\S+)
    [ \t]* \R?
"""mx

#=
### *Remarks* : For ATOMIC_POSITIONS Card
=#

const ATOMIC_POSITIONS_BLOCK = r"""
^ \s* ATOMIC_POSITIONS \s*                      # Atomic positions start with that string
[{(]? \s* (?P<units>\S+?)? \s* [)}]? \s* $\R    # The units are after the string in optional brackets
(?P<block>                                      # This is the block of positions
    (
        (
            \s*                                 # White space in front of the element spec is ok
            (
                [A-Za-z]+[A-Za-z0-9]{0,2}       # Element spec
                (
                    \s+                         # White space in front of the number
                    [-|+]?                      # Plus or minus in front of the number (optional)
                    (
                        (
                            \d*                 # optional decimal in the beginning .0001 is ok, for example
                            [\.]                # There has to be a dot followed by
                            \d+                 # at least one decimal
                        )
                        |                       # OR
                        (
                            \d+                 # at least one decimal, followed by
                            [\.]?               # an optional dot ( both 1 and 1. are fine)
                            \d*                 # And optional number of decimals (1.00001)
                        )                       # followed by optional decimals
                    )
                    ([E|e|d|D][+|-]?\d+)?       # optional exponents E+03, e-05
                ){3}                            # I expect three float values
                ((\s+[0-1]){3}\s*)?             # Followed by optional ifpos
                \s*                             # Followed by optional white space
                |
                \#.*                            # If a line is commented out, that is also ok
                |
                \!.*                            # Comments also with excl. mark in fortran
            )
            |                                   # OR
            \s*                                 # A line only containing white space
         )
        \R                                      # line break at the end
    )+                                          # A positions block should be one or more lines
)
"""imx

const ATOMIC_POSITIONS_ITEM = r"""
^                                       # Linestart
[ \t]*                                  # Optional white space
(?P<name>[A-Za-z]+[A-Za-z0-9]{0,2})\s+  # get the symbol, max 3 chars, starting with a char
(?P<x>                                  # Get x
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<y>                                  # Get y
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<z>                                  # Get z
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]*
(?P<fx>[01]?)                           # Get fx
[ \t]*
(?P<fy>[01]?)                           # Get fx
[ \t]*
(?P<fz>[01]?)                           # Get fx
"""mx

#=
### *Remarks* : For K_POINTS Card
=#

const K_POINTS_AUTOMATIC_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* automatic [ \t]* [)}]? [ \t]* \R
^ [ \t]* (\d+) [ \t]+ (\d+) [ \t]+ (\d+) [ \t]+ (\d+) [ \t]+ (\d+)
    [ \t]+ (\d+) [ \t]* \R?
"""imx

const K_POINTS_GAMMA_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* gamma [ \t]* [)}]? [ \t]* \R*
"""imx

const K_POINTS_SPECIAL_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]*
    [{(]? [ \t]* (?P<type>\S+?)? [ \t]* [)}]? [ \t]* \R+
^ [ \t]* \S+ [ \t]* \R+  # nks
(?P<block>
 (?:
  ^ [ \t]* \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]* \R*
 )+
)
"""imx

const K_POINTS_SPECIAL_ITEM = r"""
^ [ \t]* (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]* \R?
"""mx

#=
### *Input File Parsers*
=#

"""
    Base.parse(::Type{PWNamelist}, strs::Vector{String}, name::String)

Parse the `PWNamelist` object. The name of the namelist is specified
by argument `name`.

See also: [`PWNamelist`](@ref).
"""
function Base.parse(::Type{PWNamelist}, strs::Vector{String}, name::String)
    # Try to parse `strs` to extract the data
    #
    # Prepare necessary data structures
    group_data = []
    group_meet = false
    group_name = "unknown"
    #
    # Go through each line
    for l in eachindex(strs)
        # Get rid of the blanks and `,`
        strip_line = strip(strip(strs[l]), ',')
        # Meet a namelist
        if startswith(strip_line, "&")
            group_name = lowercase(split(strip_line, "&")[2])
            # It is the namelist what we try to locate
            if group_name == name
                group_meet = true
            end
        end

        # End of namelist
        if startswith(strip_line, "/")
            group_meet = false
        end

        # Store the data in group_data
        if group_meet
            push!(group_data, strip_line)
        end
    end
    #
    # The first element is not useful. We have to remove it.
    if isempty(group_data)
        return nothing # Fail to figure out a namelist
    else
        popfirst!(group_data)
    end

    # Try to build a PWNamelist object
    #
    # Prepare a dictionary, it will contain the namelist's data.
    NLData = Dict{AbstractString,Any}()
    #
    # Go throught each element in group_data
    for i in eachindex(group_data)
        # There are multiple entries in the element (line)
        if count(",", group_data[i]) > 0
            pairs = split(group_data[i], ",")
            for j in eachindex(pairs)
                key, value = map(x -> strip(x), split(pairs[j], "="))
                NLData[key] = value
            end
        # There is only one entry in the element (line)
        else
            key, value = map(x -> strip(x), split(group_data[i], "="))
            NLData[key] = value
        end
    end

    # Return the desired object
    return PWNamelist(name, NLData)
end

"""
    Base.parse(::Type{T}, str::AbstractString)

Parse the `PWCard` object. Now we support the following cards:

* `ATOMIC_SPECIES` (`AtomicSpeciesCard`)
* `ATOMIC_POSITIONS` (`AtomicPositionsCard`)
* `K_POINTS` (`AutoKmeshCard`, `GammaPointCard`, `SpecialKPointsCard`)

See also: [`PWCard`](@ref).
"""
function Base.parse(::Type{T}, str::AbstractString) where {T<:PWCard}
    x = tryparse(T, str)
    if x === nothing
        error("cannot find card `$T`!")
    else
        return x
    end
end

"""
    Base.tryparse(::Type{AtomicSpeciesCard}, str::AbstractString)

Try to parse the `AtomicSpeciesCard` object.

See also: [`AtomicSpeciesCard`](@ref).
"""
function Base.tryparse(::Type{AtomicSpeciesCard}, str::AbstractString)
    m = match(ATOMIC_SPECIES_BLOCK, str)

    # Function `match` only searches for the first match of the regular
    # expression, so it could be a `nothing`.
    if m !== nothing
        content = only(m.captures)
        return AtomicSpeciesCard(
            map(eachmatch(ATOMIC_SPECIES_ITEM, content)) do matched
                captured = matched.captures
                atom, mass, upf = captured[1], parse(F64, captured[2]), captured[3]
                AtomicSpecies(atom, mass, upf)
            end,
        )
    end
end

"""
    Base.tryparse(::Type{AtomicPositionsCard}, str::AbstractString)

Try to parse the `AtomicPositionsCard` object.

See also: [`AtomicPositionsCard`](@ref).
"""
function Base.tryparse(::Type{AtomicPositionsCard}, str::AbstractString)
    m = match(ATOMIC_POSITIONS_BLOCK, str)

    # Function `match` only searches for the first match of the regular
    # expression, so it could be a `nothing`.
    if m !== nothing
        if string(m.captures[1]) === nothing
            @warn "Not specifying units is DEPRECATED and will no longer be allowed in the future!"
            @info "No option is specified, 'alat' is assumed."
            option = "alat"
        else
            option = string(m.captures[1])
        end
        content = m.captures[2]
        return AtomicPositionsCard(
            map(eachmatch(ATOMIC_POSITIONS_ITEM, content)) do matched
                # The `matched` cannot be a `nothing` since we have
                # tested by the block regular expression.
                captured = matched.captures
                #
                # The `if_pos` field is optionally given by users. If
                # they do not give, we provide the default values `1`.
                if_pos = map(x -> isempty(x) ? 1 : parse(I64, x), captured[11:13])
                #
                # The `atom` and `pos` fields are mandatory. So we do
                # not need special treatment.
                atom, pos = captured[1],
                    map(x -> parse(F64, x), [captured[2], captured[5], captured[8]])
                AtomicPosition(atom, pos, if_pos)
            end,
            option,
        )
    end
end

"""
    Base.tryparse(::Type{AutoKmeshCard}, str::AbstractString)

Try to parse the `AutoKmeshCard` object.

See also: [`AutoKmeshCard`](@ref).
"""
function Base.tryparse(::Type{AutoKmeshCard}, str::AbstractString)
    m = match(K_POINTS_AUTOMATIC_BLOCK, str)

    if m !== nothing
        data = map(x -> parse(I64, x), m.captures)
        return AutoKmeshCard(MonkhorstPackGrid(data[1:3], data[4:6]))
    end
end

"""
    Base.tryparse(::Type{GammaPointCard}, str::AbstractString)

Try to parse the `GammaPointCard` object.

See also: [`GammaPointCard`](@ref).
"""
function Base.tryparse(::Type{GammaPointCard}, str::AbstractString)
    m = match(K_POINTS_GAMMA_BLOCK, str)
    return m === nothing ? nothing : GammaPointCard()
end

"""
    Base.tryparse(::Type{SpecialPointsCard}, str::AbstractString)

Try to parse the `SpecialPointsCard` object.

See also: [`SpecialPointsCard`](@ref).
"""
function Base.tryparse(::Type{SpecialPointsCard}, str::AbstractString)
    m = match(K_POINTS_SPECIAL_BLOCK, str)

    if m !== nothing
        option = m.captures[1] === nothing ? "tpiba" : m.captures[1]
        return SpecialPointsCard(
            map(eachmatch(K_POINTS_SPECIAL_ITEM, m.captures[2])) do matched
                # TODO: Match `nks`
                ReciprocalPoint(map(x -> parse(F64, x), matched.captures)...)
            end,
            option,
        )
    end
end

"""
    Base.tryparse(::Type{KPointsCard}, str::AbstractString)

Try to parse the `KPointsCard` object.

See also: [`KPointsCard`](@ref).
"""
function Base.tryparse(::Type{KPointsCard}, str::AbstractString)
    for T in (AutoKmeshCard, GammaPointCard, SpecialPointsCard)
        x = tryparse(T, str)
        if x !== nothing
            return x
        end
    end
end

#=
### *Input File Writers*
=#

"""
    Base.write(io::IO, x::PWNamelist)

Write the `PWNamelist` object to `IOStream`.

See also: [`PWNamelist`](@ref).
"""
function Base.write(io::IO, x::PWNamelist)
    println(io, " &$(x.name)")
    for key in keys(x.data)
        println(io, "    $key = ", x.data[key], ",")
    end
    println(io, " /")
end

"""
    Base.write(io::IO, x::AtomicSpeciesCard)

Write the `AtomicSpeciesCard` object to `IOStream`.

See also: [`AtomicSpeciesCard`](@ref).
"""
function Base.write(io::IO, x::AtomicSpeciesCard)
    println(io, "ATOMIC_SPECIES")
    for i = 1:length(x.data)
        AS = x.data[i]
        println(io, " $(AS.atom)  $(AS.mass)  $(AS.upf)")
    end
end

"""
    Base.write(io::IO, x::AtomicPositionsCard)

Write the `AtomicPositionsCard` object to `IOStream`.

See also: [`AtomicPositionsCard`](@ref).
"""
function Base.write(io::IO, x::AtomicPositionsCard)
    println(io, "ATOMIC_POSITIONS {$(x.option)}")
    for i =  1:length(x.data)
        AP = x.data[i]
        print(io, " $(AP.atom) ")
        @printf(io, "%6.3f %6.3f %6.3f\n", AP.pos...)
    end
end

"""
    Base.write(io::IO, x::AutoKmeshCard)

Write the `AutoKmeshCard` object to `IOStream`.

See also: [`AutoKmeshCard`](@ref).
"""
function Base.write(io::IO, x::AutoKmeshCard)
    println(io, "K_POINTS {automatic}")
    MPG = x.data
    @printf(io, "%3i%3i%3i%2i%2i%2i\n", MPG.mesh..., MPG.shift...)
end

"""
    Base.write(io::IO, x::GammaPointCard)

Write the `GammaPointCard` object to `IOStream`.

See also: [`GammaPointCard`](@ref).
"""
function Base.write(io::IO, x::GammaPointCard)
    println(io, "K_POINTS {gamma}")
end

"""
    Base.write(io::IO, x::SpecialPointsCard)

Write the `SpecialPointsCard` object to `IOStream`.

See also: [`SpecialPointsCard`](@ref).
"""
function Base.write(io::IO, x::SpecialPointsCard)
    println(io, "K_POINTS {$(x.option)}")
    nks = length(x.data)
    println(io, "  $nks")
    for i = 1:nks
        RP = x.data[i]
        @printf(io, " %11.7f%11.7f%11.7f%7.2f\n", RP.coord..., RP.weight)
    end
end

#=
### *Driver Functions*
=#

function pwscf_adaptor()
end

"""
    pwscf_init(it::IterInfo)

Check the runtime environment of `pwscf`, prepare necessary input files.

See also: [`pwscf_exec`](@ref), [`pwscf_save`](@ref).
"""
function pwscf_init(it::IterInfo)
    # Print the header
    println("Engine : PWSCF")
    println("Try to perform ab initio electronic structure calculation")
    println("Current directory: ", pwd())
    println("Prepare necessary input files for pwscf")

    # Prepare essential input files
    # Copy PWSCF.INP
    cp("../PWSCF.INP", joinpath(pwd(), "PWSCF.INP"), force = true)
    println("  > File PWSCF.INP is ready")
    #
    # Create the real input file, case.scf and case.nscf.
    pwscfc_input(it)
end

"""
    pwscf_exec(it::IterInfo)
"""
function pwscf_exec(it::IterInfo)
end

"""
    pwscf_save(it::IterInfo)
"""
function pwscf_save(it::IterInfo)
end

#=
### *Service Functions* : *Group A*
=#

"""
    pwscfc_input(it::IterInfo)

It will parse the `PWSCF.INP` file at first. Actually, `PWSCF.INP` is
a standard, but mini input file for `pwscf`. It only includes three
namelists (namely `control`, `system`, and `electrons`) and three
cards (namely `ATOMIC_SPECIES`, `ATOMIC_POSITIONS`, and `K_POINTS`).
If you want to support more input entries, please make your own
modifications here.

Then this function will try to customize these namelists and cards
according to the setup in `case.toml`.

Finally, this function will generate the input files for `pwscf`. They
are `case.scf` and `case.nscf`. As shown by their names, one is for the
self-consistent calculation, the other is for the non-self-consistent
calculation.

See also: [`PWNamelist`](@ref), [`PWCard`](@ref).
"""
function pwscfc_input(it::IterInfo)
    # Check the file status
    finput = "PWSCF.INP"
    @assert isfile(finput)

    # Parse the namelists, control, system, and electrons.
    lines = readlines(finput)
    ControlNL = parse(PWNamelist, lines, "control")
    SystemNL = parse(PWNamelist, lines, "system")
    ElectronsNL = parse(PWNamelist, lines, "electrons")

    # Parse the cards, ATOMIC_SPECIES, ATOMIC_POSITIONS, and K_POINTS.
    line = read(finput, String)
    AtomicSpeciesBlock = parse(AtomicSpeciesCard, line)
    AtomicPositionsBlock = parse(AtomicPositionsCard, line)
    KPointsBlock = parse(KPointsCard, line)

    # Customize the namelists and cards according to case.toml
    #
    # For smearing
    smear = get_d("smear")
    @cswitch smear begin
        @case "mp2"
            SystemNL["occupations"] = "'smearing'"
            SystemNL["smearing"] = "'m-p'"
            break

        @case "mp1"
            SystemNL["occupations"] = "'smearing'"
            SystemNL["smearing"] = "'m-p'"
            break

        @case "gauss"
            SystemNL["occupations"] = "'smearing'"
            SystemNL["smearing"] = "'gauss'"
            break

        @case "tetra"
            SystemNL["occupations"] = "'tetrahedra'"
            delete!(SystemNL, "smearing")
            delete!(SystemNL, "degauss")
            break

        @default
            SystemNL["occupations"] = "'smearing'"
            SystemNL["smearing"] = "'gauss'"
            break
    end

    # For kmesh density
    #
    # Note that kmesh == "file" is not supported for pwscf.
    kmesh = get_d("kmesh")
    if isa(KPointsBlock,AutoKmeshCard)
    @cswitch kmesh begin
        @case "accurate"
            write(ios, "KSPACING = 0.1 \n")
            break

        @case "medium"
            write(ios, "KSPACING = 0.2 \n")
            break

        @case "coarse"
            write(ios, "KSPACING = 0.4 \n")
            break

        @case "file"
            break

        @default # Very coarse kmesh
            write(ios, "KSPACING = 0.5 \n")
            break
    end
    end

    case = get_c("case")
    finput = "$case.scf"
    open(finput, "w") do fout
        write(fout, ControlNL)
        write(fout, SystemNL)
        write(fout, ElectronsNL)
        write(fout, AtomicSpeciesBlock)
        write(fout, AtomicPositionsBlock)
        write(fout, KPointsBlock)
    end
end

#=
### *Service Functions* : *Group B*
=#

#=
### *Service Functions* : *Group C*
=#
