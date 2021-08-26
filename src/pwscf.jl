#
# Project : Pansy
# Source  : pwscf.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/08/26
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
    Base.getindex(pnl::PWNamelist, key::AbstractString)

Return an entry (specified by `key`) in the namelist object (`pnl`).

See also: [`PWNamelist`](@ref).
"""
Base.getindex(pnl::PWNamelist, key::AbstractString) = pnl.data[key]

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
    AutoKmeshCard(mesh::Vector{I64}, shift::Vector{I64})

Constructor for `AutoKmeshCard`.

See also: [`KPointsCard`](@ref).
"""
function AutoKmeshCard(mesh::Vector{I64}, shift::Vector{Bool})
    return AutoKmeshCard(MonkhorstPackGrid(mesh, shift))
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

"""
    SpecialPointsCard(nkx::I64, nky::I64, nkz::I64, option::String = "crystal")

Constructor for `SpecialPointsCard`.

See also: [`KPointsCard`](@ref).
"""
function SpecialPointsCard(nkx::I64, nky::I64, nkz::I64, option::String = "crystal")
    # Sanity check
    @assert nkx >= 1
    @assert nky >= 1
    @assert nkz >= 1

    # Calculate total number of k-points
    nkpt = nkx * nky * nkz

    # Calculate weight per k-point
    w = 1.0 / nkpt

    # Generate uniform k-mesh. Note that the k-point is represented
    # by ReciprocalPoint struct.
    data = Vector{ReciprocalPoint}(undef,nkpt)
    c = 0 # Counter
    #
    for x = 0:nkx - 1
        for y = 0:nky - 1
            for z = 0:nkz - 1
                c = c + 1
                kx = float(x) / nkx
                ky = float(y) / nky
                kz = float(z) / nkz
                data[c] = ReciprocalPoint(kx, ky, kz, w)
            end
        end
    end
    #
    @assert c == nkpt

    # Call the default constructor
    return SpecialPointsCard(data, option)
end

"""
    SpecialPointsCard(nkx::I64, option::String = "crystal")

Constructor for `SpecialPointsCard`.

See also: [`KPointsCard`](@ref).
"""
function SpecialPointsCard(nkx::I64, option::String = "crystal")
    return SpecialPointsCard(nkx, nkx, nkx, option)
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
        @printf(io, " %11.7f%11.7f%11.7f%16.12f\n", RP.coord..., RP.weight)
    end
end

#=
### *Driver Functions*
=#

"""
    pwscf_adaptor(D::Dict{Symbol,Any})

Adaptor support for pwscf code. It will parse the output files of pwscf
code, extract the Kohn-Sham dataset, and then fulfill the `DFTData`
dict (i.e `D`).

The following pwscf's output files are needed:

* `scf.out`

Note that in the input file of pwscf, the verbosity parameter must be
set to 'high'.

See also: [`plo_adaptor`](@ref), [`ir_adaptor`](@ref).
"""
function pwscf_adaptor(D::Dict{Symbol,Any})
    # P01: Print the header
    println("Adaptor : PWSCF")
    println("Try to extract the Kohn-Sham dataset")
    println("Current directory: ", pwd())

    # P02: Read in lattice structure
    D[:latt] = pwscfio_lattice(pwd(), false)
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
    ControlNL, AtomicSpeciesBlock = pwscfc_input(it)
    case = get_c("case")
    println("  > File $case.scf is ready")
    println("  > File $case.nscf is ready")
    #
    # Check the pseudopotentials
    pdir = strip(ControlNL["pseudo_dir"],''')
    upf = map(x -> joinpath(pdir, x.upf), AtomicSpeciesBlock.data)
    for f in upf
        @assert isfile(f)
        println("  > File $f is ready")
    end
end

"""
    pwscf_exec(it::IterInfo, scf::Bool = true)

Execute the pwscf program, monitor the convergence progress, and output
the relevant information. The argument `scf` determines which input file
should be used. If scf == true, then the input file is case.scf, or else
it is case.nscf.

In order to execute this function correctly, you have to setup the
following environment variables:

* PWSCF_HOME

and make sure the file `MPI.toml` is available.

See also: [`pwscf_init`](@ref), [`pwscf_save`](@ref).
"""
function pwscf_exec(it::IterInfo, scf::Bool = true)
    # Print the header
    println("Detect the runtime environment for pwscf")

    # Determine mpi prefix (whether the pwscf is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "dft", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  > Using $numproc processors (MPI)")

    # Get the home directory of pwscf
    dft_home = query_dft("pwscf")
    println("  > Home directory for pwscf: ", dft_home)

    # Select suitable pwscf program
    # We use the same code pw.x for with or without spin-orbit coupling
    if get_d("lspinorb")
        pwscf_exe = "$dft_home/pw.x"
    else
        pwscf_exe = "$dft_home/pw.x"
    end
    @assert isfile(pwscf_exe)
    println("  > Executable program is available: ", basename(pwscf_exe))

    # Assemble command
    if isnothing(mpi_prefix)
        pwscf_cmd = pwscf_exe
    else
        pwscf_cmd = split("$mpi_prefix $pwscf_exe", " ")
    end
    println("  > Assemble command: $(prod(x -> x * ' ', pwscf_cmd))")

    # Determine suitable input and output files
    case = get_c("case")
    finp = "$case.scf"
    fout = "scf.out"
    if !scf
        finp = "$case.nscf"
        fout = "nscf.out"
    end
    println("  > Self-consistent DFT calculation: $scf")
    println("  > Using input file: $finp")
    println("  > Using output file: $fout")

    # Print the header
    println("Launch the computational engine pwscf")

    # Create a task, but do not run it immediately
    t = @task begin
        run(pipeline(`$pwscf_cmd`, stdin = finp, stdout = fout))
    end
    println("  > Create a task")

    # Launch it, the terminal output is redirected to `fout`.
    # Note that the task runs asynchronously. It will not block
    # the execution.
    schedule(t)
    println("  > Add the task to the scheduler's queue")
    println("  > Waiting ...")

    # To ensure that the task is executed
    while true
        sleep(2)
        istaskstarted(t) && break
    end

    # Analyze the `fout` file during the calculation
    #
    # `c` is a time counter
    c = 0
    #
    # Enter infinite loop
    while true
        # Sleep five seconds
        sleep(5)

        # Increase the counter
        c = c + 1

        # For self-consistent DFT calculation mode
        if scf

            # Parse the `fout` file
            iters = readlines(fout)
            filter!(x -> contains(x, "iteration #"), iters)
            ethrs = readlines(fout)
            filter!(x -> contains(x, "ethr ="), ethrs)

            # Figure out the number of iterations (`ni`) and deltaE (`dE`)
            if length(ethrs) > 0
                arr = line_to_array(iters[end])
                ni = parse(I64, arr[3])
                arr = line_to_array(ethrs[end])
                dE = strip(arr[3],',')
            else # The first iteration has not been finished
                ni = 0
                dE = "unknown"
            end

            # Print the log to screen
            @printf("  > Elapsed %4i seconds, %3i iterations (dE = %12s)\r", 5*c, ni, dE)

        # For non-self-consistent DFT calculation mode
        else

            # Parse the `fout` file
            lines = readlines(fout)
            filter!(x -> contains(x, "Computing kpt #"), lines)

            # Figure out how many k-points are finished
            if length(lines) > 0
                arr = line_to_array(lines[end])
                ikpt = parse(I64, arr[4])
                nkpt = parse(I64, arr[6])
            else # The first k-point has not been finished
                ikpt = 0
                nkpt = 0
            end

            # Print the log to screen
            @printf("  > Elapsed %4i seconds, %4i of %4i k-points\r", 5*c, ikpt, nkpt)

        end

        # Break the loop
        istaskdone(t) && break
    end
    #
    # Keep the last output
    println()

    # Wait for the pwscf task to finish
    wait(t)

    # Extract how many iterations are executed
    if scf
        iters = readlines(fout)
        filter!(x -> contains(x, "iteration #"), iters)
        println("  > Converged after $(length(iters)) iterations")
    # Extract how many k-points are finished
    else
        lines = readlines(fout)
        filter!(x -> contains(x, "Computing kpt #"), lines)
        println("  > Calculated eigenvalues for $(length(lines)) k-points")
    end
end

"""
    pwscf_save(it::IterInfo)

Backup the output files of pwscf if necessary. Furthermore, the DFT fermi
level in `IterInfo` struct is also updated (`IterInfo.μ₀`).

See also: [`pwscf_init`](@ref), [`pwscf_exec`](@ref).
"""
function pwscf_save(it::IterInfo)
    # Print the header
    println("Finalize the computational task")

    # Store the data files
    #
    # Create list of files
    case = get_c("case")
    fl = ["$case.scf", "$case.nscf", "scf.out", "nscf.out"]
    #
    # Go through the file list, backup the files one by one.
    for i in eachindex(fl)
        f = fl[i]
        cp(f, "$f.$(it.I₃)", force = true)
    end
    println("  > Save the key output files")

    # Anyway, the DFT fermi level is extracted from scf.out, and its
    # value will be saved at IterInfo.μ₀.
    it.μ₀ = pwscfio_fermi(pwd())
    println("  > Extract the fermi level from scf.out: $(it.μ₀) eV")

    # We also try to read the DFT band energy from scf.out, and its
    # value will be saved at IterInfo.et.
    it.et.dft = pwscfio_energy(pwd())
    println("  > Extract the DFT band energy from scf.out: $(it.et.dft) eV")
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

The return values of this function are namelist (`control`) and card
(`ATOMIC_SPECIES`), which will be used to check the pseudopotentials.

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
            break

        @default
            SystemNL["occupations"] = "'smearing'"
            SystemNL["smearing"] = "'gauss'"
            break
    end

    # For kmesh density
    #
    # Note that if kmesh == "file", the original setup in PWSCF.INP is
    # kept. In other words, KPointsBlock will not be changed.
    kmesh = get_d("kmesh")
    if isa(KPointsBlock,AutoKmeshCard)
        shift = copy(KPointsBlock.data.shift)
        @cswitch kmesh begin
            @case "accurate"
                KPointsBlock = AutoKmeshCard([10, 10, 10], shift)
                break

            @case "medium"
                KPointsBlock = AutoKmeshCard([08, 08, 08], shift)
                break

            @case "coarse"
                KPointsBlock = AutoKmeshCard([06, 06, 06], shift)
                break

            @case "file"
                break

            @default # Very coarse kmesh
                KPointsBlock = AutoKmeshCard([04, 04, 04], shift)
                break
        end
    end

    # For magnetic moment
    magmom = get_d("magmom")
    if !isa(magmom, Missing)
        SystemNL["starting_magnetization"] = magmom
    end

    # For symmetry
    lsymm = get_d("lsymm")
    if lsymm
        SystemNL["nosym"] = ".false."
    else # Ignore the symmetry completely
        SystemNL["nosym"] = ".true."
    end

    # For spin polarizations
    lspins = get_d("lspins")
    if lspins
        SystemNL["nspin"] = 2
    else
        SystemNL["nspin"] = 1
    end

    # For spin-orbit coupling
    lspinorb = get_d("lspinorb")
    if lspinorb
        SystemNL["noncolin"] = ".true."
        if lspins
            SystemNL["nspin"] = 4
        end
    else
        SystemNL["noncolin"] = ".false."
    end

    # For optimized projectors
    # SKIP

    # For local orbitals and projectors
    # SKIP

    # For number of bands
    # SKIP

    # Special treatment for verbosity
    ControlNL["verbosity"] = "'high'"

    # Special treatment for pseudo_dir
    pseudo_dir = strip(ControlNL["pseudo_dir"],''') # Get rid of `
    pseudo_dir = joinpath("..", pseudo_dir)
    ControlNL["pseudo_dir"] = "'$pseudo_dir'" # Add ' back

    # Build input files for pwscf
    #
    # Get case's name
    case = get_c("case")
    #
    # Setup filenames
    fscf = "$case.scf"
    fnscf = "$case.nscf"
    #
    # For case.scf
    ControlNL["calculation"] = "'scf'"
    open(fscf, "w") do fout
        write(fout, ControlNL)
        write(fout, SystemNL)
        write(fout, ElectronsNL)
        write(fout, AtomicSpeciesBlock)
        write(fout, AtomicPositionsBlock)
        write(fout, KPointsBlock)
    end
    #
    # For case.nscf
    ControlNL["calculation"] = "'nscf'"
    delete!(ControlNL, "restart_mode")
    #
    # We have to specify the k-points explicitly during the
    # non-self-consistent calculations.
    begin
        # At the same time, the tetrahedron algorithm will fail.
        SystemNL["occupations"] = "'smearing'"
        @cswitch kmesh begin
            @case "accurate"
                KPointsBlock = SpecialPointsCard(10)
                break

            @case "medium"
                KPointsBlock = SpecialPointsCard(08)
                break

            @case "coarse"
                KPointsBlock = SpecialPointsCard(06)
                break

            @case "file"
                @assert isa(KPointsBlock, SpecialPointsCard)
                break

            @default # Very coarse kmesh
                KPointsBlock = SpecialPointsCard(04)
                break
        end
    end
    open(fnscf, "w") do fout
        write(fout, ControlNL)
        write(fout, SystemNL)
        write(fout, ElectronsNL)
        write(fout, AtomicSpeciesBlock)
        write(fout, AtomicPositionsBlock)
        write(fout, KPointsBlock)
    end

    # Return the namelist and the card, which will be used to check
    # whether the pseudopotential files are ready.
    return ControlNL, AtomicSpeciesBlock
end

#=
### *Service Functions* : *Group B*
=#

"""
    pwscfq_files(f::String)

Check the essential output files by pwscf. Here `f` means only the
directory that contains the desired files.

See also: [`adaptor_run`](@ref).
"""
function pwscfq_files(f::String)
    fl = ["scf.out", "nscf.out"]
    for i in eachindex(fl)
        @assert isfile( joinpath(f, fl[i]) )
    end
end

"""
    pwscfq_files()

Check the essential output files by pwscf in the current directory.

See also: [`adaptor_run`](@ref).
"""
pwscfq_files() = pwscfq_files(pwd())

#=
### *Service Functions* : *Group C*
=#

"""
    pwscfio_energy(f::String)

Reading pwscf's `scf.out` file, return DFT total energy, which will
be used to determine the DFT + DMFT energy. Here `f` means only the
directory that contains `scf.out`.
"""
function pwscfio_energy(f::String)
    # Try to figure out whether the scf.out file is valid
    lines = readlines(joinpath(f, "scf.out"))
    filter!(x -> contains(x, "!    total energy"), lines)
    @assert length(lines) == 1

    # Extract the total energy
    etot = parse(F64, line_to_array(lines[end])[5])

    # Return the desired value
    return etot
end

"""
    pwscfio_energy()

Reading pwscf's `scf.out` file, return DFT total energy, which will
be used to determine the DFT + DMFT energy.
"""
pwscfio_energy() = pwscfio_energy(pwd())

"""
    vaspio_lattice(f::String, silent::Bool = true)

Reading vasp's `POSCAR` file, return crystallography information. Here `f`
means only the directory that contains `POSCAR`.

See also: [`Lattice`](@ref), [`irio_lattice`](@ref).
"""
function vaspio_lattice(f::String, silent::Bool = true)
    # Print the header
    !silent && println("Parse lattice")
    !silent && println("  > Open and read POSCAR")

    # Open the iostream
    fin = open(joinpath(f, "POSCAR"), "r")

    # Get the case
    _case = string(strip(readline(fin)))

    # Get the scaling factor
    scale = parse(F64, readline(fin))

    # Get the lattice vectors
    lvect = zeros(F64, 3, 3)
    lvect[1, :] = parse.(F64, line_to_array(fin))
    lvect[2, :] = parse.(F64, line_to_array(fin))
    lvect[3, :] = parse.(F64, line_to_array(fin))

    # Get the symbol list
    symbols = line_to_array(fin)

    # Get the number of sorts of atoms
    nsort = length(symbols)

    # Get the number list
    numbers = parse.(I64, line_to_array(fin))

    # Get the total number of atoms
    natom = sum(numbers)

    # Now all the parameters are ready, we would like to create
    # `Lattice` struct here.
    latt = Lattice(_case, scale, nsort, natom)

    # Update latt using the available data
    latt.lvect = lvect
    for i = 1:nsort
        latt.sorts[i, 1] = string(symbols[i])
        latt.sorts[i, 2] = numbers[i]
    end

    # Get the atom list
    k = 0
    for i = 1:nsort
        for j = 1:numbers[i]
            k = k + 1
            latt.atoms[k] = symbols[i]
        end
    end
    # Sanity check
    @assert k == natom

    # Get the coordinates of atoms
    readline(fin)
    for i = 1:natom
        latt.coord[i, :] = parse.(F64, line_to_array(fin)[1:3])
    end

    # Close the iostream
    close(fin)

    # Print some useful information to check
    !silent && println("  > System: ", latt._case)
    !silent && println("  > Atoms: ", latt.atoms)

    # Return the desired struct
    return latt
end

"""
    pwscfio_lattice()

Reading pwscf's `scf.out` file, return crystallography information.

See also: [`Lattice`](@ref), [`irio_lattice`](@ref).
"""
pwscfio_lattice() = pwscfio_lattice(pwd())

"""
    pwscfio_fermi(f::String, silent::Bool = true)

Reading pwscf's `scf.out` file, return the fermi level. Here `f` means
only the directory that contains `scf.out`.

See also: [`irio_fermi`](@ref).
"""
function pwscfio_fermi(f::String, silent::Bool = true)
    # Print the header
    !silent && println("Parse fermi level")
    !silent && println("  > Open and read scf.out")

    # Try to figure out whether the scf.out file is valid
    lines = readlines(joinpath(f, "scf.out"))
    filter!(x -> contains(x, "Fermi energy"), lines)
    @assert length(lines) == 1

    # Extract the fermi level
    fermi = parse(F64, line_to_array(lines[end])[5])

    # Print some useful information to check
    !silent && println("  > Fermi level: $fermi eV")

    # Return the desired data
    return fermi
end

"""
    pwscfio_fermi()

Reading pwscf's `scf.out` file, return the fermi level.

See also: [`irio_fermi`](@ref).
"""
pwscfio_fermi() = pwscfio_fermi(pwd())
