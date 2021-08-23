#
# Project : Pansy
# Source  : pwscf.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/08/23
#

#=
### *Abstract Types*
=#

"""
    PWInput

An abstract type representing an input object of quantum espresso. All
other input types should subtype `PWInput`.  It is used to build the
internal type system.
"""
abstract type PWInput end

"""
    PWInputEntry

Represent any component of an `PWInput`. The fields of an `PWInput`
should all be either `PWInputEntry` or `Nothing` (no value provided).
It is used to build the internal type system.
"""
abstract type PWInputEntry end

"""
    PWCard

Represent abstract cards of an `PWInput` in pwscf. It is used to build
the internal type system.
"""
abstract type PWCard <: PWInputEntry end

"""
    KPointsCard

Represent abstract ``k``-mesh or ``k``-path in pwscf.
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
        #
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
        #
        if eltype(shift) != Bool
            shift = Bool.(shift)
        end
        #
        return new(mesh, shift)
    end
end

"""
    MonkhorstPackGrid(k1, k2, k3, s1, s2, s3)

Constructor for `MonkhorstPackGrid`.

See also: [`ReciprocalPoint`](@ref).
"""
function MonkhorstPackGrid(k1, k2, k3, s1, s2, s3)
    k = [k1, k2, k3]
    s = [s1, s2, s3]
    return MonkhorstPackGrid(k, s)
end

#=
### *Customized Structs : Input Entry*
=#

"""
    AtomicSpecies

Represent each line of the `ATOMIC_SPECIES` card in pwscf. The `atom`
field accepts at most 3 characters.

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

Represent each line of the `ATOMIC_POSITIONS` card in pwscf. The `atom`
field accepts at most 3 characters.

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

Represent a component of an `PWInput`, a basic Fortran data structure.
It is used to build the internal type system.

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
    AtomicSpeciesCard

Represent the `ATOMIC_SPECIES` card in pwscf.

### Members

* data -> A vector containing `AtomicSpecies`.

See also: [`AtomicSpecies`](@ref).
"""
struct AtomicSpeciesCard <: Card
    data :: Vector{AtomicSpecies}
end

"""
    AtomicPositionsCard

Represent the `ATOMIC_POSITIONS` card in pwscf.

### Members

* data   -> A vector containing `AtomicPosition`s.
* option -> The scheme about how to define atomic positions.

See also: [`AtomicPosition`](@ref).
"""
struct AtomicPositionsCard <: Card
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

Represent the `K_POINTS` card in pwscf (`automatic` mode).

See also: [`KPointsCard`](@ref).
"""
struct AutoKmeshCard <: KPointsCard
    data :: MonkhorstPackGrid
end

"""
    GammaPointCard

Represent the `K_POINTS` card in pwscf (`gamma` mode).

See also: [`KPointsCard`](@ref).
"""
struct GammaPointCard <: KPointsCard end

"""
    SpecialKPointsCard

Represent the `K_POINTS` card in pwscf.

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
### *Parsers*
=#

"""
    tryparse(::Type{PWNamelist}, strs::Vector{String}, name::String)

Try to parse the `PWNamelist` object.

See also: [`PWNamelist`](@ref).
"""
function Base.tryparse(::Type{PWNamelist}, strs::Vector{String}, name::String)
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
    popfirst!(group_data)

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

function Base.tryparse(::Type{AutoKmeshCard}, str::AbstractString)
    m = match(K_POINTS_AUTOMATIC_BLOCK, str)

    if m !== nothing
        data = map(x -> parse(I64, x), m.captures)
        return AutoKmeshCard(MonkhorstPackGrid(data[1:3], data[4:6]))
    end
end

function Base.tryparse(::Type{GammaPointCard}, str::AbstractString)
    m = match(K_POINTS_GAMMA_BLOCK, str)
    return m === nothing ? nothing : GammaPointCard()
end

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

function Base.tryparse(::Type{KPointsCard}, str::AbstractString)
    for T in (AutoKmeshCard, GammaPointCard, SpecialPointsCard)
        x = tryparse(T, str)
        if x !== nothing
            return x
        end
    end
end

function Base.parse(::Type{T}, strs::Vector{String}) where {T <: Namelist}
    x = tryparse(T, strs)
    return T(x)
end

function Base.parse(::Type{T}, str::AbstractString) where {T<:Card}
    x = tryparse(T, str)
    if x === nothing
        throw(Meta.ParseError("cannot find card `$(block_name(T))`!"))
    else
        return x
    end
end

function Base.write(io::IO, x::T) where {T <: Namelist}
    println(io, " @$(block_name(T))")
    for key in keys(x.data)
        println(io, "    $key = ", x.data[key], ",")
    end
    println(io, " /")
end

function Base.write(io::IO, x::AtomicSpeciesCard)
    println(io, "ATOMIC_SPECIES")
    for i = 1:length(x.data)
        AS = x.data[i]
        println(io, " $(AS.atom)  $(AS.mass)  $(AS.upf)")
    end
end

function Base.write(io::IO, x::AtomicPositionsCard)
    println(io, "ATOMIC_POSITIONS {$(x.option)}")
    for i =  1:length(x.data)
        AP = x.data[i]
        print(io, " $(AP.atom) ")
        @printf(io, "%6.3f %6.3f %6.3f\n", AP.pos...)
    end
end

function Base.write(io::IO, x::AutoKmeshCard)
    println(io, "K_POINTS {automatic}")
    MPG = x.data
    @printf(io, "%3i%3i%3i%2i%2i%2i\n", MPG.mesh..., MPG.shift...)
end

function Base.write(io::IO, x::GammaPointCard)
    println(io, "K_POINTS {gamma}")
end

function Base.write(io::IO, x::SpecialPointsCard)
    println(io, "K_POINTS {$(x.option)}")
    nks = length(x.data)
    println(io, "  $nks")
    for i = 1:nks
        RP = x.data[i]
        @printf(io, " %11.7f%11.7f%11.7f%7.2f\n", RP.coord..., RP.weight)
    end
end

function Base.setindex!(nml::T, value, key::AbstractString) where {T <: Namelist}
    nml.data[key] = value
end

function Base.delete!(nml::T, key::AbstractString) where {T <: Namelist}
    delete!(nml.data, key)
end

#=
### *Driver Functions*
=#

function pwscf_adaptor()
end

"""
    pwscf_init(it::IterInfo)

Check the runtime environment of pwscf, prepare necessary input files.

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
    # Parse PWSCF.INP file, get the PWInput struct.
    PWINP = pwscfio_input()
    #
    # Create the real input file
    pwscfc_input(PWINP, it)
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
    pwscfc_input(PWINP::PWInput, it::IterInfo)
"""
function pwscfc_input(PWINP::PWInput, it::IterInfo)
    println(PWINP)

    # Customize your PWInput according to the case.toml
    #
    # For smearing
    smear = get_d("smear")
    @cswitch smear begin
        @case "mp2"
            PWINP.SystemNL["occupations"] = "'smearing'"
            PWINP.SystemNL["smearing"] = "'m-p'"
            break

        @case "mp1"
            PWINP.SystemNL["occupations"] = "'smearing'"
            PWINP.SystemNL["smearing"] = "'m-p'"
            break

        @case "gauss"
            PWINP.SystemNL["occupations"] = "'smearing'"
            PWINP.SystemNL["smearing"] = "'gauss'"
            break

        @case "tetra"
            PWINP.SystemNL["occupations"] = "'tetrahedra'"
            delete!(PWINP.SystemNL, "smearing")
            delete!(PWINP.SystemNL, "degauss")
            break

        @default
            PWINP.SystemNL["occupations"] = "'smearing'"
            PWINP.SystemNL["smearing"] = "'gauss'"
            break
    end

    case = get_c("case")
    finput = "$case.scf"
    open(finput, "w") do fout
        write(fout, PWINP.ControlNL)
        write(fout, PWINP.SystemNL)
        write(fout, PWINP.ElectronsNL)
        write(fout, PWINP.AtomicSpeciesBlock)
        write(fout, PWINP.AtomicPositionsBlock)
        write(fout, PWINP.KPointsBlock)
    end
end

#=
### *Service Functions* : *Group B*
=#

#=
### *Service Functions* : *Group C*
=#

"""
    pwscfio_input()

Parse the `PWSCF.INP` file, and return the PWInput struct. Actually,
`PWSCF.INP` is a standard, but mini input file for pwscf. It should
contains the `control`, `system`, `electrons` (namelists) and the
`ATOMIC_SPECIES`, `ATOMIC_POSITIONS`, `K_POINTS` (cards) input blocks
only. If you want to support more input entries, please make your
own modifications.

See also: [`PWInput`](@ref).
"""
function pwscfio_input()
    # Check the file status
    finput = "PWSCF.INP"
    @assert isfile(finput)

    # Parse the namelists
    lines = readlines(finput)
    ControlNL = parse(ControlNamelist, lines)
    SystemNL = parse(SystemNamelist, lines)
    ElectronsNL = parse(ElectronsNamelist, lines)

    # Parse the cards
    line = read(finput, String)
    AtomicSpeciesBlock = parse(AtomicSpeciesCard, line)
    AtomicPositionsBlock = parse(AtomicPositionsCard, line)
    KPointsBlock = parse(KPointsCard, line)

    # Return a PWInput struct
    return PWInput(ControlNL, 
                   SystemNL,
                   ElectronsNL,
                   AtomicSpeciesBlock,
                   AtomicPositionsBlock,
                   KPointsBlock)
end
