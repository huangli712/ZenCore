#
# Project : Pansy
# Source  : pwscf.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/08/21
#


function pwscf_adaptor()
end

function pwscf_parser()
    lines = readlines("diamond.scf")
    ControlNL = parse(ControlNamelist, lines)
    SystemNL = parse(SystemNamelist, lines)
    ElectronsNL = parse(ElectronsNamelist, lines)

    str = read("diamond.scf", String)
    AtomicSpeciesBlock = parse(AtomicSpeciesCard, str)
    AtomicPositionsBlock = parse(AtomicPositionsCard, str)
    KPointsBlock = parse(KPointsCard, str)

    open("case.scf", "w") do fout
        write(fout, ControlNL)
        write(fout, SystemNL)
        write(fout, ElectronsNL)
        write(fout, AtomicSpeciesBlock)
    end
    #println(AtomicSpeciesBlock)

    return PWInput(ControlNL, SystemNL, ElectronsNL, AtomicSpeciesBlock, AtomicPositionsBlock, KPointsBlock)
end

#=
### *Abstract Types*
=#

"""
    Input

An abstract type representing an input object of ab initio software.
All other input types should subtype `Input`.  It is used to build
the internal type system.
"""
abstract type Input end

"""
    InputEntry

Represent any component of an `Input`. The fields of an `Input` should
all be either `InputEntry` or `Nothing` (no value provided). It is used
to build the internal type system.
"""
abstract type InputEntry end

"""
    Namelist

Represent a component of an `Input`, a basic Fortran data structure. It
is used to build the internal type system.
"""
abstract type Namelist <: InputEntry end

"""
    Card

Represent abstract cards of an `Input` in Quantum ESPRESSO. It is used
to build the internal type system.
"""
abstract type Card <: InputEntry end

"""
    KPointsCard

Represent abstract ``k``-mesh or ``k``-path in Quantum ESPRESSO.
"""
abstract type KPointsCard <: Card end

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

#=
### *Customized Structs : Input Entry*
=#

"""
    AtomicSpecies

Represent each line of the `ATOMIC_SPECIES` card in Quantum ESPRESSO.
The `atom` field accepts at most 3 characters.

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

Represent each line of the `ATOMIC_POSITIONS` card in Quantum ESPRESSO.
The `atom` field accepts at most 3 characters.

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
AtomicPosition(atom, pos) = AtomicPosition(atom, pos, trues(3))
AtomicPosition(x::AtomicSpecies, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)

#=
### *Customized Structs : Input Blocks*
=#

"""
    ControlNamelist

Represent the `control` namelist in Quantum ESPRESSO.

### Members

* data -> A dict containing pairs of key and value.

See also: [`Namelist`](@ref).
"""
mutable struct ControlNamelist <: Namelist
    data :: Dict{AbstractString,Any}
end

"""
    SystemNamelist

Represent the `control` namelist in Quantum ESPRESSO.

### Members

* data -> A dict containing pairs of key and value.

See also: [`Namelist`](@ref).
"""
mutable struct SystemNamelist <: Namelist
    data :: Dict{AbstractString,Any}
end

"""
    ElectronsNamelist

Represent the `control` namelist in Quantum ESPRESSO.

### Members

* data -> A dict containing pairs of key and value.

See also: [`Namelist`](@ref).
"""
mutable struct ElectronsNamelist <: Namelist
    data :: Dict{AbstractString,Any}
end

block_name(::Type{ControlNamelist}) = "control"
block_name(::Type{SystemNamelist}) = "system"
block_name(::Type{ElectronsNamelist}) = "electrons"
block_vars(::Type{ControlNamelist}) = VAR_CONTROL
block_vars(::Type{SystemNamelist}) = VAR_SYSTEM
block_vars(::Type{ElectronsNamelist}) = VAR_ELECTRONS

"""
    AtomicSpeciesCard

Represent the `ATOMIC_SPECIES` card in Quantum ESPRESSO.

### Members

* data -> A vector containing `AtomicSpecies`.

See also: [`AtomicSpecies`](@ref).
"""
struct AtomicSpeciesCard <: Card
    data :: Vector{AtomicSpecies}
end

"""
    AtomicPositionsCard

Represent the `ATOMIC_POSITIONS` card in Quantum ESPRESSO.

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

Represent the `K_POINTS` card in Quantum ESPRESSO (`automatic` mode).

See also: [`KPointsCard`](@ref).
"""
struct AutoKmeshCard <: KPointsCard
    data :: MonkhorstPackGrid
end

"""
    GammaPointCard

Represent the `K_POINTS` card in Quantum ESPRESSO (`gamma` mode).

See also: [`KPointsCard`](@ref).
"""
struct GammaPointCard <: KPointsCard end

"""
    SpecialKPointsCard

Represent the `K_POINTS` card in Quantum ESPRESSO.

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
    PWInput

Represent the input file of Quantum ESPRESSO.
"""
mutable struct PWInput <: Input
    ControlNL   :: ControlNamelist
    SystemNL    :: SystemNamelist
    ElectronsNL :: ElectronsNamelist
    #
    AtomicSpeciesBlock   :: AtomicSpeciesCard
    AtomicPositionsBlock :: AtomicPositionsCard
    KPointsBlock         :: KPointsCard
end

#=
### *Constants Tuples*
=#

"""
    VAR_CONTROL

Represent the `CONTROL` namelist of `pwscf`.
"""
const VAR_CONTROL = (
    :calculation     ,
    :title           ,
    :verbosity       ,
    :restart_mode    ,
    :wf_collect      ,
    :nstep           ,
    :iprint          ,
    :tstress         ,
    :tprnfor         ,
    :dt              ,
    :outdir          ,
    :wfcdir          ,
    :prefix          ,
    :lkpoint_dir     ,
    :max_seconds     ,
    :etot_conv_thr64 ,
    :forc_conv_thr64 ,
    :disk_io         ,
    :pseudo_dir      ,
    :tefield         ,
    :dipfield        ,
    :lelfield        ,
    :nberrycyc       ,
    :lorbm           ,
    :lberry          ,
    :gdir            ,
    :nppstr          ,
    :lfcp            ,
    :gate
)

"""
    VAR_SYSTEM

Represent the `SYSTEM` namelist of `pwscf`.
"""
const VAR_SYSTEM = (
    :ibrav                     ,
    :celldm                    ,
    :A                         ,
    :B                         ,
    :C                         ,
    :cosAB                     ,
    :cosAC                     ,
    :cosBC                     ,
    :nat                       ,
    :ntyp                      ,
    :nbnd                      ,
    :tot_charge                ,
    :starting_charge           ,
    :tot_magnetization         ,
    :starting_magnetization    ,
    :ecutwfc                   ,
    :ecutrho                   ,
    :ecutfock                  ,
    :nr1                       ,
    :nr2                       ,
    :nr3                       ,
    :nr1s                      ,
    :nr2s                      ,
    :nr3s                      ,
    :nosym                     ,
    :nosym_evc                 ,
    :noinv                     ,
    :no_t_rev                  ,
    :force_symmorphic          ,
    :use_all_frac              ,
    :occupations               ,
    :one_atom_occupations      ,
    :starting_spin_angle       ,
    :degauss                   ,
    :smearing                  ,
    :nspin                     ,
    :noncolin                  ,
    :ecfixed                   ,
    :qcutz                     ,
    :q2sigma                   ,
    :input_dft                 ,
    :ace                       ,
    :exx_fraction              ,
    :screening_parameter       ,
    :exxdiv_treatment          ,
    :x_gamma_extrapolation     ,
    :ecutvcut                  ,
    :nqx1                      ,
    :nqx2                      ,
    :nqx3                      ,
    :localization_thr          ,
    :lda_plus_u                ,
    :lda_plus_u_kind           ,
    :Hubbard_U                 ,
    :Hubbard_J0                ,
    :Hubbard_alpha             ,
    :Hubbard_beta              ,
    :Hubbard_J                 ,
    :starting_ns_eigenvalue    ,
    :U_projection_type         ,
    :Hubbard_parameters        ,
    :ensemble_energies         ,
    :edir                      ,
    :emaxpos                   ,
    :eopreg                    ,
    :eamp                      ,
    :angle1                    ,
    :angle2                    ,
    :lforcet                   ,
    :constrained_magnetization ,
    :fixed_magnetization       ,
    :lambda                    ,
    :report                    ,
    :lspinorb                  ,
    :assume_isolated           ,
    :esm_bc                    ,
    :esm_w                     ,
    :esm_efield                ,
    :esm_nfit                  ,
    :lgcscf                    ,
    :gcscf_mu                  ,
    :gcscf_conv_thr            ,
    :gcscf_beta                ,
    :vdw_corr                  ,
    :london                    ,
    :london_s6                 ,
    :london_c6                 ,
    :london_rvdw               ,
    :london_rcut               ,
    :dftd3_version             ,
    :dftd3_threebody           ,
    :ts_vdw_econv_thr          ,
    :ts_vdw_isolated           ,
    :xdm                       ,
    :xdm_a1                    ,
    :xdm_a2                    ,
    :space_group               ,
    :uniqueb                   ,
    :origin_choice             ,
    :rhombohedral              ,
    :zgate                     ,
    :relaxz                    ,
    :block                     ,
    :block_1                   ,
    :block_2                   ,
    :block_height
)

"""
    VAR_ELECTRONS

Represent the `ELECTRONS` namelist of `pwscf`.
"""
const VAR_ELECTRONS = (
    :electron_maxstep  ,
    :scf_must_converge ,
    :conv_thr          ,
    :adaptive_thr      ,
    :conv_thr_init     ,
    :conv_thr_multi    ,
    :mixing_mode       ,
    :mixing_beta       ,
    :mixing_ndim       ,
    :mixing_fixed_ns   ,
    :diagonalization   ,
    :diago_thr_init    ,
    :diago_cg_maxiter  ,
    :diago_david_ndim  ,
    :diago_full_acc    ,
    :efield            ,
    :efield_cart       ,
    :efield_phase      ,
    :startingpot       ,
    :startingwfc       ,
    :tqr               ,
    :real_space
)

#=
### *Constants Regex*
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

function namelists(nml::Vector{Any}, vars::Tuple)
    NLData = Dict{AbstractString,Any}()

    for i in eachindex(nml)
        if count(",", nml[i]) > 0
            pairs = split(nml[i], ",")
            for j in eachindex(pairs)
                key, value = map(x -> strip(x), split(pairs[j], "="))
                ind = findfirst('(', key)
                if ind isa Nothing
                    @assert Symbol(key) in vars
                else
                    subkey = SubString(key, 1:ind-1)
                    @assert Symbol(subkey) in vars
                end
                NLData[key] = value
            end
        else
            key, value = map(x -> strip(x), split(nml[i], "="))
            ind = findfirst('(', key)
            if ind isa Nothing
                @assert Symbol(key) in vars
            else
                subkey = SubString(key, 1:ind-1)
                @assert Symbol(subkey) in vars
            end
            NLData[key] = value
        end
    end

    return NLData
end

function Base.tryparse(::Type{T}, strs::Vector{String}) where {T <: Namelist}
    group_data = []
    group_meet = false
    group_name = "unknown"

    for l in eachindex(strs)
        strip_line = strip(strip(strs[l]), ',')
        if startswith(strip_line, "&")
            group_name = lowercase(split(strip_line, "&")[2])
            if group_name == block_name(T)
                group_meet = true
            end
        end

        if startswith(strip_line, "/")
            group_meet = false
        end

        if group_meet
            push!(group_data, strip_line)
        end
    end

    popfirst!(group_data)

    return namelists(group_data, block_vars(T))
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
end