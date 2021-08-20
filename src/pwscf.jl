#
# Project : Pansy
# Source  : pwscf.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/08/20
#


function pwscf_adaptor()
end

function pwscf_parser()
    Namelists = Dict{Symbol,Any}()
    group_data = []
    group_start = false
    group_end = false
    group_name = "unknown"
    for line in eachline("diamond.scf")
        strip_line = strip(strip(line), ',')
        if startswith(strip_line, "&")
            group_name = lowercase(split(strip_line, "&")[2])
            group_start = true
            group_end = false
        end
        if startswith(strip_line, "/")
            group_start = false
            group_end = true
        end

        if group_start
            push!(group_data, strip_line)
        elseif !isempty(group_data)
            Namelists[Symbol(group_name)] = copy(group_data)
            empty!(group_data)
        end
    end

    ControlNL = process_namelists!(Namelists[:control], _CONTROL)
    SystemNL = process_namelists!(Namelists[:system], _SYSTEM)
    ElectronsNL = process_namelists!(Namelists[:electrons], _ELECTRONS)

    str = read("diamond.scf", String)
    AtomicSpeciesBlock = parse(AtomicSpeciesCard, str)
    AtomicPositionsBlock = parse(AtomicPositionsCard, str)
    KPointsBlock = parse(KPointsCard, str)

    return PWInput(ControlNL, SystemNL, ElectronsNL, AtomicSpeciesBlock, AtomicPositionsBlock, KPointsBlock)
end

function process_namelists!(nml::Vector{Any}, keylist::Tuple)
    popfirst!(nml)
    NLData = Dict{AbstractString,Any}()
    for i in eachindex(nml)
        if count(",", nml[i]) > 0
            pairs = split(nml[i], ",")
            for j in eachindex(pairs)
                key, value = map(x -> strip(x), split(pairs[j], "="))
                ind = findfirst('(', key)
                if ind isa Nothing
                    @assert Symbol(key) in keylist
                else
                    subkey = SubString(key, 1:ind-1)
                    @assert Symbol(subkey) in keylist
                end
                NLData[key] = value
            end
        else
            key, value = map(x -> strip(x), split(nml[i], "="))
            ind = findfirst('(', key)
            if ind isa Nothing
                @assert Symbol(key) in keylist
            else
                subkey = SubString(key, 1:ind-1)
                @assert Symbol(subkey) in keylist
            end
            NLData[key] = value
        end
    end

    return NLData
end


#=
### *Customized Structs*
=#

"""
    ReciprocalPoint

Represent a special point of the 3D Brillouin zone. Each of them has
a weight `w`.

See also: [`MonkhorstPackGrid`](@ref).
"""
struct ReciprocalPoint{T}
    coord  :: Vector{T}
    weight :: F64
end

"""
    MonkhorstPackGrid(mesh, is_shift)

Represent the Monkhorst--Pack grid.

### Arguments
- `mesh`: A length-three vector specifying the k-point grid (``nk_1 × nk_2 × nk_3``) as in Monkhorst--Pack grids.
- `is_shift`: A length-three vector specifying whether the grid is displaced by half a grid step in the corresponding directions.
"""
struct MonkhorstPackGrid
    mesh::Vector{I64}
    is_shift::Vector{Bool}
    function MonkhorstPackGrid(mesh, is_shift)
        @assert all(mesh .>= 1)
        if eltype(is_shift) != Bool
            is_shift = Bool.(is_shift)
        end
        return new(mesh, is_shift)
    end
end

"""
    AtomicSpecies(atom::Union{AbstractChar,String}, mass::Float64, pseudopot::String)
    AtomicSpecies(x::AtomicPosition, mass, pseudopot)

Represent each line of the `ATOMIC_SPECIES` card in QE.
The `atom` field accepts at most 3 characters.

### Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.PWscf
julia> AtomicSpecies("C1", 12, "C.pbe-n-kjpaw_psl.1.0.0.UPF")
AtomicSpecies("C1", 12.0, "C.pbe-n-kjpaw_psl.1.0.0.UPF")
julia> AtomicSpecies(
           AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764]),
           32.066,
           "S.pz-n-rrkjus_psl.0.1.UPF",
       )
AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
```
"""
struct AtomicSpecies
    "Label of the atom. Max total length cannot exceed 3 characters."
    atom :: String
    """
    Mass of the atomic species in atomic unit.
    Used only when performing molecular dynamics (MD) run
    or structural optimization runs using damped MD.
    Not actually used in all other cases (but stored
    in data files, so phonon calculations will use
    these values unless other values are provided).
    """
    mass :: Float64
    """
    File containing pseudopotential for this species.
    See also: [`pseudoformat`](@ref)
    """
    pseudopot::String
    function AtomicSpecies(atom::Union{AbstractChar,AbstractString}, mass, pseudopot)
        @assert length(atom) <= 3 "`atom` can have at most 3 characters!"
        return new(string(atom), mass, pseudopot)
    end
end

"""
    AtomicPosition(atom::Union{AbstractChar,String}, pos::Vector{Float64}[, if_pos::Vector{Int}])
    AtomicPosition(x::AtomicSpecies, pos, if_pos)

Represent each line of the `ATOMIC_POSITIONS` card in QE.
The `atom` field accepts at most 3 characters.

### Examples
```jldoctest
julia> using QuantumESPRESSOBase.Cards.PWscf
julia> AtomicPosition('O', [0, 0, 0])
AtomicPosition("O", [0.0, 0.0, 0.0], Bool[1, 1, 1])
julia> AtomicPosition(
           AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
           [0.500000000, 0.288675130, 1.974192764],
       )
AtomicPosition("S", [0.5, 0.28867513, 1.974192764], Bool[1, 1, 1])
```
"""
struct AtomicPosition
    "Label of the atom as specified in `AtomicSpecies`."
    atom::String
    "Atomic positions. A three-element vector of floats."
    pos::Vector{Float64}
    """
    Component `i` of the force for this atom is multiplied by `if_pos(i)`,
    which must be either `0` or `1`.  Used to keep selected atoms and/or
    selected components fixed in MD dynamics or structural optimization run.
    With `crystal_sg` atomic coordinates the constraints are copied in all equivalent
    atoms.
    """
    if_pos::Vector{Bool}
    function AtomicPosition(atom::Union{AbstractChar,AbstractString}, pos, if_pos)
        @assert length(atom) <= 3 "`atom` can have at most 3 characters!"
        return new(string(atom), pos, if_pos)
    end
end

abstract type Card end

"""
    AtomicSpeciesCard <: Card
Represent the `ATOMIC_SPECIES` card in QE. It does not have an "option".
"""
struct AtomicSpeciesCard <: Card
    data::Vector{AtomicSpecies}
end

"""
    AtomicPositionsCard <: Card

Represent the `ATOMIC_POSITIONS` card in QE.

### Arguments
- `data::AbstractVector{AtomicPosition}`: A vector containing `AtomicPosition`s.
- `option::String="alat"`: allowed values are: "alat", "bohr", "angstrom", "crystal", and "crystal_sg".
"""
struct AtomicPositionsCard <: Card
    data::Vector{AtomicPosition}
    option::String
    function AtomicPositionsCard(data, option = "alat")
        @assert option in optionpool(AtomicPositionsCard)
        return new(data, option)
    end
end

abstract type KPointsCard <: Card end

struct KMeshCard <: KPointsCard
    data::MonkhorstPackGrid
end

struct GammaPointCard <: KPointsCard end

"""
    SpecialKPointsCard(data, option)
Represent the `K_POINTS` card in QE.
# Arguments
- `data::Union{MonkhorstPackGrid,GammaPoint,AbstractVector{SpecialKPoint}}`: A Γ point, a Monkhorst--Pack grid or a vector containing `SpecialKPoint`s.
- `option::String="tpiba"`: allowed values are: "tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c" and "crystal_c".
"""
struct SpecialPointsCard <: KPointsCard
    data::Vector{ReciprocalPoint}
    option::String
    function SpecialPointsCard(data, option = "tpiba")
        @assert option in optionpool(SpecialPointsCard)
        return new(data, option)
    end
end

mutable struct PWInput
    ControlNL
    SystemNL
    ElectronsNL
    AtomicSpeciesBlock
    AtomicPositionsBlock
    KPointsBlock
end

"""
    _CONTROL

Represent the `CONTROL` namelist of `pwscf`.
"""
const _CONTROL = (
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
    _SYSTEM

Represent the `SYSTEM` namelist of `pwscf`.
"""
const _SYSTEM = (
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
    _ELECTRONS

Represent the `ELECTRONS` namelist of `pwscf`.
"""
const _ELECTRONS = (
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

# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/aedee19/qe_tools/parsers/_input_base.py
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

# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/aedee19/qe_tools/parsers/_input_base.py
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

# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_AUTOMATIC_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* automatic [ \t]* [)}]? [ \t]* \R
^ [ \t]* (\d+) [ \t]+ (\d+) [ \t]+ (\d+) [ \t]+ (\d+) [ \t]+ (\d+)
    [ \t]+ (\d+) [ \t]* \R?
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_GAMMA_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* gamma [ \t]* [)}]? [ \t]* \R*
"""imx
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
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
# This regular expression is taken from https://github.com/aiidateam/qe-tools/blob/develop/qe_tools/parsers/pwinputparser.py
const K_POINTS_SPECIAL_ITEM = r"""
^ [ \t]* (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]* \R?
"""mx

function Base.tryparse(::Type{AtomicSpeciesCard}, str::AbstractString)
    m = match(ATOMIC_SPECIES_BLOCK, str)
    # Function `match` only searches for the first match of the regular expression, so it could be a `nothing`
    if m !== nothing
        content = only(m.captures)
        return AtomicSpeciesCard(
            map(eachmatch(ATOMIC_SPECIES_ITEM, content)) do matched
                captured = matched.captures
                #println(captured)
                atom, mass, pseudopotential =
                    captured[1], parse(Float64, captured[2]), captured[3]
                AtomicSpecies(atom, mass, pseudopotential)
            end,
        )
    end
end # function Base.tryparse

function Base.tryparse(::Type{AtomicPositionsCard}, str::AbstractString)
    m = match(ATOMIC_POSITIONS_BLOCK, str)
    # Function `match` only searches for the first match of the regular expression, so it could be a `nothing`
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
                # The `matched` cannot be a `nothing` since we have tested by the block regular expression
                captured = matched.captures
                # The `if_pos` field is optionally given by users. If they do not give, we provide the default values `1`.
                if_pos = map(x -> isempty(x) ? 1 : parse(Int, x), captured[11:13])
                # The `atom` and `pos` fields are mandatory. So we do not need special treatment.
                atom, pos = captured[1],
                map(x -> parse(Float64, x), [captured[2], captured[5], captured[8]])
                AtomicPosition(atom, pos, if_pos)
            end,
            option,
        )
    end
end # function Base.tryparse

function Base.tryparse(::Type{KMeshCard}, str::AbstractString)
    m = match(K_POINTS_AUTOMATIC_BLOCK, str)
    if m !== nothing
        data = map(x -> parse(Int, x), m.captures)
        return KMeshCard(MonkhorstPackGrid(data[1:3], data[4:6]))
    end
end # function Base.tryparse

function Base.tryparse(::Type{GammaPointCard}, str::AbstractString)
    m = match(K_POINTS_GAMMA_BLOCK, str)
    return m === nothing ? nothing : GammaPointCard()
end # function Base.tryparse

function Base.tryparse(::Type{SpecialPointsCard}, str::AbstractString)
    m = match(K_POINTS_SPECIAL_BLOCK, str)
    if m !== nothing
        option = m.captures[1] === nothing ? "tpiba" : m.captures[1]
        return SpecialPointsCard(
            map(eachmatch(K_POINTS_SPECIAL_ITEM, m.captures[2])) do matched
                # TODO: Match `nks`
                ReciprocalPoint(map(x -> fparse(Float64, x), matched.captures)...)
            end,
            option,
        )
    end
end # function Base.tryparse

function Base.tryparse(::Type{KPointsCard}, str::AbstractString)
    for T in (GammaPointCard, KMeshCard, SpecialPointsCard)
        x = tryparse(T, str)
        if x !== nothing
            return x
        end
    end
end # function Base.tryparse

function Base.parse(::Type{T}, str::AbstractString) where {T<:Card}
    x = tryparse(T, str)
    if x === nothing
        throw(Meta.ParseError("cannot find card `$(groupname(T))`!"))
    else
        return x
    end
end # function Base.parse

optionpool(::Type{KMeshCard}) = ("automatic",)
optionpool(::Type{GammaPointCard}) = ("gamma",)
optionpool(::Type{SpecialPointsCard}) = ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
optionpool(::Type{AtomicPositionsCard}) = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")

AtomicPosition(atom, pos) = AtomicPosition(atom, pos, trues(3))
AtomicPosition(x::AtomicSpecies, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)
# Introudce mutual constructors since they share the same atoms.
AtomicSpecies(x::AtomicPosition, mass, pseudopot) = AtomicSpecies(x.atom, mass, pseudopot)