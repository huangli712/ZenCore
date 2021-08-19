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

    println(ControlNL)
    println(SystemNL)
    println(ElectronsNL)
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


"""
    AtomicSpecies(atom::Union{AbstractChar,String}, mass::Float64, pseudopot::String)
    AtomicSpecies(x::AtomicPosition, mass, pseudopot)
Represent each line of the `ATOMIC_SPECIES` card in QE.
The `atom` field accepts at most 3 characters.
# Examples
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

abstract type Card end

"""
    AtomicSpeciesCard <: Card
Represent the `ATOMIC_SPECIES` card in QE. It does not have an "option".
"""
struct AtomicSpeciesCard <: Card
    data::Vector{AtomicSpecies}
end

export AtomicSpeciesCard
export AtomicSpecies
function Base.tryparse(::Type{AtomicSpeciesCard}, str::AbstractString)
    m = match(ATOMIC_SPECIES_BLOCK, str)
    # Function `match` only searches for the first match of the regular expression, so it could be a `nothing`
    if m !== nothing
        content = only(m.captures)
        return AtomicSpeciesCard(
            map(eachmatch(ATOMIC_SPECIES_ITEM, content)) do matched
                captured = matched.captures
                println(captured)
                atom, mass, pseudopotential =
                    captured[1], parse(Float64, captured[2]), captured[3]
                AtomicSpecies(atom, mass, pseudopotential)
            end,
        )
    end
end # function Base.tryparse
