#
# Project : Pansy
# Source  : pwscf.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/08/19
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

        #println(strip_line, " ", group_name, " ", group_start, " ", group_end)

        if group_start
            push!(group_data, strip_line)
        elseif !isempty(group_data)
            Namelists[Symbol(group_name)] = copy(group_data)
            empty!(group_data)
        end
    end

    #println(Namelists)

    #for key in keys(Namelists)
    #    println(typeof(Namelists[key]), typeof(_CONTROL))
    #end

    ControlNL = process_namelists!(Namelists[:control], _CONTROL)
    SystemNL = process_namelists!(Namelists[:system], _SYSTEM)
    ElectronsNL = process_namelists!(Namelists[:electrons], _ELECTRONS)
end

function process_namelists!(nml::Vector{Any}, valid::Tuple)
    NLData = Dict{Symbol,Any}()

    popfirst!(nml)
    for i in eachindex(nml)
        if count(",", nml[i]) > 0
            pairs = split(nml[i], ",")
            for j in eachindex(pairs)
                key, value = split(pairs[j], "=")
                NLData[Symbol(strip(key))] = strip(value)
            end
        else
            key, value = split(nml[i], "=")
            NLData[Symbol(strip(key))] = strip(value)
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
