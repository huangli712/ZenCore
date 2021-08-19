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
const ElectronsNamelist = (
    electron_maxstep::UInt
    scf_must_converge::Bool
    conv_thr::Float64
    adaptive_thr::Bool
    conv_thr_init::Float64
    conv_thr_multi::Float64
    mixing_mode::String
    mixing_beta::Float64
    mixing_ndim::UInt
    mixing_fixed_ns::UInt
    diagonalization::String
    ortho_para::UInt
    diago_thr_init::Float64
    diago_cg_maxiter::UInt
    diago_david_ndim::UInt
    diago_full_acc::Bool
    efield::Float64
    efield_cart::Vector{Maybe{Float64}}
    efield_phase::String
    startingpot::String  # This depends on `calculation`
    startingwfc::String  # This depends on `calculation`
    tqr::Bool
)
