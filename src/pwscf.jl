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
    :Hubbard_U                 r{Maybe{Float64}}
    :Hubbard_J0                r{Maybe{Float64}}
    :Hubbard_alpha             r{Maybe{Float64}}
    :Hubbard_beta              r{Maybe{Float64}}
    :Hubbard_J                 x{Maybe{Float64}}
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
    SystemNamelist(; kwargs...)
"""
function SystemNamelist(;
    ibrav = 127,
    celldm = zeros(6),  # Must specify
    A = 0.0,
    B = 0.0,
    C = 0.0,
    cosAB = 0.0,
    cosAC = 0.0,
    cosBC = 0.0,
    nat = 0,
    ntyp = 0,
    nbnd = 0,
    tot_charge = 0.0,
    starting_charge = [],
    tot_magnetization = -1.0,
    starting_magnetization = [],
    ecutwfc = 0.0,
    ecutrho = 4ecutwfc,
    ecutfock = ecutrho,
    nr1 = 0,
    nr2 = 0,
    nr3 = 0,
    nr1s = 0,
    nr2s = 0,
    nr3s = 0,
    nosym = false,
    nosym_evc = false,
    noinv = false,
    no_t_rev = false,
    force_symmorphic = false,
    use_all_frac = false,
    occupations = "fixed",
    one_atom_occupations = false,
    starting_spin_angle = false,
    degauss = 0.0,
    smearing = "gaussian",
    nspin = 1,
    noncolin = false,
    ecfixed = 0.0,
    qcutz = 0.0,
    q2sigma = 0.1,  # The default value in QE's source code is 0.01
    input_dft = "none",
    exx_fraction = 0.25,
    screening_parameter = 0.106,
    exxdiv_treatment = "gygi-baldereschi",
    x_gamma_extrapolation = true,
    ecutvcut = 0.0,
    nqx1 = 1,
    nqx2 = 1,
    nqx3 = 1,
    localization_thr = 0.0,  # This is only for QE 6.4
    lda_plus_u = false,
    lda_plus_u_kind = 0,
    Hubbard_U = [],
    Hubbard_J0 = [],
    Hubbard_alpha = [],
    Hubbard_beta = [],
    # Hubbard_J = [zeros(ntyp)]  ,  # The default value in QE's source code is just one 0.0
    starting_ns_eigenvalue = -1.0,  # It's actually a multidimensional array.
    U_projection_type = "atomic",
    edir = 1,
    emaxpos = 0.5,
    eopreg = 0.1,
    eamp = 0.001,  # The default value in QE's source code is 0.0
    angle1 = [],
    angle2 = [],
    constrained_magnetization = "none",
    fixed_magnetization = zeros(3),  # The default value in QE's source code is just one 0.0
    lambda = 1.0,
    report = 100,
    lspinorb = false,
    assume_isolated = "none",
    esm_bc = "pbc",
    esm_w = 0.0,
    esm_efield = 0.0,
    esm_nfit = 4,
    fcp_mu = 0.0,
    vdw_corr = "none",
    london = false,
    london_s6 = 0.75,
    london_c6 = [],
    london_rvdw = [],
    london_rcut = 200.0,
    ts_vdw_econv_thr = 1e-06,
    ts_vdw_isolated = false,
    xdm = false,
    xdm_a1 = 0.6836,  # The default value in QE's source code is 0.0
    xdm_a2 = 1.5045,  # The default value in QE's source code is 0.0
    space_group = 0,
    uniqueb = false,
    origin_choice = 1,
    rhombohedral = true,
    zgate = 0.5,
    relaxz = false,
    block = false,
    block_1 = 0.45,
    block_2 = 0.55,
    block_height = 0.1,  # The default value in QE's source code is 0.0
)
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1378-L1499.
    @assert ibrav in union(0:1:14, (-3, -5, -9, 91, -12, -13), 127)
    @assert ntyp <= 10 && ntyp <= nat
    @assert smearing in (
        "gaussian",
        "gauss",
        "methfessel-paxton",
        "m-p",
        "mp",
        "marzari-vanderbilt",
        "cold",
        "m-v",
        "mv",
        "fermi-dirac",
        "f-d",
        "fd",
    )
    @assert nspin in (1, 2, 4)
    @assert ecutwfc >= 0  # From `set_cutoff` in https://github.com/QEF/q-e/blob/7573301/PW/src/input.f90#L1597-L1639
    @assert ecutrho >= ecutwfc # From `set_cutoff`
    @assert ecfixed >= 0
    @assert qcutz >= 0
    @assert q2sigma >= 0
    @assert lda_plus_u_kind in 0:1
    @assert edir in 1:3
    @assert space_group in 0:230
    @assert origin_choice in 1:2
    @assert length(starting_charge) <= ntyp
    @assert length(starting_magnetization) <= ntyp
    @assert length(Hubbard_U) <= ntyp
    @assert length(Hubbard_J0) <= ntyp
    @assert length(Hubbard_alpha) <= ntyp
    @assert length(Hubbard_beta) <= ntyp
    # @assert all(length(x) <= ntyp for x in Hubbard_J)
    @assert length(angle1) <= ntyp
    @assert length(angle2) <= ntyp
    @assert length(fixed_magnetization) <= 3
    @assert length(london_c6) <= ntyp
    @assert length(london_rvdw) <= ntyp
    @assert exxdiv_treatment in
            ("gygi-baldereschi", "gygi-bald", "g-b", "vcut_ws", "vcut_spherical", "none")
    @assert !(x_gamma_extrapolation && exxdiv_treatment in ("vcut_ws", "vcut_spherical")) "`x_gamma_extrapolation` cannot be used with `vcut`!"
    return SystemNamelist(
        ibrav,
        celldm,
        A,
        B,
        C,
        cosAB,
        cosAC,
        cosBC,
        nat,
        ntyp,
        nbnd,
        tot_charge,
        starting_charge,
        tot_magnetization,
        starting_magnetization,
        ecutwfc,
        ecutrho,
        ecutfock,
        nr1,
        nr2,
        nr3,
        nr1s,
        nr2s,
        nr3s,
        nosym,
        nosym_evc,
        noinv,
        no_t_rev,
        force_symmorphic,
        use_all_frac,
        occupations,
        one_atom_occupations,
        starting_spin_angle,
        degauss,
        smearing,
        nspin,
        noncolin,
        ecfixed,
        qcutz,
        q2sigma,
        input_dft,
        exx_fraction,
        screening_parameter,
        exxdiv_treatment,
        x_gamma_extrapolation,
        ecutvcut,
        nqx1,
        nqx2,
        nqx3,
        localization_thr,
        lda_plus_u,
        lda_plus_u_kind,
        Hubbard_U,
        Hubbard_J0,
        Hubbard_alpha,
        Hubbard_beta,
        starting_ns_eigenvalue,
        U_projection_type,
        edir,
        emaxpos,
        eopreg,
        eamp,
        angle1,
        angle2,
        constrained_magnetization,
        fixed_magnetization,
        lambda,
        report,
        lspinorb,
        assume_isolated,
        esm_bc,
        esm_w,
        esm_efield,
        esm_nfit,
        fcp_mu,
        vdw_corr,
        london,
        london_s6,
        london_c6,
        london_rvdw,
        london_rcut,
        ts_vdw_econv_thr,
        ts_vdw_isolated,
        xdm,
        xdm_a1,
        xdm_a2,
        space_group,
        uniqueb,
        origin_choice,
        rhombohedral,
        zgate,
        relaxz,
        block,
        block_1,
        block_2,
        block_height,
    )
end

"""
    ElectronsNamelist <: Namelist
    ElectronsNamelist(; kwargs...)

Represent the `ELECTRONS` namelist of `pw.x`.
"""
struct ElectronsNamelist <: Namelist
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
end # struct ElectronsNamelist

function ElectronsNamelist(;
    electron_maxstep = 100,
    scf_must_converge = true,
    conv_thr = 1e-6,
    adaptive_thr = false,
    conv_thr_init = 1e-3,
    conv_thr_multi = 0.1,
    mixing_mode = "plain",
    mixing_beta = 0.7,
    mixing_ndim = 8,
    mixing_fixed_ns = 0,
    diagonalization = "david",
    ortho_para = 0,
    diago_thr_init = 0.0,
    diago_cg_maxiter = 20,
    diago_david_ndim = 4,
    diago_full_acc = false,
    efield = 0.0,
    efield_cart = zeros(3),
    efield_phase = "none",
    startingpot = "atomic",  # This depends on `calculation`
    startingwfc = "atomic+random",  # This depends on `calculation`
    tqr = false,
)
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1508-L1543.
    @assert mixing_mode in ("plain", "TF", "local-TF")
    @assert diagonalization in ("david", "cg", "cg-serial", "david-serial", "ppcg")  # Different from docs
    @assert efield_phase in ("read", "write", "none")
    @assert startingpot in ("atomic", "file")
    @assert startingwfc in ("atomic", "atomic+random", "random", "file")
    return ElectronsNamelist(
        electron_maxstep,
        scf_must_converge,
        conv_thr,
        adaptive_thr,
        conv_thr_init,
        conv_thr_multi,
        mixing_mode,
        mixing_beta,
        mixing_ndim,
        mixing_fixed_ns,
        diagonalization,
        ortho_para,
        diago_thr_init,
        diago_cg_maxiter,
        diago_david_ndim,
        diago_full_acc,
        efield,
        efield_cart,
        efield_phase,
        startingpot,
        startingwfc,
        tqr,
    )
end
