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

Represent the `CONTROL` namelist of `pw.x`.
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
    ControlNL(; kwargs...)
"""
function ControlNL(;
    calculation   = "scf",
    title         = " ",
    verbosity     = "low",
    restart_mode  = "from_scratch",
    wf_collect    = true,
    nstep         = 50,
    iprint        = 100000,
    tstress       = false,
    tprnfor       = false,
    dt            = 20.0,
    outdir        = "./",
    wfcdir        = "./",
    prefix        = "pwscf",
    lkpoint_dir   = true,
    max_seconds   = 10000000.0,
    etot_conv_thr = 1e-4,
    forc_conv_thr = 1e-3,
    disk_io       = ifelse(calculation == "scf", "low", "medium"),
    pseudo_dir    = raw"$HOME/espresso/pseudo/",
    tefield       = false,
    dipfield      = false,
    lelfield      = false,
    nberrycyc     = 1,
    lorbm         = false,
    lberry        = false,
    gdir          = 1,
    nppstr        = 0,
    lfcp          = false,
    gate          = false
)

    # Sanity check
    @assert calculation in ("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md")
    @assert verbosity in ("high", "low", "debug", "medium", "default", "minimal")
    @assert restart_mode in ("from_scratch", "restart")
    @assert iprint >= 1
    @assert disk_io in ("high", "medium", "low", "none", "default")
    @assert dt >= 0
    @assert max_seconds >= 0
    @assert etot_conv_thr >= 0
    @assert forc_conv_thr >= 0
    @assert gdir in 1:3
    @assert !all((gate, tefield, !dipfield))
    @assert !all((gate, dipfield, !tefield))

    # Call the default constructor
    return ControlNL(
        calculation,
        title,
        verbosity,
        restart_mode,
        wf_collect,
        nstep,
        iprint,
        tstress,
        tprnfor,
        dt,
        outdir,
        wfcdir,
        prefix,
        lkpoint_dir,
        max_seconds,
        etot_conv_thr,
        forc_conv_thr,
        disk_io,
        pseudo_dir,
        tefield,
        dipfield,
        lelfield,
        nberrycyc,
        lorbm,
        lberry,
        gdir,
        nppstr,
        lfcp,
        gate
    )
end

export xmldir, wfcfiles
xmldir(nml::ControlNL) = expanduser(joinpath(nml.outdir, nml.prefix * ".save"))
wfcfiles(nml::ControlNL, n = 1) =
    [joinpath(xmldir(nml), nml.prefix * ".wfc$i") for i in 1:n]

export SystemNL

"""
    SystemNL <: Namelist

Represent the `SYSTEM` namelist of `pw.x`.
"""
struct SystemNL <: Namelist
    ibrav                     :: Int8
    celldm                    :: Vector{Maybe{Float64}}
    A                         :: Float64
    B                         :: Float64
    C                         :: Float64
    cosAB                     :: Float64
    cosAC                     :: Float64
    cosBC                     :: Float64
    nat                       :: UInt
    ntyp                      :: UInt8
    nbnd                      :: UInt
    tot_charge                :: Float64
    starting_charge           :: Vector{Maybe{Float64}}
    tot_magnetization         :: Float64
    starting_magnetization    :: Vector{Maybe{Float64}}
    ecutwfc                   :: Float64
    ecutrho                   :: Float64
    ecutfock                  :: Float64
    nr1                       :: UInt
    nr2                       :: UInt
    nr3                       :: UInt
    nr1s                      :: UInt
    nr2s                      :: UInt
    nr3s                      :: UInt
    nosym                     :: Bool
    nosym_evc                 :: Bool
    noinv                     :: Bool
    no_t_rev                  :: Bool
    force_symmorphic          :: Bool
    use_all_frac              :: Bool
    occupations               :: String
    one_atom_occupations      :: Bool
    starting_spin_angle       :: Bool
    degauss                   :: Float64
    smearing                  :: String
    nspin                     :: UInt8
    noncolin                  :: Bool
    ecfixed                   :: Float64
    qcutz                     :: Float64
    q2sigma                   :: Float64
    input_dft                 :: String
    ace                       :: Bool
    exx_fraction              :: Float64
    screening_parameter       :: Float64
    exxdiv_treatment          :: String
    x_gamma_extrapolation     :: Bool
    ecutvcut                  :: Float64
    nqx1                      :: UInt
    nqx2                      :: UInt
    nqx3                      :: UInt
    localization_thr          :: Float64
    lda_plus_u                :: Bool
    lda_plus_u_kind           :: UInt8
    Hubbard_U                 :: Vector{Maybe{Float64}}
    Hubbard_J0                :: Vector{Maybe{Float64}}
    Hubbard_alpha             :: Vector{Maybe{Float64}}
    Hubbard_beta              :: Vector{Maybe{Float64}}
    Hubbard_J                 :: Matrix{Maybe{Float64}}
    starting_ns_eigenvalue    :: Float64
    U_projection_type         :: String
    Hubbard_parameters        :: String
    ensemble_energies         :: Bool
    edir                      :: UInt8
    emaxpos                   :: Float64
    eopreg                    :: Float64
    eamp                      :: Float64
    angle1                    :: Vector{Maybe{Float64}}
    angle2                    :: Vector{Maybe{Float64}}
    lforcet                   :: Bool
    constrained_magnetization :: String
    fixed_magnetization       :: Vector{Maybe{Float64}}
    lambda                    :: Float64
    report                    :: UInt
    lspinorb                  :: Bool
    assume_isolated           :: String
    esm_bc                    :: String
    esm_w                     :: Float64
    esm_efield                :: Float64
    esm_nfit                  :: UInt
    lgcscf                    :: Bool
    gcscf_mu                  :: Float64
    gcscf_conv_thr            :: Float64
    gcscf_beta                :: Float64
    vdw_corr                  :: String
    london                    :: Bool
    london_s6                 :: Float64
    london_c6                 :: Vector{Maybe{Float64}}
    london_rvdw               :: Vector{Maybe{Float64}}
    london_rcut               :: Float64
    dftd3_version             :: UInt8
    dftd3_threebody           :: Bool
    ts_vdw_econv_thr          :: Float64
    ts_vdw_isolated           :: Bool
    xdm                       :: Bool
    xdm_a1                    :: Float64
    xdm_a2                    :: Float64
    space_group               :: UInt8
    uniqueb                   :: Bool
    origin_choice             :: UInt8
    rhombohedral              :: Bool
    zgate                     :: Float64
    relaxz                    :: Bool
    block                     :: Bool
    block_1                   :: Float64
    block_2                   :: Float64
    block_height              :: Float64
end

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
