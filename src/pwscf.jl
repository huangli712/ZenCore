#
# Project : Pansy
# Source  : pwscf.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/08/17
#

function pwscf_adaptor()
end




"""
    Input

An abstract type representing an input object of ab initio software.
All other input types should subtype `Input`.
"""
abstract type Input end

"""
    InputEntry

Represent any component of an `Input`. The fields of an `Input` should
all be either `InputEntry` or `Nothing` (no value provided).
"""
abstract type InputEntry end

"""
    Namelist <: InputEntry

Represent a component of an `Input`, a basic Fortran data structure.
"""
abstract type Namelist <: InputEntry end

"""
    Card <: InputEntry

Represent cards of an `Input` in Quantum ESPRESSO.
"""
abstract type Card <: InputEntry end

# The default values are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90.
"""
    ControlNamelist <: Namelist
    ControlNamelist(; kwargs...)

Represent the `CONTROL` namelist of `pw.x`.
"""
@auto_hash_equals struct ControlNamelist <: Namelist
    calculation::String
    title::String
    verbosity::String
    restart_mode::String
    wf_collect::Bool
    nstep::UInt
    iprint::UInt
    tstress::Bool
    tprnfor::Bool
    dt::Float64
    outdir::String
    wfcdir::String
    prefix::String
    lkpoint_dir::Bool
    max_seconds::Float64
    etot_conv_thr::Float64
    forc_conv_thr::Float64
    disk_io::String
    pseudo_dir::String
    tefield::Bool
    dipfield::Bool
    lelfield::Bool
    nberrycyc::UInt
    lorbm::Bool
    lberry::Bool
    gdir::UInt8
    nppstr::UInt
    lfcpopt::Bool
    gate::Bool
end # struct ControlNamelist

function ControlNamelist(;
    calculation = "scf",
    title = " ",
    verbosity = "low",
    restart_mode = "from_scratch",
    wf_collect = true,
    nstep = 50,
    iprint = 100000,
    tstress = false,
    tprnfor = false,
    dt = 20.0,
    outdir = "./",
    wfcdir = "./",
    prefix = "pwscf",
    lkpoint_dir = true,
    max_seconds = 10000000.0,
    etot_conv_thr = 1e-4,
    forc_conv_thr = 1e-3,
    disk_io = ifelse(calculation == "scf", "low", "medium"),
    pseudo_dir = raw"$HOME/espresso/pseudo/",
    tefield = false,
    dipfield = false,
    lelfield = false,
    nberrycyc = 1,
    lorbm = false,
    lberry = false,
    gdir = 1,  # The QE default value is `0`!
    nppstr = 0,
    lfcpopt = false,
    gate = false,
)
    # These checks are from https://github.com/QEF/q-e/blob/4132a64/Modules/read_namelists.f90#L1282-L1369.
    @assert calculation in ("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md")
    @assert verbosity in ("high", "low", "debug", "medium", "default", "minimal")
    @assert restart_mode in ("from_scratch", "restart")
    @assert iprint >= 1
    @assert disk_io in ("high", "medium", "low", "none", "default")
    @assert dt >= 0
    # @assert !(lkpoint_dir && wf_collect) "`lkpoint_dir` currently doesn't work together with `wf_collect`!"
    @assert max_seconds >= 0
    @assert etot_conv_thr >= 0
    @assert forc_conv_thr >= 0
    @assert gdir in 1:3
    @assert !all((gate, tefield, !dipfield)) "`gate` cannot be used with `tefield` if dipole correction is not active!"
    @assert !all((gate, dipfield, !tefield)) "dipole correction is not active if `tefield = false`!"
    return ControlNamelist(
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
        lfcpopt,
        gate,
    )
end

xmldir(nml::ControlNamelist) = expanduser(joinpath(nml.outdir, nml.prefix * ".save"))
wfcfiles(nml::ControlNamelist, n = 1) =
    [joinpath(xmldir(nml), nml.prefix * ".wfc$i") for i in 1:n]
