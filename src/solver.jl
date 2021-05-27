#
# Project : Pansy
# Source  : solver.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/27
#

#
# CT-HYB₁ Quantum Impurity Solver
#

"""
    s_qmc1_init(it::IterInfo, imp::Impurity)

Check runtime environment of the CT-HYB₁ quantum impurity solver. Prepare
the necessary input files.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_exec`](@ref), [`s_qmc1_save`](@ref).
"""
function s_qmc1_init(it::IterInfo, imp::Impurity)
    ctqmc_setup(imp)
    ctqmc_hyb_l()
    ctqmc_eimpx()
end

"""
    s_qmc1_exec(it::IterInfo)

Launch the CT-HYB₁ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_init`](@ref), [`s_qmc1_save`](@ref).
"""
function s_qmc1_exec(it::IterInfo)
    # Print the header
    println("Engine : CT-HYB₁")

    # Get the home directory of quantum impurity solver
    solver_home = query_solver("ct_hyb1")

    # Determine mpi prefix (whether the solver is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "solver", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  Para : Using $numproc processors")

    # Select suitable solver program
    solver_exe = "$solver_home/ctqmc"
    @assert isfile(solver_exe)
    println("  Exec : $solver_exe")

    # Assemble command
    if isnothing(mpi_prefix)
        solver_cmd = solver_exe
    else
        solver_cmd = split("$mpi_prefix $solver_exe", " ")
    end

    # Launch it, the terminal output is redirected to solver.out
    run(pipeline(`$solver_cmd`, stdout = "solver.out"))

    # Print the footer for a better visualization
    println()
end

"""
    s_qmc1_save(it::IterInfo)

Backup output files of the CT-HYB₁ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_init`](@ref), [`s_qmc1_exec`](@ref).
"""
function s_qmc1_save(it::IterInfo)
end

#
# CT-HYB₂ Quantum Impurity Solver
#

"""
    s_qmc2_init(it::IterInfo)

Check runtime environment of the CT-HYB₂ quantum impurity solver. Prepare
the necessary input files.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_exec`](@ref), [`s_qmc2_save`](@ref).
"""
function s_qmc2_init(it::IterInfo)
    sorry()
end

"""
    s_qmc2_exec(it::IterInfo)

Launch the CT-HYB₂ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_save`](@ref).
"""
function s_qmc2_exec(it::IterInfo)
    # Print the header
    println("Engine : CT-HYB₂")

    # Print the footer for a better visualization
    println()
end

"""
    s_qmc2_save(it::IterInfo)

Backup output files of the CT-HYB₂ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_exec`](@ref).
"""
function s_qmc2_save(it::IterInfo)
    sorry()
end

#
# HUB-I Quantum Impurity Solver
#

"""
    s_hub1_init(it::IterInfo)

Check runtime environment of the HIA quantum impurity solver. Prepare
the necessary input files.

See also: [`s_hub1_exec`](@ref), [`s_hub1_save`](@ref).
"""
function s_hub1_init(it::IterInfo)
    sorry()
end

"""
    s_hub1_exec(it::IterInfo)

Launch the HIA quantum impurity solver.

See also: [`s_hub1_init`](@ref), [`s_hub1_save`](@ref).
"""
function s_hub1_exec(it::IterInfo)
    # Print the header
    println("Engine : HIA")

    # Print the footer for a better visualization
    println()
end

"""
    s_hub1_save(it::IterInfo)

Backup output files of the HIA quantum impurity solver.

See also: [`s_hub1_init`](@ref), [`s_hub1_exec`](@ref).
"""
function s_hub1_save(it::IterInfo)
    sorry()
end

#
# NORG Quantum Impurity Solver
#

"""
    s_norg_init(it::IterInfo)

Check runtime environment of the NORG quantum impurity solver. Prepare
the necessary input files.

See also: [`s_norg_exec`](@ref), [`s_norg_save`](@ref).
"""
function s_norg_init(it::IterInfo)
    sorry()
end

"""
    s_norg_exec(it::IterInfo)

Launch the NORG quantum impurity solver.

See also: [`s_norg_init`](@ref), [`s_norg_save`](@ref).
"""
function s_norg_exec(it::IterInfo)
    # Print the header
    println("Engine : NORG")

    # Print the footer for a better visualization
    println()
end

"""
    s_norg_save(it::IterInfo)

Backup output files of the NORG quantum impurity solver.

See also: [`s_norg_init`](@ref), [`s_norg_exec`](@ref).
"""
function s_norg_save(it::IterInfo)
    sorry()
end

#
# Service Functions
#

"""
    ctqmc_setup(imp::Impurity)

Generate configuration file (`solver.ctqmc.in`) for CT-QMC quantum
impurity solvers automatically (according to the information encoded
in the Impurity struct).

See also: [`Impurity`](@ref).
"""
function ctqmc_setup(imp::Impurity)
    # File name of configuration file
    fctqmc = "solver.ctqmc.in"

    # Extract parameters from `Impurity` struct 
    #
    # mune is fixed to zero, because the chemical potential is adsorbed
    # in the local impurity level.
    nband = imp.nband
    norbs = 2 * nband
    ncfgs = 2^norbs
    Uc    = imp.upara
    Jz    = imp.jpara
    mune  = 0.0
    beta  = imp.beta

    # Create the configuration file
    open(fctqmc, "w") do fout
        # Print some warning message 
        println(fout, "# Generated by $__LIBNAME__ at ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
        println(fout, "# PLEASE DO NOT MODIFY IT MANUALLY!")
        println(fout)

        println(fout, "# Standard  parameters: Running mode")
        println(fout, "isscf = 1")
        println(fout)

        println(fout, "# Standard  parameters: System's size")
        println(fout, "nband = imp.nband")
        println(fout, "norbs = 10")
        println(fout, "ncfgs = 1024")
        println(fout)

        println(fout, "# Standard  parameters: Interaction")
        println(fout, "Uc    = 4.0")
        println(fout, "Jz    = 0.7")
        println(fout)

        println(fout, "# Standard  parameters: Others")
        println(fout, "mune  = 0.0")
        println(fout, "beta  = 40.0")
        println(fout)

    end
end

function ctqmc_hyb_l()
    fmesh = []
    Delta = []
    cdim  = 0

    open("dmft.hyb_l", "r") do fin

        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        nmesh = parse(I64, line_to_array(fin)[3])
        qdim = parse(I64, line_to_array(fin)[4])

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for frequency mesh
        fmesh = zeros(F64, nmesh)

        # Create an array for hybridization functions
        Delta = zeros(C64, qdim, qdim, nmesh, nspin)

        for s = 1:nspin
            strs = readline(fin)
            _t = parse(I64, line_to_array(strs)[3])
            _s = parse(I64, line_to_array(strs)[5])
            cdim = parse(I64, line_to_array(strs)[7])
            @assert _t == 1 && _s == s
            for m = 1:nmesh
                fmesh[m] = parse(F64, line_to_array(fin)[3])
                # Parse hybridization functions
                for q = 1:cdim
                    for p = 1:cdim
                        _re, _im = parse.(F64, line_to_array(fin)[3:4])
                        Delta[p,q,m,s] = _re + _im * im
                    end
                end
            end
            # Skip two lines
            readline(fin)
            readline(fin)
        end
    end

    _, nband, nmesh, nspin = size(Delta)
    @assert nband >= cdim

    open("solver.hyb.in", "w") do fout
        for s = 1:nspin
            for p = 1:cdim
                orb = (s - 1) * cdim + p
                for m = 1:nmesh
                    z = Delta[p,p,m,s]
                    @printf(fout, "%6i%16.8f%16.8f%16.8f%16.8f%16.8f\n", orb, fmesh[m], real(z), imag(z), 0.0, 0.0)
                end
                println(fout)
                println(fout)
            end
        end

        if nspin == 1
        for s = 1:nspin
            for p = 1:cdim
                orb = (2 - 1) * cdim + p
                for m = 1:nmesh
                    z = Delta[p,p,m,s]
                    @printf(fout, "%6i%16.8f%16.8f%16.8f%16.8f%16.8f\n", orb, fmesh[m], real(z), imag(z), 0.0, 0.0)
                end
                println(fout)
                println(fout)
            end
        end
        end
    end
end

function ctqmc_eimpx()
    Eimpx = []
    cdim  = 0

    open("dmft.eimpx", "r") do fin
    
        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        qdim = parse(I64, line_to_array(fin)[4])

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for local impurity levels
        Eimpx = zeros(C64, qdim, qdim, nspin)

        for s = 1:nspin
            strs = readline(fin)
            _t = parse(I64, line_to_array(strs)[3])
            _s = parse(I64, line_to_array(strs)[5])
            cdim = parse(I64, line_to_array(strs)[7])
            @assert _t == 1 && _s == s

            # Parse local impurity levels
            for q = 1:cdim
                for p = 1:cdim
                    _re, _im = parse.(F64, line_to_array(fin)[3:4])
                    Eimpx[p,q,s] = _re + _im * im
                end
            end

            # Skip two lines
            readline(fin)
            readline(fin)
        end

    end

    _, nband, nspin = size(Eimpx)
    @assert nband >= cdim

    open("solver.eimp.in", "w") do fout
        for s = 1:nspin
            for p = 1:cdim
                orb = (s - 1) * cdim + p
                z = Eimpx[p, p, s]
                @printf(fout, "%4i%16.8f%4i\n", orb, real(z), orb)
            end
        end

        if nspin == 1
        for s = 1:nspin
            for p = 1:cdim
                orb = (2 - 1) * cdim + p
                z = Eimpx[p, p, s]
                @printf(fout, "%4i%16.8f%4i\n", orb, real(z), orb)
            end
        end
        end
    end
end

function ctqmc_sig_l()
end

function ctqmc_nimpx()
end

"""
    GetImpurity()

Return an array of Impurity struct, which encapsulates useful information
about the quantum impurity problems.
  
See also: [`Impurity`](@ref).
"""
function GetImpurity()
    # Initialize an array of Impurity struct
    AI = Impurity[]

    # Go through each quantum impurity problems
    for i = 1:get_i("nsite")

        # Extract essential parameters
        str = get_i("atoms")[i]
        index = i
        atoms = String(line_to_array(str)[1])
        sites = parse(I64, line_to_array(str)[3])
        shell = get_i("shell")[i]
        ising = get_i("ising")[i]
        occup = get_i("occup")[i]
        upara = get_i("upara")[i]
        jpara = get_i("jpara")[i]
        lpara = get_i("lpara")[i]
        beta  = get_m("beta")

        # Call the constructor
        Im = Impurity(index, atoms, sites, shell, ising, occup, upara, jpara, lpara, beta)

        # Save Im in AI
        push!(AI, Im)

    end

    # Return the desired array
    return AI
end
