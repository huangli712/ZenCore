#
# Project : Pansy
# Source  : solver.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/28
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
    # Generate configuration file for quantum impurity solver
    ctqmc_setup(imp)

    # Extract frequency mesh and hybridization function from `dmft.hyb_l`
    fmesh, Delta = GetHyb_l(imp)

    # Write frequency mesh and hybridization function to `solver.hyb.in`
    ctqmc_hyb_l(fmesh, Delta)

    # Extract local impurity levels from `dmft.eimpx`
    Eimpx = GetEimpx(imp)

    # Write local impurity levels to `solver.eimp.in`
    ctqmc_eimpx(Eimpx)
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
# Service Functions: For I/O Operations
#

"""
    ctqmc_setup(imp::Impurity)

Generate configuration file (`solver.ctqmc.in`) for CT-QMC quantum
impurity solvers automatically (according to the information encoded
in the Impurity struct).

See also: [`Impurity`](@ref), [`ctqmc_atomx`](@ref).
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

        # Print the standard parameters
        println(fout, "# Standard parameters: Running mode")
        println(fout, "isscf = 1")
        println(fout)
        #
        println(fout, "# Standard parameters: System's size")
        println(fout, "nband = $nband")
        println(fout, "norbs = $norbs")
        println(fout, "ncfgs = $ncfgs")
        println(fout)
        #
        println(fout, "# Standard parameters: Interaction")
        println(fout, "Uc    = $Uc")
        println(fout, "Jz    = $Jz")
        println(fout)
        #
        println(fout, "# Standard parameters: Others")
        println(fout, "mune  = $mune")
        println(fout, "beta  = $beta")
        println(fout)

        # Print the user-supplied parameters
        #
        # The correctness of these auxiliary parameters should be checked
        # in config.jl.
        println(fout, "# Auxiliary parameters: By users")
        foreach(x -> println(fout, x), get_s("params"))

    end
end

"""
    ctqmc_atomx(imp::Impurity)

Generate configuration file for the atomic problem solver.

See also: [`Impurity`](@ref), [`ctqmc_setup`](@ref).
"""
function ctqmc_atomx(imp::Impurity)
    sorry()
end

"""
    ctqmc_hyb_l(fmesh::Array{F64,1}, Delta::Array{C64,4})

Write the hybridization functions to the `solver.hyb.in` file, which is
suitable for the CT-QMC quantum impurity solver.

See also: [`ctqmc_eimpx`](@ref).
"""
function ctqmc_hyb_l(fmesh::Array{F64,1}, Delta::Array{C64,4})
    # Extract key parameters from `Delta`
    _, nband, nmesh, nspin = size(Delta)

    # Write the hybridization functions to `solver.hyb.in`
    open("solver.hyb.in", "w") do fout
        # If nspin is 2, everything is OK. However, if nspin is 1, the
        # spin-down part is missed.
        for s = 1:nspin
            for p = 1:nband
                orb = (s - 1) * nband + p
                for m = 1:nmesh
                    z = Delta[p,p,m,s]
                    @printf(fout, "%6i%16.8f%16.8f%16.8f%16.8f%16.8f\n",
                        orb, fmesh[m], real(z), imag(z), 0.0, 0.0)
                end
                println(fout)
                println(fout)
            end
        end # END OF S LOOP

        # Supplement the spin-down part if it is missed. Here, we just
        # assume that the spin-up and spin-down parts are degenerated.
        if nspin == 1
            for s = 1:nspin
                for p = 1:nband
                    orb = (2 - 1) * nband + p
                    for m = 1:nmesh
                        z = Delta[p,p,m,s]
                        @printf(fout, "%6i%16.8f%16.8f%16.8f%16.8f%16.8f\n",
                            orb, fmesh[m], real(z), imag(z), 0.0, 0.0)
                    end
                    println(fout)
                    println(fout)
                end
            end # END OF S LOOP
        end
    end
end

"""
    ctqmc_eimpx(Eimpx::Array{C64,3})

Write the local impurity levels to the `solver.eimp.in` file, which is
suitable for the CT-QMC quantum impurity solver.

See also: [`ctqmc_hyb_l`](@ref).
"""
function ctqmc_eimpx(Eimpx::Array{C64,3})
    # Extract key parameters from `Eimpx`
    _, nband, nspin = size(Eimpx)

    # Analyze the orbital degeneracy
    symm = GetSymmetry(Eimpx)
    println(symm)

    # Write the local impurity levels to `solver.eimp.in`
    open("solver.eimp.in", "w") do fout
        # If nspin is 2, everything is OK. However, if nspin is 1, the
        # spin-down part is missed.
        for s = 1:nspin
            for p = 1:nband
                orb = (s - 1) * nband + p
                z = Eimpx[p, p, s]
                @printf(fout, "%4i%16.8f%4i\n", orb, real(z), orb)
            end
        end # END OF S LOOP

        # Supplement the spin-down part if it is missed. Here, we just
        # assume that the spin-up and spin-down parts are degenerated.
        if nspin == 1
            for s = 1:nspin
                for p = 1:nband
                    orb = (2 - 1) * nband + p
                    z = Eimpx[p, p, s]
                    @printf(fout, "%4i%16.8f%4i\n", orb, real(z), orb)
                end
            end # END OF S LOOP
        end
    end
end

"""
    ctqmc_sig_l(imp::Impurity)

Parse the `solver.sgm.dat` file to extract the bare self-energy functions.

In the sigma_gather() function, these data will be combined to generate
the `sigma.bare` file, which is essential for the DMFT engine.

See also: [`Impurity`](@ref), [`GetSig_l`](@ref).
"""
function ctqmc_sig_l(imp::Impurity)
    # File name for self-energy functions
    fsgm = "solver.sgm.dat"

    # To make sure the data file is present
    @assert isfile(fsgm)

    # Extract parameters from Impurity struct
    nband = imp.nband
    norbs = nband * 2
    nspin = 2 # In the CT-QMC impurity solver, nspin is fixed to 2.

    # We don't know how many frequency points are used a priori.
    # Here, we used a trick to determine `nmesh`.
    lines = readline(fsgm)
    nmesh = length(lines) / norbs - 2

    # Create an array for frequency mesh
    fmesh = zeros(F64, nmesh)

    # Create an array for self-energy functions
    Σ = zeros(C64, nband, nband, nmesh, nspin)

    # Parse the data file for the self-energy functions
    open(fsgm, "r") do fin

        # Go through each spin orientation
        for s = 1:nspin

            # Go through each band
            for b = 1:nband

                # Go through each frequency point
                for m = 1:nmesh

                    # Split the line
                    arr = line_to_array(fin)

                    # Extract frequency
                    fmesh[m] = parse(F64, arr[2])

                    # Extract self-energy functions
                    _re = parse(F64, arr[3])
                    _im = parse(F64, arr[4])
                    Σ[b,b,m,s] = _re + _im * im

                end # END OF M LOOP

                # Skip two lines
                readline(fin)
                readline(fin)

            end # END OF B LOOP

        end # END OF S LOOP

    end

    # Return the desired arrays
    return (fmesh, Σ)
end

"""
    ctqmc_nimpx(imp::Impurity)

Parse the `solver.nmat.dat` file to extract the impurity occupancy. Then
the field `occup` in Impurity struct will be updated.

In this function, only the total impurity occupancy of the current site
is return. However, sometimes we need to known the spin-up and spin-down
components. Later, we will expand the Impurity struct and this function
to fulfill this requirement.

See also: [`Impurity`](@ref), [`GetNimpx`](@ref).
"""
function ctqmc_nimpx(imp::Impurity)
    # File name for impurity occupancy
    fnmat = "solver.nmat.dat"

    # To make sure the data file is present
    @assert isfile(fnmat)

    # Parse the data file to extract total impurity occupancy
    lines = readlines(fnmat)
    filter!(x -> contains(x, "sum"), lines)
    @assert length(lines) == 1
    arr = line_to_array(iters[end])
    occup = parse(F64, arr[2])

    # Update Impurity struct
    imp.occup = occup
end

#
# Service Functions: For I/O Operations
#

"""
    GetHyb_l(imp::Impurity)

Extract hybridization functions from `dmft.hyb_l` file, which is generated
by sigma_split(). The data are essential for quantum impurity solvers.

The frequency mesh is also extracted in this function.

See also: [`Impurity`](@ref), [`GetEimpx`](@ref).
"""
function GetHyb_l(imp::Impurity)
    # Get index of the quantum impurity problem
    index = imp.index

    # Get number of orbitals for the quantum impurity problem
    nband = imp.nband

    # Examine the current directory
    dirname = basename(pwd())
    dirvect = split(dirname, ".")
    @assert dirvect[1] == "impurity"
    @assert parse(I64, dirvect[2]) == index

    # Declare empty arrays for frequency mesh and hybridization functions
    fmesh = []
    Delta = []

    # Parse the `dmft.hyb_l` file
    open("dmft.hyb_l", "r") do fin

        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        nmesh = parse(I64, line_to_array(fin)[3])
        qdim  = parse(I64, line_to_array(fin)[4])
        @assert qdim ≥ nband

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for frequency mesh
        fmesh = zeros(F64, nmesh)

        # Create an array for hybridization functions
        Delta = zeros(C64, nband, nband, nmesh, nspin)

        # Go through each spin orientation
        for s = 1:nspin

            # Analyze the important parameters
            strs = readline(fin)
            _t = parse(I64, line_to_array(strs)[3])
            _s = parse(I64, line_to_array(strs)[5])
            _d = parse(I64, line_to_array(strs)[7])
            @assert _t == index
            @assert _s == s
            @assert _d == nband

            # Go through each frequency point
            for m = 1:nmesh

                # Extract frequency
                fmesh[m] = parse(F64, line_to_array(fin)[3])

                # Parse the hybridization functions
                for q = 1:nband
                    for p = 1:nband
                        _re, _im = parse.(F64, line_to_array(fin)[3:4])
                        Delta[p,q,m,s] = _re + _im * im
                    end
                end

            end # END OF M LOOP

            # Skip two lines
            readline(fin)
            readline(fin)

        end # END OF S LOOP

    end

    # Return the desired arrays
    return (fmesh, Delta)
end

"""
    GetEimpx(imp::Impurity)

Extract local impurity levels from `dmft.eimpx` file, which is generated
by sigma_split(). The data are essential for quantum impurity solvers.

See also: [`Impurity`](@ref), [`GetHyb_l`](@ref).
"""
function GetEimpx(imp::Impurity)
    # Get index of the quantum impurity problem
    index = imp.index

    # Get number of orbitals for the quantum impurity problem
    nband = imp.nband

    # Examine the current directory
    dirname = basename(pwd())
    dirvect = split(dirname, ".")
    @assert dirvect[1] == "impurity"
    @assert parse(I64, dirvect[2]) == index

    # Declare an empty array for local impurity levels
    Eimpx = []

    # Parse the `dmft.eimpx` file
    open("dmft.eimpx", "r") do fin

        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        qdim  = parse(I64, line_to_array(fin)[4])
        @assert qdim ≥ nband

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for local impurity levels
        Eimpx = zeros(C64, nband, nband, nspin)

        # Go through each spin orientation
        for s = 1:nspin

            # Analyze the important parameters
            strs = readline(fin)
            _t = parse(I64, line_to_array(strs)[3])
            _s = parse(I64, line_to_array(strs)[5])
            _d = parse(I64, line_to_array(strs)[7])
            @assert _t == index
            @assert _s == s
            @assert _d == nband

            # Parse local impurity levels
            for q = 1:nband
                for p = 1:nband
                    _re, _im = parse.(F64, line_to_array(fin)[3:4])
                    Eimpx[p,q,s] = _re + _im * im
                end
            end

            # Skip two lines
            readline(fin)
            readline(fin)
        end

    end

    # Return the desired array
    return Eimpx
end

"""
    GetSig_l(imp::Impurity)

Extract self-energy functions from the output files of various quantum
impurity solvers. The data will be combined in sigma_gather() function.
Then they will be fed back to the DMFT engine.

See also: [`Impurity`](@ref), [`ctqmc_sig_l`](@ref).
"""
function GetSig_l(imp::Impurity)
    # Get the index for current quantum impurity problem
    index = imp.index

    # Change the directory
    cd("impurity.$index")

    # Determine the chosen solver
    engine = get_s("engine")

    # Activate the corresponding solver_sig_l() functions for various
    # quantum impurity solvers
    sig_l = nothing
    @cswitch engine begin
        @case "ct_hyb1"
            sig_l = ctqmc_sig_l(imp)
            break

        @case "ct_hyb2"
            sig_l = ctqmc_sig_l(imp)
            break

        @case "hub1"
            sorry()
            break

        @case "norg"
            sorry()
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")

    # Return the desired array
    return sig_l
end

"""
    GetNimpx(imp::Impurity)

Extract impurity occupancy from the output files of various quantum
impurity solvers. Then the field `occup` in Impurity struct will be
updated, which will then be used to update the double counting term
for self-energy functions.

The argument `imp` will be modified in this function.

See also: [`Impurity`](@ref), [`ctqmc_nimpx`](@ref).
"""
function GetNimpx(imp::Impurity)
    # Get the index for current quantum impurity problem
    index = imp.index

    # Change the directory
    cd("impurity.$index")

    # Determine the chosen solver
    engine = get_s("engine")

    # Activate the corresponding solver_nimpx() functions for various
    # quantum impurity solvers
    @cswitch engine begin
        @case "ct_hyb1"
            ctqmc_nimpx(imp)
            break

        @case "ct_hyb2"
            ctqmc_nimpx(imp)
            break

        @case "hub1"
            sorry()
            break

        @case "norg"
            sorry()
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")
end

#
# Service Functions: For Impurity Problems
#

"""
    GetSymmetry(Eimpx::Array{C64,3})

Analyze the symmetry according to the diagonal elements of the matrix of
the local impurity levels.

See also: [`GetEimpx`](@ref).
"""
function GetSymmetry(Eimpx::Array{C64,3})
    # Extract some key parameters from Eimpx
    _, nband, nspin = size(Eimpx)

    # Extract the diagonal elements of Eimpx
    eimp = zeros(F64, nband, nspin)
    for s = 1:nspin
        for b = 1:nband
            eimp[b,s] = round(real(Eimpx[b,b,s]), digits = 4)
        end
    end

    # Create an array for symmetry
    symm = zeros(I64, nband, nspin)

    # Go through each spin orientation
    for s = 1:nspin

        # Get the unique energy levels
        E = unique(eimp[:,s])

        # Scan all the impurity levels, find out its corresponding value in E
        for b = 1:nband
            symm[b,s] = findfirst(isequal(eimp[b]), E)
        end

    end # END OF S LOOP

    # Return the desired array
    return symm
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

    # Go through each quantum impurity problem
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
