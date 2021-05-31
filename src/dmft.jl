#
# Project : Pansy
# Source  : dmft.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/31
#

"""
    dmft_init(it::IterInfo, task::I64)

Initialize the dynamical mean-field theory engine. Prepare the necessary
files, and generate the configuration file.

See also: [`dmft_exec`](@ref), [`dmft_save`](@ref).
"""
function dmft_init(it::IterInfo, task::I64)
    # Check the task
    @assert task in (1, 2)

    # Well, determine which files are necessary. They are defined in
    # `fsig`, `fir1`, `fir2`, and `fdmft`.
    #
    # Self-energy functions and their double counting terms
    fsig = ["sigma.bare", "sigma.dc"]
    #
    # Parameter sets within the IR format
    fir1 = ["params.ir", "maps.ir", "groups.ir", "windows.ir"]
    #
    # Kohn-Sham data within the IR format
    fir2 = ["lattice.ir", "kmesh.ir", "eigen.ir", "projs.ir"]
    #
    # Tetrahedron data are available
    if get_d("smear") === "tetra"
        push!(fir2, "tetra.ir")
    end
    #
    # Main control file for the DMFT engine
    fdmft = ("dmft.in")

    # Next, we have to copy Kohn-Sham data from `dft` to `dmft1`.
    foreach( x ->
        begin
            file_src = joinpath("../dft", x)
            file_dst = x
            cp(file_src, file_dst, force = true)
        end,
    union(fir1, fir2) )

    # Extract key parameters
    axis = get_m("axis")
    beta = get_m("beta")
    mc = ( get_m("mc") isa Missing ? 0.001 : get_m("mc") )
    lfermi = ( get_m("lfermi") isa Missing ? true : get_m("lfermi") )
    ltetra = ( get_d("smear") === "tetra" )

    # Generate essential input files, such as dmft.in, dynamically.
    # If the `dmft.in` file exists already, it will be overwritten.
    open("dmft.in", "w") do fout
        println(fout, "task = $task")
        println(fout, "axis = $axis")
        println(fout, "beta = $beta")
        println(fout, "mc   = $mc")
        println(fout, "lfermi = $lfermi")
        println(fout, "ltetra = $ltetra")
    end

    # Check essential input files
    flist = (fdmft, fsig..., fir1..., fir2...)
    for i in eachindex(flist)
        filename = flist[i]
        if !isfile(filename)
            error("Please make sure the file $filename is available")
        end
    end
end

"""
    dmft_exec(it::IterInfo, task::I64)

Execute the dynamical mean-field theory engine.

See also: [`dmft_init`](@ref), [`dmft_save`](@ref).
"""
function dmft_exec(it::IterInfo, task::I64)
    # Check the task
    @assert task in (1, 2)

    # Print the header
    println("Engine : DMFT$(subscript(task))")

    # Get the home directory of DMFT engine
    dmft_home = query_dmft()

    # Determine mpi prefix (whether the dmft is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "dmft", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  Para : Using $numproc processors")

    # Select suitable dmft program
    dmft_exe = "$dmft_home/dmft"
    @assert isfile(dmft_exe)
    println("  Exec : $dmft_exe")

    # Assemble command
    if isnothing(mpi_prefix)
        dmft_cmd = dmft_exe
    else
        dmft_cmd = split("$mpi_prefix $dmft_exe", " ")
    end

    # Launch it, the terminal output is redirected to dmft.out
    run(pipeline(`$dmft_cmd`, stdout = "dmft.out"))

    # Print the footer for a better visualization
    println()
end

"""
    dmft_save(it::IterInfo, task::I64)

Backup the files outputed by the dynamical mean-field theory engine.

See also: [`dmft_init`](@ref), [`dmft_exec`](@ref).
"""
function dmft_save(it::IterInfo, task::I64)
    # Check the task
    @assert task in (1, 2)

    # Create a list of files that need to be backup
    fdmf1 = ["dmft.out"]
    fdmf2 = ["dmft.eimps", "dmft.eimpx", "dmft.fermi"]
    fdmf3 = ["dmft.grn_l", "dmft.hyb_l", "dmft.wss_l"]

    # Be careful, the final file list depends on the task
    if task == 1
        file_list = union(fdmf1, fdmf2, fdmf3)
    else
        file_list = fdmf1
    end

    # Store the data files
    for i in eachindex(file_list)
        f = file_list[i]
        cp(f, "$f.$(it.dmft_cycle).$(it.dmft1_iter)", force = true)
    end
end

"""
    read_hyb_l(imp::Impurity)

Extract hybridization functions from `dmft.hyb_l` file, which is generated
by sigma_split(). The data are essential for quantum impurity solvers.

The frequency mesh is also extracted in this function.

See also: [`Impurity`](@ref), [`read_eimpx`](@ref).
"""
function read_hyb_l(imp::Impurity)
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

function read_hyb_l()
end

"""
    read_eimpx(imp::Impurity)

Extract local impurity levels from `dmft.eimpx` file, which is generated
by sigma_split(). The data are essential for quantum impurity solvers.

See also: [`Impurity`](@ref), [`read_hyb_l`](@ref).
"""
function read_eimpx(imp::Impurity)
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

function read_eimpx()
end
