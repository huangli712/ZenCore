#
# Project : Pansy
# Source  : wannier.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/09/10
#

#=
### *Driver Functions*
=#

"""
    wannier_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})

Adaptor support. It will firstly launch wannier90 + pw2wannier90 codes to
generate maximally localized wannier functions and related transformation
matrix. Then it will read and parse the outputs, convert the data into
IR format. The data contained in `D` dict will be modified.

Be careful, now this adaptor only supports `pwscf`.

See also: [`pwscf_adaptor`](@ref), [`ir_adaptor`](@ref).
"""
function wannier_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})
    # Print the header
    println("Adaptor : WANNIER")
    println("Try to process the Kohn-Sham dataset")
    println("Current directory: ", pwd())

    # Extract key parameters
    case = get_c("case") # Prefix for pwscf
    sp = get_d("lspins") # Is it a spin-polarized system

    # Now this feature require pwscf as a dft engine
    @assert get_d("engine") == "pwscf"

#=
    # W01: Execute wannier90 to generate w90.nnkp.
    if sp # For spin-polarized system
        # Spin up
        wannier_inir(D, "up")
        wannier_exec("up", op = "-pp")
        wannier_save("up", op = "-pp")
        # Spin down
        wannier_init(D, "dn")
        wannier_exec("dn", op = "-pp")
        wannier_save("dn", op = "-pp")
    else # For spin-unpolarized system
        wannier_init(D)
        wannier_exec(op = "-pp")
        wannier_save(op = "-pp")
    end

    # W02: Execute pw2wannier90 to generate w90.amn, w90.mmn, etc.
    if sp # For spin-polarized system
        # Spin up
        pw2wan_init(case, "up")
        pw2wan_exec(case, "up")
        pw2wan_save("up")
        # Spin down
        pw2wan_init(case, "dn")
        pw2wan_exec(case, "dn")
        pw2wan_save("dn")
    else # For spin-unpolarized system
        pw2wan_init(case)
        pw2wan_exec(case)
        pw2wan_save()
    end

    # W03: Execute wannier90 again to generate wannier functions.
    if sp # For spin-polarized system
        # Spin up
        wannier_init(D, "up")
        wannier_exec("up")
        wannier_save("up")
        # Spin down
        wannier_init(D, "dn")
        wannier_exec("dn")
        wannier_save("dn")
    else # For spin-unpolarized system
        wannier_init(D)
        wannier_exec()
        wannier_save()
    end
=#

    # Read accurate eigenvalues from w90.eig
    if sp # For spin-polarized system
        # Get eigenvalues for spin up and down components
        eigs_up = w90_read_eigs("up")
        eigs_dn = w90_read_eigs("dn")
        # Sanity check
        @assert size(eigs_up) == size(eigs_dn)
        # Convert eigs_up and eigs_dn to 3D array (nband, nkpt, 1)
        nband, nkpt = size(eigs_up)
        eigs_up = reshape(eigs_up, (nband, nkpt, 1))
        eigs_dn = reshape(eigs_dn, (nband, nkpt, 1))
        # Concatenate eigs_up and eigs_dn
        D[:enk] = cat(eigs_up, eigs_dn, dims = 3)
        # Sanity check
        @assert size(D[:enk]) == (nband, nkpt, 2)
    else # For spin-unpolarized system
        eigs = w90_read_eigs()
        nband, nkpt = size(eigs)
        D[:enk] = reshape(eigs, (nband, nkpt, 1))
    end

    latt =D[:latt]
    if sp
        D[:PT], D[:PG] = w90_make_group(latt, "up")
        PT, PG = w90_make_group(latt, "dn")
        append!(D[:PT], PT)
        append!(D[:PG], PG)
    else
        D[:PT], D[:PG] = w90_make_group(latt)
        @show typeof(D[:PT]), typeof(D[:PG])
    end

    sorry()
end

#=
### *Service Functions* : *Group A*
=#

"""
    wannier_init(D::Dict{Symbol,Any}, sp::String = "")

Try to generate the `w90.win` file, which is the essential input for
the `wannier90` code. Here, we always use `w90` as the seedname. If
the system is spin polarized (the argument `sp` is `up` or `dn`), then
the seednames will be `w90up` and `w90dn`, respectively.

See also: [`wannier_exec`](@ref), [`wannier_save`](@ref).
"""
function wannier_init(D::Dict{Symbol,Any}, sp::String = "")
    # Print the header
    println("Generate input files for wannier90")

    # Extract necessary data from D
    latt  = D[:latt]
    kmesh = D[:kmesh]
    enk   = D[:enk]

    # Extract the nband parameter
    nband, _, _ = size(enk)

    # Try to prepare control parameters
    w90c = w90_make_ctrl(latt, nband)

    # Try to prepare projections
    proj = w90_make_proj()

    # Try to write w90.win
    #
    # Setup filename correctly. It depends on the spin.
    fwin = "w90" * sp * ".win"
    #
    # Write seedname.win
    open(fwin, "w") do fout
        w90_write_win(fout, w90c)
        w90_write_win(fout, proj)
        w90_write_win(fout, latt)
        w90_write_win(fout, kmesh)
    end
    #
    println("  > File $fwin is created")
end

"""
    wannier_exec(sp::String = ""; op::String = "")

Execute the wannier90 program, monitor the convergence progress, and
output the relevant information. The argument `sp` denotes the spin
component, while `op` specifies the running mode for wannier90.

See also: [`wannier_init`](@ref), [`wannier_save`](@ref).
"""
function wannier_exec(sp::String = ""; op::String = "")
    # Print the header
    println("Detect the runtime environment for wannier90")

    # Get the home directory of wannier90
    #
    # We can not guarantee that the wannier90 code is always installed
    # within the directory of pwscf.
    wannier90_home = query_dft("wannier90")
    println("  > Home directory for wannier90: ", wannier90_home)

    # Select suitable wannier90 program
    wannier90_exe = "$wannier90_home/wannier90.x"
    @assert isfile(wannier90_exe)
    println("  > Executable program is available: ", basename(wannier90_exe))

    # Assemble command
    seedname = "w90" * sp
    if op == "-pp"
        wannier90_cmd = split("$wannier90_exe $op $seedname", " ")
    else
        wannier90_cmd = split("$wannier90_exe $seedname", " ")
    end
    println("  > Assemble command: $(prod(x -> x * ' ', wannier90_cmd))")

    # Determine suitable output file
    finp = "w90" * sp * ".win"
    fout = "w90" * sp * ".out"
    println("  > Applying input file: $finp")
    println("  > Applying output file: $fout")

    # Print the header
    println("Launch the computational engine wannier90")

    # Create a task, but do not run it immediately
    t = @task begin
        run(pipeline(`$wannier90_cmd`, stdout = fout))
    end
    println("  > Create a task")

    # Launch it, the terminal output is redirected to `fout`.
    # Note that the task runs asynchronously. It will not block
    # the execution.
    schedule(t)
    println("  > Add the task to the scheduler's queue")
    println("  > Waiting ...")

    # To ensure that the task is executed
    while true
        sleep(2)
        istaskstarted(t) && break
    end

    # Wait for the wannier90 task to finish
    wait(t)
end

"""
    wannier_save(sp::String = ""; op::String = "")

Backup the output files of wannier90 if necessary.

See also: [`wannier_init`](@ref), [`wannier_exec`](@ref).
"""
function wannier_save(sp::String = ""; op::String = "")
    # Print the header
    println("Finalize the computational task")

    # Check the output files of wannier90
    if op == "-pp"
        fwan = "w90" * sp * ".nnkp"
        if isfile(fwan)
            println("  > File $fwan is created")
        else
            error("File $fwan is not created")
        end
    else
        fhr = "w90" * sp * "_hr.dat"
        fu = "w90" * sp * "_u.mat"
        fudis = "w90" * sp * "_u_dis.mat"
        #
        flist = (fhr, fu, fudis)
        for i in eachindex(flist)
            filename = flist[i]
            if isfile(filename)
                println("  > File $filename is created")
            else
                error("File $filename is not created")
            end
        end
    end
end

#=
### *Service Functions* : *Group B*
=#

"""
    w90_make_ctrl(latt:Lattice, nband::I64)

Try to make the control parameters for the `w90.win` file. The `latt`
object represent the crystallography information, and `nband` is the
number of Kohn-Sham states outputed by the dft code.

See also: [`w90_make_proj`](@ref).
"""
function w90_make_ctrl(latt::Lattice, nband::I64)
    # Create a dict, which will be returned.
    w90c = Dict{String,Any}()

    # Generate a dict, which defines a mapping from orbital definition
    # to number of orbitals.
    #
    # Perhaps you have to extend it to support your cases.
    orb_dict = Dict{String,I64}(
                 "s"   => 1,
                 "l=0" => 1,
                 "p"   => 3,
                 "l=1" => 3,
                 "d"   => 5,
                 "l=2" => 5,
                 "f"   => 7,
                 "l=3" => 7,
             )

    # Get number of wanniers, `num_wann`.
    #
    # Step 1, get the string for projection.
    sproj = get_d("sproj")
    @assert length(sproj) ≥ 2
    # The first element of sproj should specify the type of wannier
    # function. Now only the MLWF and SAWF modes are supported.
    @assert sproj[1] in ("mlwf", "sawf")
    #
    # Step 2, calculate num_wann.
    num_wann = 0
    # Go throught each element of sproj
    for i = 2:length(sproj)
        # Extract atomic symbol and orbital specification
        str_atm, str_orb = strip.(split(sproj[i], ":"))
        # We have to remove all white spaces in str_orb. It is necessary.
        str_orb = replace(str_orb, " " => "")
        # Extract the corresponding number of orbitals
        num_orb = orb_dict[str_orb]
        # Find out how many atoms are there in the lattice.
        num_atm = 0
        for j = 1:latt.natom
            if latt.atoms[j] == str_atm
                num_atm = num_atm + 1
            end
        end
        @assert num_atm ≥ 1
        # Update num_wann
        num_wann = num_wann + num_orb * num_atm
    end
    #
    # Step 3, store num_wann in the dict.
    w90c["num_wann"] = num_wann

    # Get number of bands, `num_bands`.
    @assert nband ≥ num_wann
    w90c["num_bands"] = nband

    # Deal with the disentanglement setup
    #
    # Step 1, get the string for disentanglement.
    window = get_d("window")
    @assert length(window) ≥ 2
    # The first element of window should specify the scheme for
    # disentanglement. Now only the exclude_bands and disentanglement
    # modes are supported.
    @assert window[1] in ("exc", "dis")
    #
    # Step 2, determine the disentanglement parameters and store them.
    if window[1] == "exc"
        w90c["exclude_bands"] = join(window[2:end], ", ")
    else
        if length(window) == 3
            w90c["dis_win_min"]  = window[2]
            w90c["dis_win_max"]  = window[3]
        elseif length(window) == 5
            w90c["dis_win_min"]  = window[2]
            w90c["dis_win_max"]  = window[3]
            w90c["dis_froz_min"] = window[4]
            w90c["dis_froz_max"] = window[5]
        else
            error("Wrong window's definition.")
        end
    end

    # Some additional but necessary parameters for wannier90
    #
    # We should write the hamiltonian
    w90c["write_hr"] = ".true."
    #
    # We should write the transform matrix
    w90c["write_u_matrices"] = ".true."
    #
    # Support symmetry-adapted wannier functions
    if sproj[1] == "sawf"
        w90c["site_symmetry"] = ".true."
    end
    #
    # Special treatment for spinors
    if get_d("lspinorb")
        w90c["num_wann"] = num_wann * 2
        w90c["spinors"] = ".true."
    end

    # Return the required object.
    return w90c
end

"""
    w90_make_proj()

Try to make the projection block for the `w90.win` file. We will not
check the validness of these projections here.

See also: [`w90_make_ctrl`](@ref).
"""
function w90_make_proj()
    # Create a string array to store the definitions for projections
    proj = String[]

    # Get the string for projection
    sproj = get_d("sproj")

    # Check the length of sproj
    @assert length(sproj) ≥ 2

    # The first element of sproj should specify the type of wannier
    # function. Now only the MLWF and SAWF modes are supported.
    @assert sproj[1] in ("mlwf", "sawf")

    # Transfer sproj to proj
    for i = 2:length(sproj)
        push!(proj, sproj[i])
    end

    # Return the desired array
    return proj
end

#=
### *Service Functions* : *Group C*
=#

"""
    w90_make_map()
"""
function w90_make_map()
end

"""
    w90_make_group(latt::Lattice, sp::String = "")

Try to read the `w90.nnkp` file, parse the `projections` block. Finally,
it will return arrays of PrTrait and PrGroup objects, which contain the
definitions of projectors. The argument `latt` is essential. It includes
the atomic coordinates for all lattice sites, which are quite useful to
distinguish these projectors. And the argument `sp` is optional. It is
used when the system is spin-polarized.

See also: [`PrTrait`](@ref), [`PrGroup`](@ref).
"""
function w90_make_group(latt::Lattice, sp::String = "")
    # Read and parse the `w90.nnkp` file
    #
    # Build the filename
    fnnkp = "w90" * sp * ".nnkp"
    #
    # Read it and figure out the projections block
    lines = readlines(fnnkp)
    ind = findfirst(x -> contains(x, "begin projections"), lines)
    @assert ind > 0
    start = ind + 1
    #
    # Create some arrays to store the projectors
    nproj = parse(I64, lines[start]) # Number of projectors
    coord = zeros(F64, nproj, 3)     # Atomic coordinates
    l_vec = zeros(I64, nproj) # Quantum number 𝑙
    m_vec = zeros(I64, nproj) # Quantum number 𝑚
    #
    # Parse the projectors
    for i = 1:nproj
        start = start + 1
        arr = line_to_array(lines[start])
        coord[i,:] = parse.(F64, arr[1:3])
        l_vec[i] = parse(I64, arr[4])
        m_vec[i] = parse(I64, arr[5])
        start = start + 1
    end

    # Try to build PrTrait struct. The raw information about projectors
    # should be encapsulated in it.
    #
    # Define all possible specifications for projectors. Be careful, the
    # following orbital sequence is different from that defined in the
    # constructor of PrTrait struct.
    spec = ("s",
            "pz", "px", "py",
            "dz2", "dxz", "dyz", "dx2-y2", "dxy",
            "fz3", "fxz2", "fyz2", "fz(x2-y2)", "fxyz", "fx(x2-3y2)", "fy(3x2-y2)"
    )
    #
    # Generate PrTrait struct one by one
    PT = PrTrait[]
    for i = 1:nproj
        # Try to figure out which site it is. We just compare the atomic
        # coordinates in coord with those saved in latt.coord.
        site = -1
        for j = 1:latt.natom
            if latt.coord[j,1:3] == coord[i,1:3]
                site = j # That is it
                break
            end
        end
        # Sanity check
        @assert site > 0
        # Get 𝑙, 𝑚, and desc.
        l = l_vec[i]
        m = m_vec[i]
        desc = spec[m + l*l]
        # Save the PrTrait struct
        push!(PT, PrTrait(site, l, m, desc))
    end

    # Try to split these projectors into groups.
    #
    # At first, we collect the tuple (site, l) for all projectors.
    site_l = Tuple[]
    for i in eachindex(PT)
        push!(site_l, (PT[i].site, PT[i].l))
    end
    #
    # Second, we figure out the unique (site, l) pairs. Here, we use the
    # union! function, which is much more efficient than the unique! one.
    union!(site_l)
    #
    # Third, we create a array of PrGroup struct (except for site and
    # l, most of its member variables need to be corrected). Note
    # that for a given PrGroup, the projectors indexed by PrGroup.Pr
    # should share the same site and l.
    PG = PrGroup[]
    for i in eachindex(site_l)
        push!(PG, PrGroup(site_l[i]...))
    end
    #
    # Fourth, for each PrGroup, we scan all of the projectors to find
    # out those with correct site and l; record their indices; and
    # save them at PrGroup.Pr array.
    for i in eachindex(PG)
        site, l = PG[i].site, PG[i].l
        PG[i].Pr = findall(x -> (x.site, x.l) == (site, l), PT)
    end
    #
    # Finally, check correctness
    @assert nproj == sum(x -> length(x.Pr), PG)

    # Print some useful information to check
    println("  > Number of projectors: ", nproj)
    for i in eachindex(PT)
        println("    [ PrTrait $i ]")
        println("      site -> ", PT[i].site)
        println("      l -> ", PT[i].l)
        println("      m -> ", PT[i].m)
        println("      desc -> ", PT[i].desc)
    end
    println("  > Number of groups: ", length(PG))
    for i in eachindex(PG)
        println("    [ PrGroup $i ]")
        println("      site -> ", PG[i].site)
        println("      l -> ", PG[i].l)
        println("      corr -> ", PG[i].corr)
        println("      shell -> ", PG[i].shell)
        println("      Pr -> ", PG[i].Pr)
    end

    # Return the desired arrays
    # Note: PG should be further setup at plo_group() function.
    return PT, PG
end

"""
    w90_make_window(PG::Array{PrGroup,1}, enk::Array{F64,2})
"""
function w90_make_window(PG::Array{PrGroup,1}, enk::Array{F64,2})
    # Extract the key parameters
    nband, nkpt = size(enk)

    # Initialize an array of PrWindow struct
    PW = PrWindow[]

    for p in eachindex(PG)
        bwin = (1, nband)

        kwin = zeros(I64, nkpt, 1, 2)
        fill!(view(kwin, :, :, 1), 1)
        fill!(view(kwin, :, :, 1), nband)
        push!(PW, PrWindow(kwin, bwin))

        # Print some useful information
        println("  > Create window $p: $bwin <--> ($(PW[p].bmin), $(PW[p].bmax))")
    end # END OF P LOOP
end

#=
### *Service Functions* : *Group D*
=#

"""
    w90_read_amat(sp::String = "")

Try to read and parse the `w90.amn` file to get the ``A_{mn}`` matrix,
which represents the projection of the Bloch states onto trail localized
orbitals. The argument `sp` denotes the spin component.
"""
function w90_read_amat(sp::String = "")
    # Build the filename
    famn = "w90" * sp * ".amn"

    # Read the `w90.amn` file
    lines = readlines(famn)
    @assert length(lines) ≥ 3

    # Determine the key parameters
    nband, nkpt, nproj = parse.(I64, line_to_array(lines[2]))

    # Create an array for ``A_{mn}``
    Amn = zeros(C64, nproj, nband, nkpt)

    # We use `start` to record the line index
    start = 2

    # Parse the data and fill in the Amn array
    for k = 1:nkpt
        for p = 1:nproj
            for b = 1:nband
                # Increase the counter
                start = start + 1
                # Convert string to array
                arr = line_to_array(lines[start])
                # Determine indices for bands, projectors, and 𝑘-points
                _b, _p, _k = parse.(I64, arr[1:3])
                # Sanity check
                @assert _b == b
                @assert _p == p
                @assert _k == k
                # Fill in the array
                _re, _im = parse.(F64, arr[4:5])
                Amn[p,b,k] = _re + im*_im
            end # END OF B LOOP
        end # END OF P LOOP
    end  # END OF K LOOP
    #
    @assert start - 2 == nband * nproj * nkpt

    # Return the desired array
    return Amn
end

"""
    w90_read_eigs(sp::String = "")

Try to read and parse the `w90.eig` file to get the band eigenvalues.
Note that the eigenvalues from `nscf.out` are not accurate enough. We
will use the data extracted from `w90.eig` instead. The argument `sp`
denotes the spin component.
"""
function w90_read_eigs(sp::String = "")
    # Build the filename
    feig = "w90" * sp * ".eig"

    # Read the `w90.eig` file
    lines = readlines(feig)
    @assert length(lines) ≥ 1

    # Determine the key parameters from the last line!
    nband, nkpt = parse.(I64, line_to_array(lines[end])[1:2])

    # Create an array for ``E_{nk}``
    eigs = zeros(F64, nband, nkpt)

    # We use `start` to record the line index
    start = 0

    # Parse the data and fill in the eigs array
    for k = 1:nkpt
        for b = 1:nband
            # Increase the counter
            start = start + 1
            # Convert string to array
            arr = line_to_array(lines[start])
            # Determine indices for bands and 𝑘-points
            _b, _k = parse.(I64, arr[1:2])
            # Sanity check
            @assert _b == b
            @assert _k == k
            # Fill in the array
            _v = parse(F64, arr[3])
            eigs[b,k] = _v
        end # END OF B LOOP
    end # END OF K LOOP
    #
    @assert start == nkpt * nband

    # Return the desired array
    return eigs
end

"""
    w90_read_hmat(sp::String = "")

Try to read and parse the `w90_hr.dat` file, return the hamiltonian
matrix in WF representation, the Wigner-Seitz grid points, and their
weights (degeneracies). The argument `sp` denotes the spin component.
"""
function w90_read_hmat(sp::String = "")
    # Build the filename
    fhr = "w90" * sp * "_hr.dat"

    # Read the `w90_hr.dat` file
    lines = readlines(fhr)
    @assert length(lines) ≥ 3

    # Determine the key parameters
    nproj = parse(I64, lines[2]) # Number of wannier functions
    nrpt = parse(I64, lines[3])  # Number of Wigner-Seitz grid points

    # Create arrays
    hamr = zeros(C64, nproj, nproj, nrpt) # Hamiltonian matrix
    rvec = zeros(I64, nrpt, 3) # Wigner-Seitz grids
    rdeg = zeros(I64, nrpt)    # Degeneracies or weights

    # We use `start` to record the line index
    start = 3

    # Read the degeneracies
    #
    # Determine how many lines are there for this block
    nrow = div(nrpt, 15)
    nrem = rem(nrpt, 15)
    @assert nrow ≥ 1
    #
    # Read the first nrow lines
    # There are 15 elements per line.
    for r = 1:nrow
        start = start + 1
        rs = (r - 1) * 15 + 1
        re = (r - 1) * 15 + 15
        rdeg[rs:re] = parse.(I64, line_to_array(lines[start]))
    end
    #
    # Read the final line
    if nrem > 0
        start = start + 1
        rs = nrow * 15 + 1
        re = nrpt
        @assert nrem == re - rs + 1
        rdeg[rs:re] = parse.(I64, line_to_array(lines[start]))
    end

    # Try to read the hamiltonian
    for r = 1:nrpt
        for i = 1:nproj
            for j = 1:nproj
                # Increase the counter
                start = start + 1
                # Convert string to line
                arr = line_to_array(lines[start])
                # Get Wigner-Seitz grid
                rvec[r,:] = parse.(I64, arr[1:3])
                # Get indices of wannier functions
                _j, _i = parse.(I64, arr[4:5])
                # Sanity check
                @assert _j == j
                @assert _i == i
                # Fill in the hamiltonian
                _re, _im = parse.(F64, arr[6:7])
                hamr[j, i, r] = _re + im * _im
            end # END OF J LOOP
        end # END OF I LOOP
    end # END OF R LOOP

    # Return the desired arrays
    return rdeg, rvec, hamr
end

"""
    w90_read_umat(sp::String = "")

Try to read and parse the `w90_u.mat` file, return the u-matrix, which
gives the unitary rotations from the optimal subspace to the optimally
smooth states. The argument `sp` denotes the spin component.

See also: [`w90_read_udis`](@ref).
"""
function w90_read_umat(sp::String = "")
    # Build the filename
    fu = "w90" * sp * "_u.mat"

    # Read the `w90_u.mat` file
    lines = readlines(fu)
    @assert length(lines) ≥ 3

    # Determine the key parameters
    nkpt, nproj = parse.(I64, line_to_array(lines[2])[1:2])

    # Create array
    umat = zeros(C64, nproj, nproj, nkpt)

    # We use `start` to record the line index
    start = 2

    # Try to read the u-matrix
    for k = 1:nkpt
        # Increase the counter
        start = start + 2 # Skip one empty line and one 𝑘-point line
        for j = 1:nproj
            for i = 1:nproj
                # Increase the counter
                start = start + 1
                # Parse the line and fill in the array
                _re, _im = parse.(F64, line_to_array(lines[start]))
                umat[i, j, k] = _re + im * _im
            end # END OF I LOOP
        end # END OF J LOOP
    end # END OF K LOOP

    # Return the desired array
    return umat
end

"""
    w90_read_udis(sp::String = "")

Try to read and parse the `w90_u_dis.mat` file. Return the udis-matrix,
which gives the nproj dimension optimal subspace from the original
bloch states. Actually, it is the transform matrix for disentanglement.
The argument `sp` denotes the spin component.

See also: [`w90_read_udis`](@ref).
"""
function w90_read_udis(sp::String = "")
    # Build the filename
    fu = "w90" * sp * "_u_dis.mat"

    # Read the `w90_u_dis.mat` file
    lines = readlines(fu)
    @assert length(lines) ≥ 3

    # Determine the key parameters
    nkpt, nproj, nband = parse.(I64, line_to_array(lines[2]))

    # Create array
    udis = zeros(C64, nband, nproj, nkpt)

    # We use `start` to record the line index
    start = 2

    # Try to read the udis-matrix
    for k = 1:nkpt
        # Increase the counter
        start = start + 2 # Skip one empty line and one 𝑘-point line
        for j = 1:nproj
            for i = 1:nband
                # Increase the counter
                start = start + 1
                # Parse the line and fill in the array
                _re, _im = parse.(F64, line_to_array(lines[start]))
                udis[i,j,k] = _re + im * _im
            end # END OF I LOOP
        end # END OF J LOOP
    end # END OF K LOOP

    # Return the desired array
    return udis
end

#=
### *Service Functions* : *Group E*
=#

"""
    w90_write_win(io::IOStream, w90c::Dict{String,Any})

Write control parameters into `w90.win`.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, w90c::Dict{String,Any})
    for key in keys(w90c)
        val = w90c[key]
        println(io, key, " = ", val)
    end
    # Add an empty line
    println(io)
end

"""
    w90_write_win(io::IOStream, proj::Array{String,1})

Write projection block into `w90.win`.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, proj::Array{String,1})
    println(io, "begin projections")
    #
    for i = 1:length(proj)
        println(io, proj[i])
    end
    #
    println(io, "end projections\n")
end

"""
    w90_write_win(io::IOStream, latt::Lattice)

Write crystallography information into `w90.win`.

See also: [`Lattice`](@ref), [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, latt::Lattice)
    # Extract parameters
    natom = latt.natom

    # Convert atomiclength (bohr) to angstrom
    # 1 bohr = 0.529177249 angstrom
    lvect = latt.lvect * (latt.scale * 0.529177249)

    # Print atoms_frac block
    println(io, "begin atoms_frac")
    #
    for i = 1:natom
        @printf(io, "%4s%12.8f%12.8f%12.8f\n", latt.atoms[i], latt.coord[i,:]...)
    end
    #
    println(io, "end atoms_frac\n")

    # Print unit_cell_cart block
    println(io, "begin unit_cell_cart")
    #
    for i = 1:3
        @printf(io, "%12.8f%12.8f%12.8f\n", lvect[i,:]...)
    end
    #
    println(io, "end unit_cell_cart\n")
end

"""
    w90_write_win(io::IOStream, kmesh::Array{F64,2})

Write 𝑘-mesh block into `w90.win`.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, kmesh::Array{F64,2})
    # Extract parameters
    nkpt, ndir = size(kmesh)

    # Sanity check
    @assert ndir == 3

    # Write the block for k-points
    ndiv = ceil(I64, nkpt^(1/3))
    #
    @printf(io, "mp_grid :%3i%3i%3i\n", ndiv, ndiv, ndiv)
    println(io, "begin kpoints")
    #
    for k = 1:nkpt
        @printf(io, "%12.8f%12.8f%12.8f\n", kmesh[k,:]...)
    end
    #
    println(io, "end kpoints\n")
end

#=
### *Service Functions* : *Group X*
=#

"""
    pw2wan_init(case::String, sp::String = "")

Check the runtime environment of `pw2wannier90`, prepare necessary input
files (`case.pw2wan`). The argument `case` means the prefix for `pwscf`,
and `sp` determines the spin component which can be empty string, `up`,
or `dn`.

See also: [`pw2wan_exec`](@ref), [`pw2wan_save`](@ref).
"""
function pw2wan_init(case::String, sp::String = "")
    # Print the header
    println("Generate input files for pw2wannier90")

    # Try to create a PWNamelist object.
    #
    # Setup name of namelist. It is always fixed.
    name = "inputpp"
    #
    # Create an empty dict
    NLData = Dict{AbstractString,Any}()
    #
    # Update the dict
    NLData["outdir"] = "'./'"
    NLData["prefix"] = "'$case'"
    if sp == ""
        NLData["seedname"] = "'w90'"
        NLData["spin_component"] = "'none'"
    elseif sp == "up"
        NLData["seedname"] = "'w90up'"
        NLData["spin_component"] = "'up'"
    elseif sp == "dn"
        NLData["seedname"] = "'w90dn'"
        NLData["spin_component"] = "'down'"
    end
    NLData["write_mmn"] = ".true."
    NLData["write_amn"] = ".true."
    NLData["write_dmn"] = ".true."
    NLData["write_unk"] = ".false."
    #
    # Create PWNamelist
    PWN = PWNamelist(name, NLData)

    # Try to write case.pw2wan
    fwan = case * sp * ".pw2wan"
    open(fwan, "w") do fout
        write(fout, PWN) # This write function is defined in pwscf.jl
    end
    #
    println("  > File $fwan is created")
end

"""
    pw2wan_exec(case::String, sp::String = "")

Execute the pw2wannier90 program, monitor the convergence progress, and
output the relevant information. Here, `case` means the prefix for input
files, and `sp` is `up` or `dn`.

See also: [`pw2wan_init`](@ref), [`pw2wan_save`](@ref).
"""
function pw2wan_exec(case::String, sp::String = "")
    # Print the header
    println("Detect the runtime environment for pw2wannier90")

    # Get the home directory of pw2wannier90
    # It is actually the same with that of pwscf.
    pwscf_home = query_dft("pwscf")
    println("  > Home directory for pw2wannier90: ", pwscf_home)

    # Select suitable pw2wannier90 program
    pwscf_exe = "$pwscf_home/pw2wannier90.x"
    @assert isfile(pwscf_exe)
    println("  > Executable program is available: ", basename(pwscf_exe))

    # Assemble command
    pwscf_cmd = split("$pwscf_exe", " ")
    println("  > Assemble command: $(prod(x -> x * ' ', pwscf_cmd))")

    # Determine suitable input and output files
    finp = case * sp * ".pw2wan"
    fout = "pw2wan.out"
    println("  > Applying input file: $finp")
    println("  > Applying output file: $fout")

    # Print the header
    println("Launch the computational engine pw2wannier90")

    # Create a task, but do not run it immediately
    t = @task begin
        run(pipeline(`$pwscf_cmd`, stdin = finp, stdout = fout))
    end
    println("  > Create a task")

    # Launch it, the terminal output is redirected to `fout`.
    # Note that the task runs asynchronously. It will not block
    # the execution.
    schedule(t)
    println("  > Add the task to the scheduler's queue")
    println("  > Waiting ...")

    # To ensure that the task is executed
    while true
        sleep(2)
        istaskstarted(t) && break
    end

    # Wait for the pw2wannier90 task to finish
    wait(t)
end

"""
    pw2wan_save(sp::String = "")

Backup the output files of `pw2wannier90` if necessary. The argument `sp`
specifies the spin component. It could be empty string, `up`, or `dn`.

See also: [`pw2wan_init`](@ref), [`pw2wan_exec`](@ref).
"""
function pw2wan_save(sp::String = "")
    # Print the header
    println("Finalize the computational task")

    # Check the output files of pw2wannier90
    #
    # Determine filenames
    seedname = "w90" * sp
    famn = seedname * ".amn"
    fmmn = seedname * ".mmn"
    feig = seedname * ".eig"
    fdmn = seedname * ".dmn"
    fsym = seedname * ".sym"
    #
    # Determine list of files
    flist = [famn, fmmn, feig]
    #
    # If the symmetry-adapted wannier function mode is chosen, more
    # output files are created by pw2wannier90.
    sproj = get_d("sproj")
    if sproj[1] == "sawf"
        push!(flist, fdmn)
        push!(flist, fsym)
    end
    #
    # Check the availability of the output files
    for i in eachindex(flist)
        filename = flist[i]
        if isfile(filename)
            println("  > File $filename is created")
        else
            error("File $filename is not created")
        end
    end
end
