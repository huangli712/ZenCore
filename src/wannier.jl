#
# Project : Pansy
# Source  : wannier.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/09/04
#

#=
### *Driver Functions*
=#

"""
    wannier_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})
"""
function wannier_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})
    # Print the header
    println("Adaptor : WANNIER")
    println("Try to process the Kohn-Sham dataset")
    println("Current directory: ", pwd())

    # Extract key parameters
    case = get_c("case")
    sp = get_d("lspins")

    # W01: Generate seedname.win.
    wannier_init(D, sp)

    # W02: Generate case.pw2wan.
    pw2wan_init(case, sp)

    # W03: Execute wannier90 to generate w90.nnkp.
    if sp # For spin-polarized system
        wannier_exec("up", op = "-pp")
        wannier_save("up", op = "-pp")
        wannier_exec("dn", op = "-pp")
        wannier_save("dn", op = "-pp")
    else # For spin-unpolarized system
        wannier_exec(op = "-pp")
        wannier_save(op = "-pp")
    end
    
    # W04: Execute pw2wannier90 to generate w90.amn, w90.mmn, etc.
    if sp # For spin-polarized system
        pw2wan_exec(case, "up")
        pw2wan_save("up")
        pw2wan_exec(case, "dn")
        pw2wan_save("dn")
    else # For spin-unpolarized system
        pw2wan_exec(case)
        pw2wan_save()
    end

    # W05: Execute wannier90 again to generate wannier functions.
    if sp # For spin-polarized system
        wannier_exec("up")
        wannier_save("up")
        wannier_exec("dn")
        wannier_save("dn")
    else # For spin-unpolarized system
        wannier_exec()
        wannier_save()
    end
end

#=
### *Service Functions* : *Group A*
=#

"""
    wannier_init(D::Dict{Symbol,Any}, sp::Bool = false)

Try to generate the `w90.win` file, which is the essential input for
the `wannier90` code. Here, we always use `w90` as the seedname. If
the system is spin polarized (`sp = true`), then the seednames will
be `w90up` and `w90dn`, respectively.

See also: [`wannier_exec`](@ref), [`wannier_save`](@ref).
"""
function wannier_init(D::Dict{Symbol,Any}, sp::Bool = false)
    # Print the header
    println("Generate input files for wannier90")

    # Extract necessary data from D
    latt  = D[:latt] 
    kmesh = D[:kmesh]
    enk   = D[:enk]

    # Extract the nband parameter
    nband, _, _ = size(enk)

    # Try to prepare control parameters
    w90c = w90_build_ctrl(latt, nband)

    # Try to prepare projections
    proj = w90_build_proj()

    # Try to write w90.win
    #
    # Setup filename correctly. It depends on the spin.
    fwin = "w90.win"
    if sp # Spin polarized system
        fwin = "w90up.win"
    end
    #
    # Write seedname.win
    open(fwin, "w") do fout
        w90_write_win(fout, w90c)
        w90_write_win(fout, proj)
        w90_write_win(fout, latt)
        w90_write_win(fout, kmesh)
    end
    println("  > File $fwin is created")
    #
    # Copy seedname.win if necessary
    if sp # Spin polarized system
        cp(fwin, "w90dn.win")
        println("  > File w90dn.win is created")
    end
end

"""
    wannier_exec(sp::String = ""; op::String = "")

Execute the wannier90 program, monitor the convergence progress, and
output the relevant information.

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
    pw2wan_init(case::String, sp::Bool = false)

Check the runtime environment of `pw2wannier90`, prepare necessary input
files (`case.pw2wan`).

See also: [`pw2wan_exec`](@ref), [`pw2wan_save`](@ref).
"""
function pw2wan_init(case::String, sp::Bool = false)
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
    NLData["seedname"] = "'w90'"
    NLData["spin_component"] = "'none'"
    NLData["write_mmn"] = ".true."
    NLData["write_amn"] = ".true."
    NLData["write_dmn"] = ".true."
    #
    # Create PWNamelist
    PWN = PWNamelist(name, NLData)

    # Try to write case.pw2wan
    # For spin unpolarized system
    if !sp
        fwan = "$case.pw2wan"
        open(fwan, "w") do fout
            write(fout, PWN)
        end
        #
        println("  > File $fwan is created")
    # For spin polarized system
    else
        # Spin up case
        fwan = case * "up.pw2wan"
        PWN["seedname"] = "'w90up'"
        PWN["spin_component"] = "'up'"
        open(fwan, "w") do fout
            write(fout, PWN)
        end
        #
        println("  > File $fwan is created")
        #
        # Spin down case
        fwan = case * "dn.pw2wan"
        PWN["seedname"] = "'w90dn'"
        PWN["spin_component"] = "'dn'"
        open(fwan, "w") do fout
            write(fout, PWN)
        end
        #
        println("  > File $fwan is created")
    end
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

Backup the output files of pw2wannier90 if necessary.

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
    # output files are created.
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

#=
### *Service Functions* : *Group C*
=#

"""
    w90_build_ctrl(latt:Lattice, nband::I64)

Try to make the control parameters for the `w90.win` file. The `latt`
object represent the crystallography information, and `nband` is the
number of Kohn-Sham states outputed by the dft code.

See also: [`w90_build_proj`](@ref).
"""
function w90_build_ctrl(latt::Lattice, nband::I64)
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

    # Return the required object.
    return w90c
end

"""
    w90_build_proj()

Try to make the projection block for the `w90.win` file.

See also: [`w90_build_ctrl`](@ref).
"""
function w90_build_proj()
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
### *Service Functions* : *Group D*
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
    #
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

Write k-mesh block into `w90.win`.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, kmesh::Array{F64,2})
    # Extract parameters
    nkpt, ndir = size(kmesh)

    # Sanity check
    @assert ndir == 3

    # Write the block for k-points
    ndiv = ceil(I64, nkpt^(1/3))
    @printf(io, "mp_grid :%3i%3i%3i\n", ndiv, ndiv, ndiv)
    println(io, "begin kpoints")
    #
    for k = 1:nkpt
        @printf(io, "%12.8f%12.8f%12.8f\n", kmesh[k,:]...)
    end
    #
    println(io, "end kpoints\n")
end
