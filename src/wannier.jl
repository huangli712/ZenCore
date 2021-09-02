#
# Project : Pansy
# Source  : wannier.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/09/02
#

#=
### *Driver Functions*
=#

function wannier_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})
    # Print the header
    println("Adaptor : WANNIER")
    println("Try to process the Kohn-Sham dataset")
    println("Current directory: ", pwd())

    wannier_init(D)

    wannier_exec()
    wannier_save()

    pw2wan_init()
    pw2wan_exec()
    pw2wan_save()

    wannier_exec()
    wannier_save()
end

#=
### *Service Functions* : *Group A*
=#

"""
    wannier_init(D::Dict{Symbol,Any})

Try to generate the `w90.win` file, which is the essential input for
the `wannier90` code. Here, we always use `w90` as the seedname. If
the system is spin polarized, then the seednames will be `w90up` and
`w90dn`, respectively.

See also: [`wannier_exec`](@ref), [`wannier_save`](@ref).
"""
function wannier_init(D::Dict{Symbol,Any})
    # Extract necessary data from D
    latt  = D[:latt] 
    kmesh = D[:kmesh]
    enk = D[:enk]

    # Extract number of spin
    _, _, nspin = size(enk)

    # Try to prepare control parameters
    w90c = w90_build_ctrl(latt, enk)

    # Try to prepare projections
    proj = w90_build_proj()

    # Try to write w90.win
    #
    # Setup filename correctly. It depends on the spin.
    fwin = "w90.win"
    if nspin == 2
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
    #
    # Copy seedname.win if necessary
    if nspin == 2
        cp(fwin, "w90dn.win")
    end

    sorry()
end

"""
    wannier_exec()
"""
function wannier_exec()
end

"""
    wannier_save()
"""
function wannier_save()
end

#=
### *Service Functions* : *Group X*
=#

"""
    pw2wan_init()
"""
function pw2wan_init()
end

"""
    pw2wan_exec()
"""
function pw2wan_exec()
end

"""
    pw2wan_save()
"""
function pw2wan_save()
end

"""
    w90_build_ctrl(latt:Lattice, enk::Array{F64,3})

Try to make the control parameters for the `w90.win` file.

See also: [`w90_build_proj`](@ref).
"""
function w90_build_ctrl(latt::Lattice, enk::Array{F64,3})
    w90c = Dict{String,Any}()

    orb_dict = Dict{String,I64}(
                 "s" => 1,
                 "p" => 3,
                 "d" => 5,
                 "f" => 7,
             )

    # Get number of wanniers, `num_wann`
    sproj = get_d("sproj")
    @assert length(sproj) ≥ 2
    @assert sproj[1] in ("mlwf", "sawf")
    num_wann = 0
    for i = 2:length(sproj)
        str_atm, str_orb = strip.(split(sproj[i], ":"))
        num_atm = 0
        for j = 1:latt.natom
            if latt.atoms[j] == str_atm
                num_atm = num_atm + 1
            end
        end
        @assert num_atm >= 1
        num_orb = orb_dict[str_orb]
        num_wann = num_wann + num_orb * num_atm
    end
    w90c["num_wann"] = num_wann

    # Get number of bands, `num_bands`
    num_bands, _, _ = size(enk)
    w90c["num_bands"] = num_bands

    #
    window = get_d("window")
    @assert length(window) >= 2
    @assert window[1] in ("exc", "dis")
    if window[1] == "exc"
        w90c["exclude_bands"] = join(window[2:end], ", ")
    else
        if length(window) == 3
            w90c["dis_win_min"] = window[2]
            w90c["dis_win_max"] = window[3]
        elseif length(window) == 5
            w90c["dis_win_min"] = window[2]
            w90c["dis_win_max"] = window[3]
            w90c["dis_froz_min"] = window[4]
            w90c["dis_froz_max"] = window[5]
        else
            error("Wrong window.")
        end
    end

    return w90c
end

"""
    w90_build_proj()

Try to make the projection block for the `w90.win` file.

See also: [`w90_build_ctrl`](@ref).
"""
function w90_build_proj()
    proj = String[]
    sproj = get_d("sproj")
    @assert length(sproj) ≥ 2
    @assert sproj[1] in ("mlwf", "sawf")
    for i = 2:length(sproj)
        push!(proj, sproj[i])
    end
    return proj
end

"""
    w90_write_win(io::IOStream, w90c::Dict{String,Any})

Write control parameters into w90.win.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, w90c::Dict{String,Any})
    for key in keys(w90c)
        val = w90c[key]
        println(io, "$key = $val")
    end
    println(io)
end

"""
    w90_write_win(io::IOStream, proj::Array{String,1})

Write projection block into w90.win.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, proj::Array{String,1})
    println(io, "begin projections")
    for i = 1:length(proj)
        println(io, proj[i])
    end
    println(io, "end projections\n")
end

"""
    w90_write_win(io::IOStream, latt::Lattice)

Write crystallography information into w90.win.
 
See also: [`Lattice`](@ref), [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, latt::Lattice)
    # Extract parameters
    natom = latt.natom

    # Convert atomiclength (bohr) to angstrom
    lvect = latt.lvect * (latt.scale * 0.52918)

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

Write k-mesh block into w90.win.

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
