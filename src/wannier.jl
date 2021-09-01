#
# Project : Pansy
# Source  : wannier.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/09/01
#

#=
### *Driver Functions*
=#

function wannier_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})
    # Print the header
    println("Adaptor : WANNIER")
    println("Try to process the Kohn-Sham dataset")
    println("Current directory: ", pwd())

    wannier_init(D, ai)

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
    wannier_init(D::Dict{Symbol,Any}, ai::Array{Impurity,1})

Try to generate the `w90.win` file, which is the essential input for
the `wannier90` code.

See also: [`wannier_exec`](@ref), [`wannier_save`](@ref).
"""
function wannier_init(D::Dict{Symbol,Any}, ai::Array{Impurity,1})
    # Extract key parameters and arrays
    latt  = D[:latt] 
    kmesh = D[:kmesh]

    w90c = w90_build_ctrl()
    w90_build_proj()

    open("w90.win", "w") do fout
        w90_write_win(fout, w90c)
        w90_write_win(fout, latt)
        w90_write_win(fout, kmesh)
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
    w90_build_ctrl()

"""
function w90_build_ctrl()
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
        str_orb = strip(split(sproj[i], ":")[2])
        num_orb = orb_dict[str_orb]
        num_wann = num_wann + num_orb
    end
    println(num_wann)

    return w90c
end

"""
    w90_build_proj()
"""
function w90_build_proj()
end

"""
    w90_write_win(io::IOStream, w90c::Dict{String,Any})

Write control parameters into case.win.

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
    w90_write_win(io::IOStream, latt::Lattice)

Write crystallography information into case.win.
 
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

Write the block for k-points into case.win.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, kmesh::Array{F64,2})
    # Extract parameters
    nkpt, ndir = size(kmesh)

    # Sanity check
    @assert ndir == 3

    # Write the block for k-points
    println(io, "begin kpoints")
    #
    for k = 1:nkpt
        @printf(io, "%12.8f%12.8f%12.8f\n", kmesh[k,:]...)
    end
    #
    println(io, "end kpoints\n")
end
