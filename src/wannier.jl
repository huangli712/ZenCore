#
# Project : Pansy
# Source  : wannier.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/08/30
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
    wannier_init()
"""
function wannier_init(D::Dict{Symbol,Any})
    latt  = D[:latt] 
    kmesh = D[:kmesh]

    case  = latt._case

    open("$case.win", "w") do fout
        w90_write_win(fout, kmesh)
    end
end

function wannier_exec()
end

function wannier_save()
end

#=
### *Service Functions* : *Group X*
=#

function pw2wan_init()
end

function pw2wan_exec()
end

function pw2wan_save()
end

function w90_write_win(io::IOStream, kmesh::Array{F64,2})
    println("begin kpoints")
    println("end kpoints")
end
