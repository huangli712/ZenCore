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
end

#=
### *Service Functions* : *Group A*
=#

function wannier_init()
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
