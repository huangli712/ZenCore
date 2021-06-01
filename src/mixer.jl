#
# Project : Pansy
# Source  : mixer.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/02
#

"""
    mixer_sigma(it::IterInfo, ai::Array{Impurity,1})

Try to mix the self-energy functions and then use the mixed values to
update the `dmft1/sigma.bare` file.

See also: [`mixer_core`](@ref).
"""
function mixer_sigma(it::IterInfo, ai::Array{Impurity,1})
    # Print the log
    println("Mixer : Sigma")

    # Get current dmft loop
    cycle = it.dmft_cycle

    # Get current iteration
    curr = it.dmft1_iter

    # Get previous iteration
    prev = it.dmft1_iter - 1
    @assert prev > 0

    # Determine filenames for hybridization functions
    fcurr = "dmft1/sigma.bare.$cycle.$curr"
    fprev = "dmft1/sigma.bare.$cycle.$prev"

    # Print blank line for better visualization
    println()
end

"""
    mixer_delta(it::IterInfo, ai::Array{Impurity,1})

Try to mix the hybridization functions and then use the mixed values to
update the `dmft1/dmft.hyb_l` file.

See also: [`mixer_core`](@ref).
"""
function mixer_delta(it::IterInfo, ai::Array{Impurity,1})
    # Print the log
    println("Mixer : Delta")

    # Get current dmft loop
    cycle = it.dmft_cycle

    # Get current iteration
    curr = it.dmft1_iter

    # Get previous iteration
    prev = it.dmft1_iter - 1
    @assert prev > 0

    # Determine filenames for hybridization functions
    fcurr = "dmft1/dmft.hyb_l.$cycle.$curr"
    fprev = "dmft1/dmft.hyb_l.$cycle.$prev"

    # Check whether these files are available
    @assert isfile(fcurr) && isfile(fprev)

    # Read in the hybridization functions (previous and current)
    fcurr, Dcurr = read_delta(ai, fcurr)
    fprev, Dprev = read_delta(ai, fprev)
    @assert size(Dcurr) == size(Dprev) && size(fcurr) == size(fprev)

    # Mix the hybridization functions
    Dnew = Dcurr * get_m("mixer") + Dprev * (1.0 - get_m("mixer"))

    # Write the new hybridization functions into `dmft1/dmft.hyb_l`
    write_delta(fcurr, Dnew, ai, "dmft1/dmft.hyb_l")

    # Print blank line for better visualization
    println()
end

"""
    mixer_eimpx(it::IterInfo, ai::Array{Impurity,1})

Try to mix the local impurity levels and then use the mixed value to
update the `dmft1/dmft.eimpx` file.

See also: [`mixer_core`](@ref).
"""
function mixer_eimpx(it::IterInfo, ai::Array{Impurity,1})
    # Print the log
    println("Mixer : Eimpx")

    # Get current dmft loop
    cycle = it.dmft_cycle

    # Get current iteration
    curr = it.dmft1_iter

    # Get previous iteration
    prev = it.dmft1_iter - 1
    @assert prev > 0

    # Determine filenames for local impurity levels
    fcurr = "dmft1/dmft.eimpx.$cycle.$curr"
    fprev = "dmft1/dmft.eimpx.$cycle.$prev"

    # Check whether these files are available
    @assert isfile(fcurr) && isfile(fprev)

    # Read in the local impurity levels (previous and current)
    Ecurr = read_eimpx(ai, fcurr)
    Eprev = read_eimpx(ai, fprev)
    @assert size(Ecurr) == size(Eprev)

    # Mix the local impurity levels
    Enew = Ecurr * get_m("mixer") + Eprev * (1.0 - get_m("mixer"))

    # Write the new local impurity levels into `dmft1/dmft.eimpx`
    write_eimpx(Enew, ai, "dmft1/dmft.eimpx")

    # Print blank line for better visualization
    println()
end

"""
    mixer_gamma(it::IterInfo)

See also: [`mixer_core`](@ref).
"""
function mixer_gamma(it::IterInfo)
    # Print the log
    println("Mixer : Gamma")

    # Print blank line for better visualization
    println()
end
