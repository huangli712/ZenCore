#
# Project : Pansy
# Source  : mixer.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/01
#

"""
    mixer_sigma(it::IterInfo)

See also: [`mixer_core`](@ref).
"""
function mixer_sigma(it::IterInfo)
    # Print the log
    println("Mixer : Sigma")

    # Print blank line for better visualization
    println()
end

"""
    mixer_delta(it::IterInfo)

See also: [`mixer_core`](@ref).
"""
function mixer_delta(it::IterInfo)
    # Print the log
    println("Mixer : Delta")

    # Print blank line for better visualization
    println()
end

"""
    mixer_eimpx(it::IterInfo)

Try to mix the local impurity levels.

See also: [`mixer_core`](@ref).
"""
function mixer_eimpx(it::IterInfo)
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
