#
# Project : Pansy
# Source  : sigma.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/31
#

"""
    sigma_reset(ai::Array{Impurity,1})

Create initial self-energy functions and write them to `sigma.bare`.

See also: [`sigma_dcount`](@ref).
"""
function sigma_reset(ai::Array{Impurity,1})
    # Print the log
    println("Sigma : Reset")

    # Extract some necessary parameters
    axis = get_m("axis")
    nmesh = get_m("nmesh")
    beta = get_m("beta")
    nsite = get_i("nsite")
    nspin = ( get_d("lspins") ? 2 : 1 )
    @assert nsite == length(ai)

    # Create frequency mesh
    fmesh = zeros(F64, nmesh)
    if axis === 1 # Imaginary axis
        for i = 1:nmesh
            fmesh[i] = (2 * i - 1) * pi / beta
        end
        println("  Create Matsubara frequency mesh")
    else # Real axis
        sorry()
    end

    # Create self-energy functions
    #
    # Initialize an array for self-energy functions
    SA = Array{C64,4}[]
    #
    # Go through the impurity problems
    for i = 1:nsite
        # Get the dimension of impurity problem
        nband = ai[i].nband

        # Create a temporary array for self-energy function
        S = zeros(C64, nband, nband, nmesh, nspin)

        # Push S into SA to save it
        push!(SA, S)
    end
    println("  Create local self-energy functions")

    # Write self-energy functions to sigma.bare
    open("dmft1/sigma.bare", "w") do fout
        # Write the header
        println(fout, "# File: sigma.bare")
        println(fout, "# Data: bare self-energy functions")
        println(fout)
        println(fout, "axis  -> $axis")
        println(fout, "beta  -> $beta")
        println(fout, "nsite -> $nsite")
        println(fout, "nmesh -> $nmesh")
        println(fout, "nspin -> $nspin")
        for i = 1:nsite
            println(fout, "ndim$i -> $(ai[i].nband)")
        end
        println(fout)

        # Write the body
        # Go through each quantum impurity problem
        for i = 1:nsite
            for s = 1:nspin
                println(fout, "# site: $i spin: $s")
                for m = 1:nmesh
                    # Write frequency grid
                    @printf(fout, "%4s %6i %16.12f\n", "w:", m, fmesh[m])
                    # Write self-energy function data
                    # There are 2 columns and nband * nband rows
                    for p = 1:ai[i].nband
                        for q = 1:ai[i].nband
                            x = SA[i][q, p, m, s]
                            @printf(fout, "%16.12f %16.12f\n", real(x), imag(x))
                        end
                    end
                end # END OF M LOOP
                # Write separator for each frequency point
                println(fout)
            end # END OF S LOOP
        end # END OF I LOOP
    end
    println("  Write self-energy functions into dmft1/sigma.bare")

    # Print blank line for better visualization
    println()
end

"""
    sigma_dcount(it::IterInfo, ai::Array{Impurity,1})

Calculate double counting terms for local self-energy functions and
write them to `sigma.dc`.

See also: [`sigma_reset`](@ref).
"""
function sigma_dcount(it::IterInfo, ai::Array{Impurity,1})
    # Print the log
    println("Sigma : Dcount")

    # Extract some necessary parameters
    nsite = get_i("nsite")
    nspin = ( get_d("lspins") ? 2 : 1 )
    @assert nsite == length(ai)

    # Create double counting terms for self-energy functions
    #
    # Initialize an array for dc
    DCA = Array{F64,3}[]
    #
    # Go through the impurity problems and calculate dc
    for i = 1:nsite
        # Get interaction parameters
        U = ai[i].upara
        J = ai[i].jpara

        # Get nominal occupation number
        N = get_i("occup")[i]

        # Get realistic occupation number
        GetNimpx(ai[i])
        occup = ai[i].occup

        # Get number of orbitals
        nband = ai[i].nband

        # Create a temporary array for the double counting terms
        DC = zeros(F64, nband, nband, nspin)

        # Choose suitable double counting scheme
        @cswitch get_m("dcount") begin
            # Fully localized limit scheme with fixed occupation number
            @case "fll1"
                sigdc = cal_dc_fll(U, J, N)
                fill!(DC, sigdc)
                break

            # Fully localized limit scheme with dynamic occupation number
            @case "fll2"
                sigdc = cal_dc_fll(U, J, occup)
                fill!(DC, sigdc)
                break

            # Around mean-field scheme
            @case "amf"
                sorry()
                break

            # Exact double counting scheme
            @case "exact"
                sorry()
                break
        end

        # Special treatment for the first iteration
        if it.dmft_cycle <= 1 && it.dmft1_iter <= 1
            fill!(DC, 0.0)
        end

        # Push DC into DCA to save it
        push!(DCA, DC)
    end # END OF I LOOP
    println("  Create double counting terms: $(get_m("dcount"))")

    # Write double counting terms to sigma.dc
    open("dmft1/sigma.dc", "w") do fout
        # Write the header
        println(fout, "# File: sigma.dc")
        println(fout, "# Data: double counting terms")
        println(fout)
        println(fout, "nsite -> $nsite")
        println(fout, "nspin -> $nspin")
        for i = 1:nsite
            println(fout, "ndim$i -> $(ai[i].nband)")
        end
        println(fout)

        # Write the body
        # Go through each impurity problem
        for i = 1:nsite
            for s = 1:nspin
                println(fout, "# site: $i spin: $s")
                # There are 2 columns and nband * nband rows
                # The double counting terms are assumed to be complex
                # numbers with zero imaginary parts.
                for p = 1:ai[i].nband
                    for q = 1:ai[i].nband
                        # Only the diagonal elements are useful
                        if p == q
                            @printf(fout, "%16.12f %16.12f\n", DCA[i][q, p, s], 0.0)
                        else
                            @printf(fout, "%16.12f %16.12f\n", 0.0, 0.0)
                        end
                    end # END OF Q LOOP
                end # END OF P LOOP
                println(fout)
            end # END OF S LOOP
        end # END OF I LOOP
    end
    println("  Write double counting terms into: dmft1/sigma.dc")

    # Print blank line for better visualization
    println()
end

"""
    sigma_split(ai::Array{Impurity,1})

Split the hybridization functions (and local impurity levels) and then
distribute them into the `impurity.i` folder.

See also: [`sigma_gather`](@ref).
"""
function sigma_split(ai::Array{Impurity,1})
    # Print the log
    println("Sigma : Split")

    # Split the hybridization functions
    split_hyb_l(ai)

    # Split the local impurity levels
    split_eimpx(ai)

    # Print blank line for better visualization
    println()
end

"""
    sigma_gather(ai::Array{Impurity,1})

Gather the self-energy functions (or similar local functions) from the
`impurity.i` folder and then combine them into a single `sigma.bare` file.

See also: [`sigma_split`](@ref).
"""
function sigma_gather(ai::Array{Impurity,1})
    # Print the log
    println("Sigma : Gather")

    # Extract some necessary parameters
    axis = get_m("axis")
    nmesh = get_m("nmesh")
    beta = get_m("beta")
    nsite = get_i("nsite")
    nspin = ( get_d("lspins") ? 2 : 1 )
    @assert nsite == length(ai)

    # Create empty array for self-energy functions
    SA = Array{C64,4}[]

    # Declare frequency mesh
    fmesh = nothing

    # Go through each quantum impurity problems
    for t = 1:nsite

        # Extract the frequency mesh and self-energy function
        fmesh, sig_l = GetSig_l(ai[t])
        println("  Read self-energy functions for impurity: $t")

        # Extract and verify the dimensional parameters
        _, _b, _m, _s = size(sig_l)
        @assert _b == ai[t].nband
        @assert _m == nmesh

        # Store sig_l in SA
        push!(SA, sig_l)
    end # END OF T LOOP

    # Now the self-energy functions for all quantum impurity problems and
    # the corresponding frequency mesh are ready. We are going to write
    # them into the `sigma.bare` file.

    # Write self-energy functions to sigma.bare
    open("dmft1/sigma.bare", "w") do fout
        # Write the header
        println(fout, "# File: sigma.bare")
        println(fout, "# Data: bare self-energy functions")
        println(fout)
        println(fout, "axis  -> $axis")
        println(fout, "beta  -> $beta")
        println(fout, "nsite -> $nsite")
        println(fout, "nmesh -> $nmesh")
        println(fout, "nspin -> $nspin")
        for i = 1:nsite
            println(fout, "ndim$i -> $(ai[i].nband)")
        end
        println(fout)

        # Write the body
        # Go through each quantum impurity problem
        for i = 1:nsite
            for s = 1:nspin
                println(fout, "# site: $i spin: $s")
                for m = 1:nmesh
                    # Write frequency grid
                    @printf(fout, "%4s %6i %16.12f\n", "w:", m, fmesh[m])
                    # Write self-energy function data
                    # There are 2 columns and nband * nband rows
                    for p = 1:ai[i].nband
                        for q = 1:ai[i].nband
                            x = SA[i][q, p, m, s]
                            @printf(fout, "%16.12f %16.12f\n", real(x), imag(x))
                        end
                    end
                end # END OF M LOOP
                # Write separator for each frequency point
                println(fout)
            end # END OF S LOOP
        end # END OF I LOOP
    end
    println("The local self-energy functions are written to sigma.bare")

    # Print blank line for better visualization
    println()
end

"""
    cal_dc_fll(U::F64, J::F64, N::F64)

Evaluate the double counting term by the fully localized limit scheme.

See also: [`cal_dc_amf`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_fll(U::F64, J::F64, N::F64)
    U * ( N - 0.5 ) - J / 2.0 * ( N - 1.0 )
end

"""
    cal_dc_amf(U::F64, J::F64, N::F64)

Evaluate the double counting term by the around mean-field scheme.

See also: [`cal_dc_fll`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_amf(U::F64, J::F64, N::F64)
    sorry()
end

"""
    cal_dc_exact(U::F64, J::F64, N::F64)

Evaluate the double counting term by the exact scheme.

See also: [`cal_dc_fll`](@ref), [`cal_dc_amf`](@ref).
"""
function cal_dc_exact(U::F64, J::F64, N::F64)
    sorry()
end

"""
    split_hyb_l(ai::Array{Impurity,1})

Split the dmft1/dmft.hyb_l file into impurity.[i]/dmft.hyb_l files, which
`i` denotes the index of quantum impurity problems.

See also: [`split_eimpx`](@ref), [`sigma_split`](@ref).
"""
function split_hyb_l(ai::Array{Impurity,1})
    # Declare the frequency mesh and hybridization function
    fmesh = []
    Delta = []

    # Filename for hybridization functions
    fhyb = "dmft1/dmft.hyb_l"

    # Make sure the existence of hybridization functions
    @assert isfile(fhyb)

    # Parse `fhyb`, extract the hybridization functions
    open(fhyb, "r") do fin
        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        nmesh = parse(I64, line_to_array(fin)[3])
        qdim  = parse(I64, line_to_array(fin)[4])
        @assert nsite == length(ai)

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for frequency mesh
        fmesh = zeros(F64, nmesh)

        # Create an array for hybridization functions
        Delta = zeros(C64, qdim, qdim, nmesh, nspin, nsite)

        # Read the data
        # Go through each impurity site and spin
        for t = 1:nsite
            for s = 1:nspin
                # Parse indices and dimensional parameter
                strs = readline(fin)
                _t = parse(I64, line_to_array(strs)[3])
                _s = parse(I64, line_to_array(strs)[5])
                _d = parse(I64, line_to_array(strs)[7])
                @assert _t == t
                @assert _s == s
                @assert _d == ai[t].nband

                # Go through each frequency point
                for m = 1:nmesh
                    # Parse frequency mesh
                    fmesh[m] = parse(F64, line_to_array(fin)[3])
                    # Parse hybridization functions
                    for q = 1:ai[t].nband
                        for p = 1:ai[t].nband
                            _re, _im = parse.(F64, line_to_array(fin)[3:4])
                            Delta[p,q,m,s,t] = _re + _im * im
                        end
                    end
                end # END OF M LOOP

                # Skip two lines
                readline(fin)
                readline(fin)
            end # END OF S LOOP
        end # END OF T LOOP
    end
    println("  Read hybridization functions from: $fhyb")

    # Next, we are going to split the hybridization functions according
    # to the quantum impurity problems

    # Extract the dimensional parameters
    _, qdim, nmesh, nspin, nsite = size(Delta)

    # Go through each quantum impurity problem
    for t = 1:nsite
        # Determine filename for hybridization functions
        fhyb = "impurity.$t/dmft.hyb_l"

        # Write the data
        open(fhyb, "w") do fout
            # Write dimensional parameters
            @printf(fout, "# nsite: %4i\n", nsite)
            @printf(fout, "# nspin: %4i\n", nspin)
            @printf(fout, "# nmesh: %4i\n", nmesh)
            @printf(fout, "# qdim : %4i\n", qdim)

            # Write separators
            println(fout)
            println(fout)

            # Go through each spin
            for s = 1:nspin
                # Write key parameters
                @printf(fout, "# site:%4i  spin:%4i  dims:%4i\n", t, s, ai[t].nband)

                # Go through each frequency point
                for m = 1:nmesh
                    @printf(fout, "w:%6i%16.8f\n", m, fmesh[m])
                    # Go through the orbital space
                    for q = 1:ai[t].nband
                        for p = 1:ai[t].nband
                            z = Delta[p,q,m,s,t]
                            @printf(fout, "%4i%4i%16.8f%16.8f\n", p, q, real(z), imag(z))
                        end
                    end
                end # END OF M LOOP

                # Write separators
                println(fout)
                println(fout)
            end # END OF S LOOP
        end
        println("  Split hybridization functions into: $fhyb")
    end # END OF T LOOP
end

"""
    split_eimpx(ai::Array{Impurity,1})

Split the dmft1/dmft.eimpx file into impurity.[i]/dmft.eimpx files, which
`i` denotes the index of quantum impurity problems.

See also: [`split_hyb_l`](@ref), [`sigma_split`](@ref).
"""
function split_eimpx(ai::Array{Impurity,1})
    # Declare the local impurity levels
    Eimpx = []

    # Filename for local impurity levels
    flev = "dmft1/dmft.eimpx"

    # Make sure the existence of local impurity levels
    @assert isfile(flev)

    # Parse `flev`, extract the local impurity levels
    open(flev, "r") do fin
        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        qdim  = parse(I64, line_to_array(fin)[4])
        @assert nsite == length(ai)

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for local impurity levels
        Eimpx = zeros(C64, qdim, qdim, nspin, nsite)

        # Read the data
        # Go through each impurity site and spin
        for t = 1:nsite
            for s = 1:nspin
                # Parse indices and dimensional parameter
                strs = readline(fin)
                _t = parse(I64, line_to_array(strs)[3])
                _s = parse(I64, line_to_array(strs)[5])
                _d = parse(I64, line_to_array(strs)[7])
                @assert _t == t
                @assert _s == s
                @assert _d == ai[t].nband

                # Parse local impurity levels
                for q = 1:ai[t].nband
                    for p = 1:ai[t].nband
                        _re, _im = parse.(F64, line_to_array(fin)[3:4])
                        Eimpx[p,q,s,t] = _re + _im * im
                    end
                end

                # Skip two lines
                readline(fin)
                readline(fin)
            end # END OF S LOOP
        end # END OF T LOOP
    end
    println("  Read local impurity levels from: $flev")

    # Next, we are going to split the local impurity levels according
    # to the quantum impurity problems

    # Extract the dimensional parameters
    _, qdim, nspin, nsite = size(Eimpx)

    # Go through each quantum impurity problems
    for t = 1:nsite
        # Determine filename for local impurity levels
        flev = "impurity.$t/dmft.eimpx"

        # Write the data
        open(flev, "w") do fout
            # Write dimensional parameters
            @printf(fout, "# nsite: %4i\n", nsite)
            @printf(fout, "# nspin: %4i\n", nspin)
            @printf(fout, "# qdim : %4i\n", qdim)

            # Write separators
            println(fout)
            println(fout)

            # Go through each spin
            for s = 1:nspin
                # Write key parameters
                @printf(fout, "# site:%4i  spin:%4i  dims:%4i\n", t, s, ai[t].nband)

                # Go through the orbital space
                for q = 1:ai[t].nband
                    for p = 1:ai[t].nband
                        z = Eimpx[p,q,s,t]
                        @printf(fout, "%4i%4i%16.8f%16.8f\n", p, q, real(z), imag(z))
                    end
                end

                # Write separators
                println(fout)
                println(fout)
            end # END OF S LOOP
        end
        println("  Split local impurity levels into: $flev")
    end # END OF T LOOP
end
