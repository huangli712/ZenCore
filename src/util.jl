#
# Project : Pansy
# Source  : util.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/25
#

#=
### *Basic Macros*
=#

"""
    @cswitch(constexpr, body)

Provides a C-like switch statement with the *falling through* behavior.
This implementation was borrowed from the following github repository:
* https://github.com/Gnimuc/CSyntax.jl

### Examples
```julia
engine = get_d("engine")
@cswitch engine begin
    @case "vasp"
        just_do_it()
        break

    @default
        sorry()
        break
end
```
"""
macro cswitch(constexpr, body)
    case2label = Dict{Any,Symbol}()
    flow = Expr(:block)
    end_label = gensym("end")
    default_label = end_label

    for arg in body.args
        if Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@case")
            label = gensym("case")
            case2label[arg.args[3]] = label
            labelexpr = Expr(:symboliclabel, label)
            push!(flow.args, labelexpr)
        elseif Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@default")
            default_label = gensym("default")
            labelexpr = Expr(:symboliclabel, default_label)
            push!(flow.args, labelexpr)
        elseif arg == Expr(:break)
            labelexpr = Expr(:symbolicgoto, end_label)
            push!(flow.args, labelexpr)
        else
            push!(flow.args, arg)
        end
    end
    push!(flow.args, Expr(:symboliclabel, end_label))

    jumptable = Expr(:block)
    for (case, label) in case2label
        condition = Expr(:call, :(==), constexpr, case)
        push!(jumptable.args, Expr(:if, condition, Expr(:symbolicgoto, label)))
    end
    push!(jumptable.args[end].args, Expr(:symbolicgoto, default_label))

    return esc(Expr(:block, jumptable, flow))
end

"""
    @ps1(str, c)

Wrapper for `printstyled` function. Here `str` is a string, and `c`
denotes color.

### Examples
```julia
@ps1 "Hello world!" :green
```

See also: [`@ps2`](@ref).
"""
macro ps1(str, c)
    return :( printstyled($str, color = $c) )
end

"""
    @ps2(str1, c1, str2, c2)

Wrapper for `printstyled` function. Here `str1` and `str2` are strings,
and `c1` and `c2` denote colors.

See also: [`@ps1`](@ref).
"""
macro ps2(str1, c1, str2, c2)
    ex = quote
        printstyled($str1, color = $c1)
        printstyled($str2, color = $c2)
    end
    return :( $(esc(ex)) )
end

#=
### *Query Runtime Environment*
=#

"""
    require()

Check the version of julia runtime environment. It should be higher
than v1.6.x. One of the most important philosophies of the `ZenCore`
package is minimizing the dependence on the third-party libraries as
far as possible. Note that the `ZenCore` package relys on the `TOML`
package to parse the *.toml file. Only in v1.6.0 and higher versions,
julia includes the `TOML` package in its standard library.
"""
function require()
    if VERSION < v"1.6-"
        error("Please upgrade your julia to v1.6.0 or higher")
    end
end

"""
    setup_args(x::Vararg{String})

Setup `ARGS` manually. This function is used only in `REPL` environment.
We can use this function to update `ARGS`, so that the `query_args()`
and the other related functions can work correctly.

### Examples
```julia-repl
julia > setup_args("SrVO3.toml")
1-element Array{String,1}:
 "SrVO3.toml"
```

See also: [`query_args`](@ref).
"""
function setup_args(x::Vararg{String})
    # Make sure it is the REPL
    @assert isinteractive()

    # Clean `ARGS`
    empty!(ARGS)

    # Convert the arguments to an array of strings
    X = collect(x)

    # Push them into `ARGS` one by one
    for i in eachindex(X)
        push!(ARGS, X[i])
    end

    # Return `ARGS`, only for debug.
    ARGS
end

"""
    query_args()

Check whether the configuration file (`case.toml`) is provided.

See also: [`setup_args`](@ref).
"""
function query_args()
    nargs = length(ARGS)
    if nargs < 1
        error("Please specify the configuration file")
    else
        ARGS[1]
    end
end

"""
    query_case()

Return `case`, in other words, the job's name.

See also: [`query_stop`](@ref).
"""
function query_case()
    basename( splitext(query_args())[1] )
end

#=
*Remarks*:

For `VASP`, the essential input files include:
* POSCAR
* POTCAR

As for the `INCAR` and `KPOINTS`, they will be generated automatically.
=#

"""
    query_inps(engine::String)

Check whether the essential input files exist. This function is designed
for the DFT engine only. The input files for the DMFT engine, quantum
impurity solver, and Kohn-Sham adaptor will be generated automatically.
The `ZenCore` package will take care of them. Do not worry.

See also: [`query_case`](@ref).
"""
function query_inps(engine::String)
    @cswitch engine begin
        @case "vasp"
            if !isfile("POSCAR") || !isfile("POTCAR")
                error("Please provide both POSCAR and POTCAR files")
            end
            break

        @default
            sorry()
            break
    end
end

"""
    query_stop()

Query whether the `case.stop` file exists.

See also: [`query_case`](@ref).
"""
function query_stop()
    isfile(query_case() * ".stop")
end

#=
*Remarks*:

In the `ZenCore` package, the following environment variables matter:

* ZEN_HOME
* ZEN_CORE
* ZEN_DMFT
* ZEN_SOLVER
* VASP_HOME

Please setup them in your `.bashrc` (Lniux) or `.profile` (macOS) files.
=#

"""
    query_home()

Query the home directory of `Zen`. Actually, the `ZEN_HOME` means the
directory that the `Zen` Framework is installed.

See also: [`query_core`](@ref).
"""
function query_home()
    # We have to setup the environment variable ZEN_HOME
    if haskey(ENV, "ZEN_HOME")
        ENV["ZEN_HOME"]
    else
        error("ZEN_HOME is undefined")
    end
end

"""
    query_core()

Query the src/core directory of `Zen`. Actually, the `ZEN_CORE` denotes
the directory that contains the ZenCore.jl file. Be careful, `ZEN_CORE`
must be included in `LOAD_PATH`.

See also: [`query_home`](@ref).
"""
function query_core()
    # We have to setup the environment variable ZEN_CORE
    if haskey(ENV, "ZEN_CORE")
        ENV["ZEN_CORE"]
    else
        error("ZEN_CORE is undefined")
    end
end

"""
    query_dft(engine::String)

Query the home directory of the chosen DFT engine.

See also: [`query_dmft`](@ref), [`query_solver`](@ref).
"""
function query_dft(engine::String)
    @cswitch engine begin
        # We have to setup the environment variable VASP_HOME
        @case "vasp"
            if haskey(ENV, "VASP_HOME")
                return ENV["VASP_HOME"]
            else
                error("VASP_HOME is undefined")
            end
            break

        @default
            sorry()
            break
    end
end

"""
    query_dmft()

Query the home directory of the DMFT engine.

See also: [`query_dft`](@ref), [`query_solver`](@ref).
"""
function query_dmft()
    # We have to setup the environment variable ZEN_DMFT
    if haskey(ENV, "ZEN_DMFT")
        ENV["ZEN_DMFT"]
    # For develop stage only
    else
        joinpath(query_home(), "src/dmft")
    end
end

#=
*Remarks*:

The `atomic` code is considered as a preprocessor of `ct_hyb2`, it is
not a valid quantum impurity solver.
=#

"""
    query_solver(engine::String)

Query the home directory of the quantum impurity solver.

See also: [`query_dft`](@ref), [`query_dmft`](@ref).
"""
function query_solver(engine::String)
    # We have to setup the environment variable ZEN_SOLVER
    if haskey(ENV, "ZEN_SOLVER")
        ENV["ZEN_SOLVER"]
    # For develop stage only
    else
        @cswitch engine begin
            @case "ct_hyb1"
                solver_dir = joinpath(query_home(), "src/solver/ct_hyb1")
                break

            @case "ct_hyb2"
                solver_dir = joinpath(query_home(), "src/solver/ct_hyb2")
                break

            @case "hub1"
                solver_dir = joinpath(query_home(), "src/solver/hub1")
                break

            @case "norg"
                solver_dir = joinpath(query_home(), "src/solver/norg")
                break

            @case "atomic"
                solver_dir = joinpath(query_home(), "src/solver/atomic")
                break

            @default
                sorry()
                break
        end
        solver_dir
    end
end

#=
### *Colorful Outputs*
=#

"""
    welcome()

Print out the welcome messages to the screen.
"""
function welcome()
    @ps1 "                                        |\n" :green
    @ps2 "ZZZZZZZZZZZZ EEEEEEEEEEEE NNNNNNNNNNNN  | "  :green "A Modern DFT + DMFT Computation Framework\n" :magenta
    @ps1 "          Z               N          N  |\n" :green
    @ps1 "         Z                N          N  |\n" :green
    @ps2 "   ZZZZZZ    EEEEEEEEEEEE N          N  | "  :green "Package: $__LIBNAME__\n" :magenta
    @ps2 "  Z                       N          N  | "  :green "Version: $__VERSION__\n" :magenta
    @ps2 " Z                        N          N  | "  :green "Release: $__RELEASE__\n" :magenta
    @ps2 "ZZZZZZZZZZZZ EEEEEEEEEEEE N          N  | "  :green "Powered by the julia programming language\n" :magenta
    @ps1 "                                        |\n" :green
    println()
end

"""
    overview()

Print out the overview of Zen to the screen.
"""
function overview()
    # Build strings
    str1 = nprocs() === 1 ? " processor " : " processors "
    str2 = "(myid = $(myid()))"

    # Write the information
    prompt("ZEN", "Overview")
    println("Time : ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
    println("Para : Using ", nprocs(), str1, str2)
    println("Dirs : ", pwd())
    println("Task : ", query_args())
    println()
end

"""
    goodbye()

Print the goodbye messages to the screen.
"""
function goodbye()
    println("See you later")
end

"""
    sorry()

Print an error message to the screen.
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

"""
    prompt(from::String, msg::String)

Print a format Zen message to the screen.
"""
function prompt(from::String, msg::String)
    @ps2 "$from > " :green msg :magenta
    println()
end

"""
    prompt(msg::String)

Print a format Zen message to the screen.
"""
function prompt(msg::String)
    @ps2 "Task -> " :blue msg :light_red
    println()
end

"""
    prompt(io::IOStream, msg::String)

Print a format Zen message to the given IOStream. This function is used
to log the events during DFT + DMFT iterations.
"""
function prompt(io::IOStream, msg::String)
    date = Dates.format(now(), "yyyy-mm-dd / HH:MM:SS")
    println(io, "$date  $msg")
    flush(io)
end

#=
### *I/O Operations*
=#

"""
    line_to_array(io::IOStream)

Convert a line (reading from an IOStream) to a string array.
"""
@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

"""
    line_to_array(str::AbstractString)

Convert a string (AbstractString) to a string array.
"""
@inline function line_to_array(str::AbstractString)
    split(str, " ", keepempty = false)
end

"""
    line_to_cmplx(io::IOStream)

Convert a line (reading from an IOStream) to a cmplx number. It is used
to parse the `LOCPROJ` file only.

See also: [`vaspio_projs`](@ref).
"""
@inline function line_to_cmplx(io::IOStream)
    str = readline(io)
    _re = str[10:28]
    _im = str[30:end]
    return parse(F64, _re) + parse(F64, _im) * im
end

#=
### *Mathematical Functions*
=#

#=
*Remarks 1*:

The definition of Gauss error function is as follows:

```math
erf(x) = \frac{2}{\sqrt{\pi}}\int^{x}_{0} e^{-\eta^2} d\eta.
```

*Remarks 2*:

We call the `erf()` function defined in the mathematical library `libm`
or `openlibm` directly, instead of implementing it again by ourselves.
For more details about `libm` and `openlibm`, please visit the following
websites:
* https://openlibm.org
* https://github.com/JuliaMath/openlibm
* https://sourceware.org/newlib/libm.html

*Remarks 3*:

This below implementation is taken from the `SpecialFunctions.jl`. See:
* https://github.com/JuliaMath/SpecialFunctions.jl

*Remarks 4*:

`Base.Math.libm` is actually a string. It denotes `libopenlibm`.
=#

"""
    erf(x::F64)

Calculate the Gauss error function.

See also: [`gauss_weight`](@ref).
"""
function erf(x::F64)
    ccall(("erf", libm), F64, (F64,), x)
end

#=
### *Utility Functions*
=#

"""
    subscript(num::I64)

Convert a number (it must be in [0,9]) to subscript.
"""
function subscript(num::I64)
    @assert num >=0 && num <= 9
    SUB = ["\u2080" "\u2081" "\u2082" "\u2083" "\u2084"
           "\u2085" "\u2086" "\u2087" "\u2088" "\u2089"]
    return SUB[num + 1]
end

"""
    time_call(f)

Evaluate a function call (`f`), and then print the elapsed time (number
of seconds) it took to execute.
"""
function time_call(f)
    t = @elapsed f
    println("Elapsed time: $t s\n")
end
