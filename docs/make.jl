haskey(ENV,"ZEN_CORE") && pushfirst!(LOAD_PATH, ENV["ZEN_CORE"])

using Documenter
using ZenCore

makedocs(
    sitename = "ZenCore",
    clean = false,
    authors = "Li Huang <huangli@caep.cn> and contributors",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
        repolink = "https://github.com/huangli712/ZenCore",
        size_threshold = 409600, # 400kb
        assets = ["assets/zencore.css"],
        collapselevel = 1,
    ),
    remotes = nothing,
    modules = [ZenCore],
    pages = [
        "Home" => "index.md",
        "Introduction" => "intro.md",
        "Installation" => "install.md",
        "Inside The Library" => Any[
            "Outline" => "library/outline.md",
            "ZenCore" => "library/zencore.md",
            "Global Constants" => "library/global.md",
            "Utilities" => "library/util.md",
            "Tetrahedron Integration" => "library/tetra.md",
            "Types" => "library/types.md",
            "Configuration" => "library/config.md",
            "Workflow" => "library/base.md",
            "Vienna Ab initio Simulation Package" => "library/vasp.md",
            "Quantum ESPRESSO" => "library/qe.md",
            "Projected Local Orbital" => "library/plo.md",
            "Wannier Function" => "library/wannier.md",
            "Intermediate Representation" => "library/ir.md",
            "Dynamical Mean-Field Theory" => "library/dmft.md",
            "Quantum Impurity Solver" => "library/solver.md",
            "Self-Energy Function" => "library/sigma.md",
            "Mixing" => "library/mixing.md",
        ],
    ],
)
