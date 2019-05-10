# push!(LOAD_PATH,"../src/")
using Documenter
include("../src/simulate.jl")

makedocs(sitename="ISRL-MBD Documentation",
         doctest = false,
         format = :html,
         pages = ["Introduction" => "index.md",
                  "User Documentation" => [
                  "Simulate" => "simulate.md"
                  "Rigid Body" => "rigidbody.md"
                  "Joint" => "joint.md"
                  "Constraint Forces" => "constraintUK.md"
                  ]])


# using Documenter, ForwardDiff
#
# makedocs(modules=[ForwardDiff],
#          doctest = false,
#          format = :html,
#          sitename = "ForwardDiff",
#          pages = ["Introduction" => "index.md",
#                   "User Documentation" => [
#                     "Limitations of ForwardDiff" => "user/limitations.md",
#                     "Differentiation API" => "user/api.md",
#                     "Advanced Usage Guide" => "user/advanced.md",
#                     "Upgrading from Older Versions" => "user/upgrade.md"],
#                   "Developer Documentation" => [
#                     "How ForwardDiff Works" => "dev/how_it_works.md",
#                     "How to Contribute" => "dev/contributing.md"]])
#
# deploydocs(repo = "github.com/JuliaDiff/ForwardDiff.jl.git",
#            osname = "linux",
#            julia = "0.7",
#            target = "build",
#            deps = nothing,
#            make = nothing)
