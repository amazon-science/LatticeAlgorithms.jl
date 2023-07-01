push!(LOAD_PATH,"../src/")

using Documenter, LatticeAlgorithms 

makedocs(sitename="LatticeAlgorithms.jl")

deploydocs(repo="https://github.com/amazon-science/LatticeAlgorithms.jl",)