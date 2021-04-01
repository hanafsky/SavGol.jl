using SavGol
using Test
using DelimitedFiles

@testset "SavGol.jl" begin
    # Write your own tests here.
    data  = readdlm("signal.txt")
    sg = SG(20,3)
    smoothed = apply_filter(sg[:,1],data)

end
