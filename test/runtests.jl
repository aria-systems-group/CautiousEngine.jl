using CautiousEngine
using Test

@testset "CautiousEngine.jl" begin
    # Write your tests here.
    include("gptests.jl")
    include("dfatests.jl")
    include("systemtests.jl")
end
