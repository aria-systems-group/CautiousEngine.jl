using CautiousEngine
using Test

@testset "CautiousEngine.jl" begin
    # Write your tests here.
    include("gptests.jl")
    include("dfatests.jl")
    include("serializationtests.jl")
    include("systemtests.jl")
    include("regiontests.jl")
end
