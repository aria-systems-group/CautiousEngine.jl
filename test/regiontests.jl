using CautiousEngine

function test_intersect_2d() 
    shape1 = Dict("vertices" => [[-1. -1.], [-1. 1.], [1. 1.], [1. -1.]])
    shape2 = Dict("vertices" => [[-.5 -1.], [-.5 1.], [-4. 1.], [-4. -1.]])
    axes = [1. 0.;0. 1.]                                 
    res = CautiousEngine.do_shapes_intersect(shape1, shape2, axes=axes)
end

function test_intersect()
    shape1 = Dict("vertices" => [[-1. -1. -1.], [-1. 1. -1.], [1. 1. -1.], [1. -1. -1.],
                                 [-1. -1. 1.], [-1. 1. 1.], [1. 1. 1.], [1. -1. 1.]])
    shape2 = Dict("vertices" => [[-5. -2. -2.], [-5. 2. -2.], [-4. 2. -2.], [-4. -2. -2.],
                                 [-5. -2. 2.], [-5. 2. 2.], [-4. 2. 2.], [-4. -2. 2.]])
    axes = [1. 0. 0.;0. 1. 0.; 0. 0. 1]                                 
    res = CautiousEngine.do_shapes_intersect(shape1, shape2, axes=axes)
return res
end

function test_containment()

    shape1 = Dict("vertices" => [[-1. -1. -1.], [-1. 1. -1.], [1. 1. -1.], [1. -1. -1.],
                                 [-1. -1. 1.], [-1. 1. 1.], [1. 1. 1.], [1. -1. 1.]])
    shape2 = Dict("vertices" => [[-2. -2. -2.], [-2. 2. -2.], [2. 2. -2.], [2. -2. -2.],
                                 [-2. -2. 2.], [-2. 2. 2.], [2. 2. 2.], [2. -2. 2.]])
    # axes = [1. 0. 0.;0. 1. 0.; 0. 0. 1]                                 
    res = CautiousEngine.does_shape2_contain_shape1(shape1, shape2, axes=nothing)
    return res
end

@testset "region intersection tests" begin
    @test test_intersect_2d()[1]
	@test !test_intersect()[1]
    @test test_containment()[1]
end