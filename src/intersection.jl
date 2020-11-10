"""
Implementation of a Separating Axis Test
Input: Two Convex Shapes as a Dict
Keyword axes allows us to calculate all axes for regions ahead of time.
"""

function do_shapes_intersect(shape1, shape2; axes=nothing)
    int_result = true
    min_distances = [-1.] 

    if isnothing(axes)
        axes = get_axes([shape1, shape2])
    end

    num_axes = size(axes)[2]
    for i=1:num_axes
        axis = axes[:, i:i]
        p1 = project(shape1, axis)
        p2 = project(shape2, axis)
        overlap = (p2[1] <= p1[1] <= p2[2] || p2[1] <= p1[2] <= p2[2] ||
                   p1[1] <= p2[1] <= p1[2] || p1[1] <= p2[2] <= p1[2] )
        push!(min_distances, minimum([abs(p1[2] - p2[1]), abs(p2[2] - p1[1])]))
        if !overlap
            # push!(min_distances, minimum([abs(p1[2] - p2[1]), abs(p2[2] - p1[1])]))
            int_result = false 
        end
    end

    return int_result, maximum(min_distances) 
end

function does_shape2_contain_shape1(shape1, shape2; axes=nothing)
    int_result = false 
    min_distance = Inf

    if isnothing(axes)
        axes = get_axes([shape1, shape2])
    end

    num_axes = size(axes)[2]
    for i=1:num_axes
        axis = axes[:, i:i]
        p1 = project(shape1, axis)
        p2 = project(shape2, axis)
        # overlap = (p2[1] <= p1[1] <= p2[2] || p2[1] <= p1[2] <= p2[2] ||
        #            p1[1] <= p2[1] <= p1[2] || p1[1] <= p2[2] <= p1[2] )
        contain = p2[1] <= p1[1] <= p2[2] && p2[1] <= p1[2] <= p2[2]
        if contain
            min_distance = minimum([abs(p1[1] - p2[1]), abs(p1[2] - p2[2]), min_distance])
            int_result = true 
        else
            int_result = false
            break
        end
    end

    return int_result, min_distance
end

"""
Function that finds the projection on the axis.
"""

function project(shape, axis)
    dots = [(vertex*axis)[1] for vertex in shape["vertices"]]
    min = minimum(dots)
    max = maximum(dots)
    return [min, max]
end

function get_axes(shapes)
    axes = []
    for shape in shapes
        num_v = length(shape["vertices"])
        for i=1:num_v
            normals = nullspace(shape["vertices"][i]-(i==num_v ? shape["vertices"][1] : shape["vertices"][i+1])) 
            if isempty(axes)
                axes = normals
            else
                # TODO: remove redundant axes
                axes = hcat([axes normals])
            end 
        end
    end
    return axes
end

function test_time()

    N = 10000

    @time for i=1:N

        # generate a random 3D object
        shape1 = Dict("vertices" => [10*rand(1,3), 10*rand(1,3), 10*rand(1,3)])
        shape2 = Dict("vertices" => [10*rand(1,3), 10*rand(1,3), 10*rand(1,3)])
        do_shapes_intersect(shape1, shape2) 
    end

end

function test_intersect_2d() 
    shape1 = Dict("vertices" => [[-1. -1.], [-1. 1.], [1. 1.], [1. -1.]])
    shape2 = Dict("vertices" => [[-.5 -1.], [-.5 1.], [-4. 1.], [-4. -1.]])
    axes = [1. 0.;0. 1.]                                 
    res = do_shapes_intersect(shape1, shape2, axes=axes)
end

function test_intersect()
    shape1 = Dict("vertices" => [[-1. -1. -1.], [-1. 1. -1.], [1. 1. -1.], [1. -1. -1.],
                                 [-1. -1. 1.], [-1. 1. 1.], [1. 1. 1.], [1. -1. 1.]])
    shape2 = Dict("vertices" => [[-5. -2. -2.], [-5. 2. -2.], [-4. 2. -2.], [-4. -2. -2.],
                                 [-5. -2. 2.], [-5. 2. 2.], [-4. 2. 2.], [-4. -2. 2.]])
    axes = [1. 0. 0.;0. 1. 0.; 0. 0. 1]                                 
    res = do_shapes_intersect(shape1, shape2, axes=axes)
return res
end

function test_containment()

    shape1 = Dict("vertices" => [[-1. -1. -1.], [-1. 1. -1.], [1. 1. -1.], [1. -1. -1.],
                                 [-1. -1. 1.], [-1. 1. 1.], [1. 1. 1.], [1. -1. 1.]])
    shape2 = Dict("vertices" => [[-2. -2. -2.], [-2. 2. -2.], [2. 2. -2.], [2. -2. -2.],
                                 [-2. -2. 2.], [-2. 2. 2.], [2. 2. 2.], [2. -2. 2.]])
    # axes = [1. 0. 0.;0. 1. 0.; 0. 0. 1]                                 
    res = does_shape2_contain_shape1(shape1, shape2, axes=nothing)
    return res
end