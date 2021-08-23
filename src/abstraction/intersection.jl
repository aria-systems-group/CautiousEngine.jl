" Implementation of a separating axis test for the intersection of two convex shapes. 
# Arguments
- `shape1::Dict` - Convex shape with entry `vertices`
- `shape2::Dict` - Convex shape with entry `vertices`
"
function do_shapes_intersect(shape1::Dict, shape2::Dict; axes=nothing)
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

" Checks if shape2 contains shape1. 
# Arguments
- `shape1::Dict` - Convex shape with entry `vertices`
- `shape2::Dict` - Convex shape with entry `vertices`
"
function does_shape2_contain_shape1(shape1::Dict, shape2::Dict; axes=nothing)
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

" Projects the shape vertices onto the axis. 
# Arguments
- `shape::Dict` - Shape with entry `vertices`
- `axis` - Axis of projection
"
function project(shape::Dict, axis)
    dots = [(vertex*axis)[1] for vertex in shape["vertices"]]
    min = minimum(dots)
    max = maximum(dots)
    return [min, max]
end

" Gets the principal normals of each shape.
# Arguments
- `shapes::Vector{Dict{Any, Any}}` - Array with shape dictionaries
"
function get_axes(shapes::Vector{Dict{Any, Any}})
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
