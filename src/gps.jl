# Methods related to all things GP regression

"""
Update the GP with the datapoint given as an input-output tuple. 
"""
function update_gp(gp, datapoint::Tuple)
    # Update the GP in place - does not return a new GP structure. 
    newx = [gp.x'; datapoint[1][:]']'
    newy = gp.y
    push!(newy, datapoint[2])
    GaussianProcesses.fit!(gp, newx, newy)
end

"""
Load a GP bin file.
"""
function load_gp(filename)
    f = open(filename)
    gp_set = deserialize(f)
    close(f)
    return gp_set
end
