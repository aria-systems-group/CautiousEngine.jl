" Saves region data.
# Arguments
- `region_dict::Dict` - Dictionary with region extents 
- `region_post_dict::Dict` - Dictionary with region image overapproximations
- `filename::String` - Location to save region data
"
function save_region_data(region_dict::Dict, region_post_dict::Dict, filename::String; lazy_flag=false)
	region_data = serialize_region_data(region_dict, region_post_dict, lazy_flag=lazy_flag)
	bson(filename, Dict(:res => region_data))
end

" Saves region data to the experiment root directory.
# Arguments
- `exp_dir::String` - Experiment root directory
- `region_data::Dict` - Dictionary containing region extent, post dictionaries 
"
function save_region_data(exp_dir::String, region_data::Dict)
    @info "Saving the region data to the experiment directory..."
    region_filename = "$exp_dir/regions.bson"
    save_region_data(region_data[:extents], region_data[:posts], region_filename)
end

" Loads region data.
# Arguments 
- `filename::String` - Region data file to load
"
function load_region_data(filename::String) 
	region_dict, region_post_dict, _ = deserialize_region_data(filename)
    region_data = Dict(:extents => region_dict, :posts => region_post_dict)
	return region_data
end

function create_region_data(space::Dict, grid_sizes::Dict)
    region_dict = discretize_set(space, grid_sizes)
    # ? THIS IS NOT CORRECT for a 3d model ? Not sure what I mean by this comment 
    region_dict[region_dict.count+1] = space # Region with index L+1 corresponds to the safe set as a whole. 
    return region_dict
end

function create_region_data_new(space::Dict, grid_size::Dict)
    region_dict, extent_dict = discretize_set_lazy(space, grid_size)
    x = space["x1"]
    y = space["x2"]
    region_dict[region_dict.count+1] = VPolygon([[x[1], y[1]], [x[1], y[2]], 
                                    [x[2], y[2]], [x[2], y[1]]])
    extent_dict[region_dict.count+1] = space

    return region_dict, extent_dict
end

function serialize_region_data(region_dict::Dict, region_post_dict::Dict; lazy_flag=false)

	if lazy_flag
		return serialize_region_data_lazy(region_dict, region_post_dict)
	end

	num_regions = region_dict.count
	# Assuming 2D square for now:
	M = 10
	region_data_array = zeros(num_regions, M)
	# Format: [minX, maxX, minY, maxY, μlX, μUX, μlY, μUY, σUX, σUY]
	for i=1:num_regions - 1
		region_data_array[i,1:2] = region_dict[i]["x1"]
		region_data_array[i,3:4] = region_dict[i]["x2"]
		region_data_array[i,5:6] = region_post_dict[i][1]["x1"]
		region_data_array[i,7:8] = region_post_dict[i][1]["x2"]
		region_data_array[i,9:10] = region_post_dict[i][2][:]
	end

	region_data_array[num_regions,1:2] = region_dict[num_regions]["x1"]
	region_data_array[num_regions,3:4] = region_dict[num_regions]["x2"]

	return region_data_array
end

function serialize_region_data_lazy(region_dict::Dict, region_post_dict)
	num_regions = region_dict.count
	# Assuming 2D square for now:
	M = 5
	region_data_array = zeros(num_regions, M)
	# Format: [minX, maxX, minY, maxY, postmap]
	for i=1:num_regions - 1
		region_data_array[i,1:2] = region_dict[i]["x1"]
		region_data_array[i,3:4] = region_dict[i]["x2"]
		region_data_array[i,5] = 0.
	end

	region_data_array[num_regions,1:2] = region_dict[num_regions]["x1"]
	region_data_array[num_regions,3:4] = region_dict[num_regions]["x2"]

	return region_data_array	
end

function deserialize_region_data(region_data_filename::String; lazy_flag=false)

	if lazy_flag
		return deserialize_region_data_lazy(region_data_filename)
	end

	region_data_array = BSON.load(region_data_filename)[:res]
	num_regions = size(region_data_array,1)
	pair_iterator = Base.product(1:num_regions, 1:num_regions)
	region_dict = Dict()
	region_post_dict = Dict()
	
	for i=1:num_regions-1
		region_dict[i] = Dict()
		region_dict[i]["x1"] = region_data_array[i,1:2]
		region_dict[i]["x2"] = region_data_array[i,3:4]
		
		region_post_dict[i] = (Dict(), region_data_array[i,9:10])
		region_post_dict[i][1]["x1"] = region_data_array[i,5:6]
		region_post_dict[i][1]["x2"] = region_data_array[i,7:8]
	end

	region_dict[num_regions] = Dict()
	region_dict[num_regions]["x1"] = region_data_array[num_regions,1:2]
	region_dict[num_regions]["x2"] = region_data_array[num_regions,3:4]
	
	return region_dict, region_post_dict, pair_iterator
end

function deserialize_region_data_lazy(region_data_filename::String)
	region_data_array = BSON.load(region_data_filename)[:res]
	num_regions = size(region_data_array,1)
	pair_iterator = Base.product(1:num_regions, 1:num_regions)
	region_dict = Dict()
	region_post_dict = Dict()
	
	for i=1:num_regions-1
		region_dict[i] = Dict()
		region_dict[i]["x1"] = region_data_array[i,1:2]
		region_dict[i]["x2"] = region_data_array[i,3:4]
		
		# TODO: Ignoring post map for now
	end

	region_dict[num_regions] = Dict()
	region_dict[num_regions]["x1"] = region_data_array[num_regions,1:2]
	region_dict[num_regions]["x2"] = region_data_array[num_regions,3:4]
	
	return region_dict, region_post_dict, pair_iterator
end