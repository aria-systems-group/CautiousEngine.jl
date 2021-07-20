function save_region_data(region_dict, region_post_dict, filename; lazy_flag=false)
	region_data = serialize_region_data(region_dict, region_post_dict, lazy_flag=lazy_flag)
	bson(filename, Dict(:res => region_data))
end

function serialize_region_data(region_dict, region_post_dict; lazy_flag=false)

	if lazy_flag
		return serialize_region_data_lazy(region_dict, region_post_dict)
	end

	num_regions = length(keys(region_dict))
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

function serialize_region_data_lazy(region_dict, region_post_dict)
	num_regions = length(keys(region_dict))
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

function deserialize_region_data(region_data_filename; lazy_flag=false)

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

function deserialize_region_data_lazy(region_data_filename)
	region_data_array = BSON.load(region_data_filename)[:res]
	num_regions = size(region_data_array,1)
	pair_iterator = Base.product(1:num_regions, 1:num_regions)
	region_dict = Dict()
	region_post_dict = Dict()
	
	for i=1:num_regions-1
		region_dict[i] = Dict()
		region_dict[i]["x1"] = region_data_array[i,1:2]
		region_dict[i]["x2"] = region_data_array[i,3:4]
		
		# Ignoring post map for now
	end

	region_dict[num_regions] = Dict()
	region_dict[num_regions]["x1"] = region_data_array[num_regions,1:2]
	region_dict[num_regions]["x2"] = region_data_array[num_regions,3:4]
	
	return region_dict, region_post_dict, pair_iterator
end