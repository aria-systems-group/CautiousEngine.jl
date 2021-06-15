function serialize_region_data(region_dict, region_post_dict)

	num_regions = length(keys(region_dict))
	# region_data[:pairs] = region_pairs	# This does not need to be serialized.
	# region_pairs = product(1..num_regions, 1..num_regions)
    region_data[:extents] = region_dict		# doesn't really need to be serialized...
    region_data[:posts] = region_post_dict
    # region_data[:gps] = region_gp_dic

	# Assuming 2D square for now:

	M = 10
	region_data_array = zeros(num_regions, M)
	
	# Format: [minX, maxX, minY, maxY, μlX, μUX, μlY, μUY, σUX, σUY]
	for i=1:keys(region_dict)
		region_data_array[i,1:2] = region_dict[i]["x1"]
		region_data_array[i,3:4] = region_dict[i]["x2"]
		region_data_array[i,5:6] = region_post_dict[i][1]["x1"]
		region_data_array[i,7:8] = region_post_dict[i][1]["x2"]
		region_data_array[i,9:10] = region_post_dict[i][2][:]
	end
	
	return region_data_array
end