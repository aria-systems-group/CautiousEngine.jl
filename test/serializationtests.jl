using CautiousEngine
using BSON

@testset "serialization and deserialization test" begin

	region_dict = Dict()
	region_post_dict = Dict()
	N_regions = 10000

	for i=1:N_regions-1
		region_dict[i] = Dict("x1" => [-0.1, 0.], "x2" => [-0.1, 0.])
		region_post_dict[i] = (Dict("x1" => [0.1, 0.2], "x2" => [0.3, 0.4]), [0.5, 0.6]) 
	end
	region_dict[N_regions] = Dict("x1" => [-0.1, 0.], "x2" => [-0.1, 0.])

	region_array = CautiousEngine.serialize_region_data(region_dict, region_post_dict)
	test_filename = "test.bson"
	bson(test_filename, Dict(:res => region_array))

	r_region_array = BSON.load(test_filename)
	@test r_region_array[:res] == region_array

	etime_new = @elapsed CautiousEngine.save_region_data(region_dict, region_post_dict, test_filename) 
	etime_old = @elapsed bson("foo.bson", Dict(:regions => region_dict, :region_posts => region_post_dict))

	println("new save time: ", etime_new)
	println("old save time: ", etime_old)

	r_region_dict, r_region_post_dict, _ = CautiousEngine.deserialize_region_data(test_filename)

	@test r_region_dict == region_dict
	@test r_region_post_dict == region_post_dict
end