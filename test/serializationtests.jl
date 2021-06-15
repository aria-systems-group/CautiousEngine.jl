using CautiousEngine
using BSON

@testset "serialization and deserialization test" begin
    N = 100 
    x = rand(2, N)
    y = [sum(x[:, i]) for i=1:N]
    test_gp = CautiousEngine.train_gp_1dim(x, y)

	region_dict = Dict(1=>Dict("x1" => [-0.1, 0.], "x2" => [-0.1, 0.]),
				   2=>Dict("x1" => [0., 0.1], "x2" => [0., 0.1]))      # Defines the compact safe set

	region_post_dict = Dict(1 => (Dict("x1" => [0.1, 0.2], "x2" => [0.3, 0.4]), [0.5, 0.6]),
							2 => (Dict("x1" => [0.1, 0.2], "x2" => [0.3, 0.4]), [0.5, 0.6]))

	region_array = CautiousEngine.serialize_region_data(region_dict, region_post_dict)
	test_filename = "test.bson"
	bson(test_filename, Dict(:res => region_array))
	r_region_array = BSON.load(test_filename)
	@test r_region_array[:res] === region_array
end