using CautiousEngine

@testset "adding and chanding GP data" begin
    N = 3
    x = rand(2, N)
    y = [sum(x[:, i]) for i=1:N]
    @test size(x, 2) == length(y)
    test_gp = CautiousEngine.train_gp_1dim(x, y)
    @test x == test_gp.x
    @test y == test_gp.y

    new_x = rand(2,1)
    new_y = sum(new_x)
    new_datapoint = (new_x, new_y)

    CautiousEngine.add_datapoint_to_gp(test_gp, new_datapoint)
    @test test_gp.nobs == N+1

    # Test Updating GP
    M = 5
    x_replace = rand(2,M)
    @test size(x_replace, 2) == M
    y_replace = [sum(x_replace[:, i]) for i=1:M]
    updated_gp = CautiousEngine.update_gp_data(test_gp, x_replace, y_replace)
    @test updated_gp.x === x_replace
    @test updated_gp.y === y_replace
    @test updated_gp.nobs == M

    # Just check if dimension is correct
    @test updated_gp.dim == 2

    # KNN Test
    P = 10
    x = rand(2,P)
    @info typeof(x)
    y = rand(1,P)
    center = [0.;0.]
    NN = 3
    x_sub, y_sub = CautiousEngine.get_local_data_knn(center, x, y; num_neighbors=NN)
    @test size(x_sub, 2) == NN
    @test length(y_sub) == NN
end