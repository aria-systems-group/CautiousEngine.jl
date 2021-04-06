using CautiousEngine

@testset "system composition tests" begin

    minPr_mats = [[0.0 0.1; 0.2 0.3], [0.4 0.5; 0.6 0.7]]
    maxPr_mats = [[0.1 0.2; 0.3 0.4], [0.5 0.6; 0.7 0.8]]

    exp_minPr_mat = [0.0 0.1; 0.4 0.5; 0.2 0.3; 0.6 0.7]
    exp_maxPr_mat = [0.1 0.2; 0.5 0.6; 0.3 0.4; 0.7 0.8]

    minPr_mat, maxPr_mat = CautiousEngine.create_switched_system_matrix(minPr_mats, maxPr_mats)

    @test minPr_mat == exp_minPr_mat
    @test maxPr_mat == exp_maxPr_mat
end