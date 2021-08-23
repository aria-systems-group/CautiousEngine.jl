# TODO: Fix this function
function generate_linear_truth(params, system_matrix; single_mode_verification=false)
    logfile = initialize_log(params)
    timing_info = Dict()
    total_runtime = 0.
    # TODO: Move this
    @info "Determining the post-images of the regions under the linear map."
    region_time = @elapsed begin 
        region_dict, extents_dict = create_region_data_new(params.domain, params.discretization_step)
        region_post_dict = Dict()

        for i=1:region_dict.count
            region = region_dict[i]
            
            # Lazy representation saves time
            # TODO: Add the process noise term here for verification
            region_post = LinearMap(0.9999*system_matrix, region) 
            region_post_dict[i] = region_post 
        end 

        region_data = Dict()
        region_data[:extents] = extents_dict 
        region_data[:posts] = region_post_dict 
    end
    timing_info["region_bound_time_s"] = region_time 
    total_runtime += region_time
    @info "Region generation time: " region_time
    save_region_data(params, region_data)
    @info "Generating the transition bounds..."
    bound_time = @elapsed begin
        @info "Calculating transition probability bounds between regions..."
        results_df = DataFrame(Set1 = Int[], Set2 = Int[], MinPr = Float64[], MaxPr = Float64[], MinPrPoint = Array[], MaxPrPoint = Array[])
        for region_pair in Base.product(1:region_dict.count, 1:region_dict.count) 
            r1 = region_dict[region_pair[1]]
            r2 = region_dict[region_pair[2]]
            r1_post = region_post_dict[region_pair[1]]
        
            # Check lazy intersection
            cap = Intersection(r1_post, r2)
            if isempty(cap)
                prange = [0., 0.]
            else
                if isequivalent(r1_post, cap)
                    prange = [1., 1.]
                else
                    prange = [0., 1.]
                end
            end   

            df_row = [region_pair[1], region_pair[2], prange[1], prange[2], [-1.], [-1.]]
            push!(results_df, df_row)
        end

        # TODO: Move this constructor to its own fcn
        num_states = region_dict.count
        minPr_mat = spzeros(num_states, num_states)
        maxPr_mat = spzeros(num_states, num_states)
        for i = 1:1:num_states
            sub_row = results_df[results_df.Set1 .== i, :]
            if i == num_states
                minPr_mat[i, i] = 1.
                maxPr_mat[i, i] = 1.
            else
                for j = 1:1:num_states
                    if j == num_states
                        subsub_row = sub_row[sub_row.Set2 .== j, :]
                        minPr_mat[i, j] = 1. - subsub_row.MaxPr[1]
                        maxPr_mat[i, j] = 1. - subsub_row.MinPr[1]
                    else
                        subsub_row = sub_row[sub_row.Set2 .== j, :]
                        minPr_mat[i, j] = subsub_row.MinPr[1]
                        maxPr_mat[i, j] = subsub_row.MaxPr[1]
                    end
                end
            end
        end

        result_mats = Dict("minPr" => minPr_mat, "maxPr" => maxPr_mat)
    end
    total_runtime += bound_time
    timing_info["transition_bound_time_s"] = bound_time 
    @info "Bound generation time: " bound_time
    save_transition_matrices(params, result_mats)

    verification_result_mat = nothing
    exp_dir = params.experiment_directory
    if single_mode_verification
        @info "Performing safety verification on single mode for 1 step..."
        verification_time = @elapsed begin 
            imdp = create_simple_imdp(result_mats["minPr"], result_mats["maxPr"])
            horizon=1
            verification_result_mat = Globally(imdp, "safe", horizon, "$exp_dir/imdp.txt")
            save_legacy_mats(verification_result_mat, exp_dir, horizon)
        end
        timing_info["single_verification_time_s"] = verification_time  
        @info "Verification time: " verification_time
    end

    @info "Total runtime: " total_runtime
    timing_info["total_runtime_s"] = total_runtime
    save_time_info(params, timing_info)
    
    flush(logfile)
    close(logfile)

    return timing_info, verification_result_mat, result_mats
end