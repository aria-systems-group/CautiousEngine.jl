"""
Run the PIMDP synthesis
"""
function run_pimdp_synthesis(pimdp::PIMDP, pimdp_filename::String; add_opt=false)
    write_pimdp_to_file(pimdp, pimdp_filename)

    @info "Synthesizing a controller... (maximize pessimistic)"
    res_mat = run_unbounded_imdp_synthesis(pimdp_filename)

    if add_opt
        @info "Synthesizing a controller... (maximize optimistic)"
        basename = pimdp_filename[1:end-3]
        opt_pimdp_filename = "$basename-optimistic.txt"
        res_mat_opt = run_imdp_synthesis(pimdp_filename, -1, mode2="optimistic", save_mats=false)
        for j in 1:length(res_mat[:,1])
            if res_mat[j, 3] == res_mat_opt[j,3] && res_mat_opt[j,4] > res_mat[j,4] 
                res_mat[j, 2] = res_mat_opt[j,2]
                res_mat[j, 4] = res_mat_opt[j,4] 
            end
        end
    end

    return res_mat
end
