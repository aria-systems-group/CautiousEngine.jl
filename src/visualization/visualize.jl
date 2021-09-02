"""
Plot the results from file using the legacy format.
"""
function plot_results_from_file(results_path; num_dfa_states=1, plot_gp=false, min_threshold=0.95, trajectories=nothing, filename=nothing, plot_gamma=false, extents_dict=nothing, filepath=nothing)
    extents, post_extents = deserialize_region_data("$results_path/regions.bson")
    num_regions = extents.count - 1
    X = extents[extents.count]
    minx = minimum(X["x1"])
    maxx = maximum(X["x1"])
    miny = minimum(X["x2"])
    maxy = maximum(X["x2"])

    mat_files = filter(x -> occursin(r"globally-safe-.*.mat", x), readdir(results_path))

    if plot_gp
        files = readdir("$results_path/gps")
        gp_file = filter(x->occursin(".bin", x), files)[1]
        f = open("$results_path/gps/$gp_file")
        gp_set = deserialize(f)
        gp_x1 = gp_set["x1"]
        close(f)
    end

    for mat_file in mat_files

        res = matread("$results_path/$mat_file")
        k = occursin("--1", mat_file) ? -1 : parse(Int, split(mat_file, "-")[3])
        base_str = mat_file[1:end-4]
        # Plot the maximum results
        # Probabilities for each cell from 0 to 1 corresponding to regions 1 to N
        # TODO: Need to somehow generalize this skipping step
        maxPrs = res["indVmax"][1:num_dfa_states:end]
        minPrs = res["indVmin"][1:num_dfa_states:end]
        min_minPr = minimum(minPrs)
        min_maxPr = minimum(maxPrs)

        function plot_cell(extent, prob_value, norm_val)
            x = [extent["x1"][1], extent["x1"][1], extent["x1"][2], extent["x1"][2]]
            y = [extent["x2"][1], extent["x2"][2], extent["x2"][2], extent["x2"][1]]
            shape = Plots.Shape(x, y)
            plot!(shape, color=:black, fillalpha=1-prob_value, linealpha=0.0, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
        end

        plt_max = plot(#title="$k Step Maximum Safety",
                       aspect_ratio=1,
                       size=(300,300), dpi=300,
                       xlims=[minx, maxx], ylims=[miny, maxy],
                       #xlabel="X Units", ylabel="Y Units",
                       xtickfont=font(10),
                       ytickfont=font(10),
                       titlefont=font(10),
                       grid=false)
        plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
        [plot_cell(extents[i], maxPrs[i], min_maxPr) for i in 1:num_regions]
        savefig(plt_max, "$results_path/$base_str-max-heatmap.png")
        # Plot the minimum results
        plt_min = plot(#title=L"$k Step Minimum Safety",
                       aspect_ratio=1,
                       size=(300,300), dpi=300,
                       xlims=[minx, maxx], ylims=[miny, maxy],
                       #xlabel="X Units", ylabel="Y Units",
                       xtickfont=font(10),
                       ytickfont=font(10),
                       titlefont=font(10),
                       grid=false,
                       backgroundcolor=128)
        plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
        [plot_cell(extents[i], minPrs[i], min_minPr) for i in 1:num_regions]
        savefig(plt_min, "$results_path/$base_str-min-heatmap.png")

        plt_min = plot(#title=L"$k Step Minimum Safety",
                       aspect_ratio=1,
                       size=(300,300), dpi=300,
                       xlims=[minx, maxx], ylims=[miny, maxy],
                       #xlabel="X Units", ylabel="Y Units",
                       xtickfont=font(10),
                       ytickfont=font(10),
                       titlefont=font(10),
                       grid=false,
                       backgroundcolor=128)

        function plot_cell_verify(extent, min_prob_value, max_prob_value, threshold)
            x = [extent["x1"][1], extent["x1"][1], extent["x1"][2], extent["x1"][2]]
            y = [extent["x2"][1], extent["x2"][2], extent["x2"][2], extent["x2"][1]]
            shape = Plots.Shape(x, y)

            if min_prob_value >= threshold
                plot!(shape, color=:purple, fillalpha=0.10, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
            elseif max_prob_value < threshold
                plot!(shape, color=:purple, fillalpha=0.99, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
            else
                plot!(shape, color=:purple, fillalpha=0.50, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
            end
        end
        # Plot the cells
        plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
        [plot_cell_verify(extents[i], minPrs[i], maxPrs[i], min_threshold) for i in 1:num_regions]
        savefig(plt_min, "$results_path/$base_str-verification.png")
        if !isnothing(filepath)
            savefig(plt_min, "$filepath-verification.png")
        end

        if !isnothing(trajectories)
            for trajectory in trajectories
                if length(trajectory) > 1
                    [plot!([trajectory[i][1], trajectory[i+1][1]], [trajectory[i][2], trajectory[i+1][2]], color=:black, label="", linewith=2) for i=1:length(trajectory)-1]
                    scatter!([trajectory[end][1]], [trajectory[end][2]], color=:black, markershape=:star5, label="")
                end
            end

            if !isnothing(extents_dict)
                for key in keys(extents_dict)
                    for extent in extents_dict[key]
                        shape = extent_dict_to_shape(extent)
                        @info shape
                        plot!(shape, label="", color=:black, fillalpha=0.0)
                    end
                end
            end
            savefig(plt_min, "$results_path/$filename")
            if !isnothing(filepath)
                savefig(plt_min, "$filepath-sim.png")
            end
        end

        if plot_gamma
            plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
            [plot_cell(extents[i], abs(minPrs[i]-maxPrs[i]), nothing) for i in 1:num_regions]
            savefig(plt_min, "$results_path/$base_str-gamma-values.png")
        end
        
        if plot_gp
            # Plot where the data was collected and a convex hull
            data_shape = VPolygon(convex_hull([gp_x1.x[:, i] for i=1:gp_x1.nobs]))
            plot!(data_shape, fillalpha=0.25)
            [scatter!([gp_x1.x[1, i]], [gp_x1.x[2,i]], label="", color=:black) for i=1:gp_x1.nobs]
            savefig(plt_min, "$results_path/$base_str-data-hull.pdf")
        end
    end
end

function plot_2d_verification_results(results_path;num_dfa_states=1, min_threshold=0.95, extents_dict=Dict(), extents_labels_colors=Dict(), prob_plots=false, lazy_flag=false, plot_tag="", state_outlines_flag=false)
    extents, _ = deserialize_region_data("$results_path/regions.bson", lazy_flag=lazy_flag)
    num_regions = extents.count - 1 
    X = extents[extents.count]
    minx = minimum(X["x1"])
    maxx = maximum(X["x1"])
    miny = minimum(X["x2"])
    maxy = maximum(X["x2"])

    mat_files = filter(x -> occursin(r"globally-safe-.*.mat", x), readdir(results_path))

    for mat_file in mat_files
        res = matread("$results_path/$mat_file")
        base_str = mat_file[1:end-4]

        # TODO: Need to somehow generalize this skipping step
        maxPrs = res["indVmax"][1:num_dfa_states:end]
        minPrs = res["indVmin"][1:num_dfa_states:end]

        if prob_plots
            # Plot the maximumprobabilities 
            plt_max = plot(aspect_ratio=1,
                        size=(300,300), dpi=300,
                        xlims=[minx, maxx], ylims=[miny, maxy],
                        xtickfont=font(10),
                        ytickfont=font(10),
                        titlefont=font(10),
                        xticks = [minx, 0, maxx],
                        yticks = [miny, 0, maxy],
                        grid=false)
            [plot_cell(extents[i], maxPrs[i], state_outlines_flag=state_outlines_flag) for i in 1:num_regions]
            plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
            savefig(plt_max, "$results_path/$base_str$plot_tag-max-heatmap.png")

            # Plot the minimum probabilities 
            plt_min = plot(aspect_ratio=1,
                        size=(300,300), dpi=300,
                        xlims=[minx, maxx], ylims=[miny, maxy],
                        xtickfont=font(10),
                        ytickfont=font(10),
                        titlefont=font(10),
                        xticks = [minx, 0, maxx],
                        yticks = [miny, 0, maxy],
                        grid=false,
                        backgroundcolor=128)
            [plot_cell(extents[i], minPrs[i], state_outlines_flag=state_outlines_flag) for i in 1:num_regions]
            plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
            savefig(plt_min, "$results_path/$base_str$plot_tag-min-heatmap.png")
        end

        plt_verification = plot(aspect_ratio=1,
                                size=(300,300), dpi=300,
                                xlims=[minx, maxx], ylims=[miny, maxy],
                                xtickfont=font(10),
                                ytickfont=font(10),
                                titlefont=font(10),
                                xticks = [minx, 0, maxx],
                                yticks = [miny, 0, maxy],
                                grid=false,
                                backgroundcolor=128)

        # Plot the cells
        [plot_cell_verify(extents[i], minPrs[i], maxPrs[i], min_threshold, state_outlines_flag=state_outlines_flag) for i in 1:num_regions]
        plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")

        # Plot extents with labels
        for key in keys(extents_dict)
            labelled_extents = extents_dict[key]
            for extent in labelled_extents
                plot_labelled_extent_outline(extent, key, color=extents_labels_colors[key])
            end
        end

        savefig(plt_verification, "$results_path/$base_str$plot_tag-verification.png")
    end
end

"""
Plot a labelled extent
"""
function plot_labelled_extent_outline(extent, label; color=:white)
    xvals = extent["x1"]
    yvals = extent["x2"]
    plot!(extent_dict_to_shape(extent), fillalpha=0, linecolor=:black, linewidth=1, label="")
    annotate!((minimum(xvals)+maximum(xvals))/2, (minimum(yvals)+maximum(yvals))/2, text(label, color, :center, 10,))
end

"""
Plot a cell with alpha associated with probability value
"""
function plot_cell(extent, prob_value; state_outlines_flag=false)
    x = [extent["x1"][1], extent["x1"][1], extent["x1"][2], extent["x1"][2]]
    y = [extent["x2"][1], extent["x2"][2], extent["x2"][2], extent["x2"][1]]
    shape = Plots.Shape(x, y)
    linealpha = state_outlines_flag ? 1.0 : 0.0
    plot!(shape, color=:black, fillalpha=1-prob_value, linewidth=1, linecolor=:black, linealpha=linealpha, foreground_color_border=:white, foreground_color_axis=:white, label="")
end

"""
Determine the cell color based on the verification probability threshold.
"""
function plot_cell_verify(extent, min_prob_value, max_prob_value, threshold; state_outlines_flag=false)
    x = [extent["x1"][1], extent["x1"][1], extent["x1"][2], extent["x1"][2]]
    y = [extent["x2"][1], extent["x2"][2], extent["x2"][2], extent["x2"][1]]
    shape = Plots.Shape(x, y)

    linealpha = state_outlines_flag ? 1.0 : 0.0

    if min_prob_value >= threshold
        plot!(shape, color=:purple, fillalpha=0.10, linewidth=1, linecolor=:black, linealpha=linealpha, foreground_color_border=:white, foreground_color_axis=:white, label="")
    elseif max_prob_value < threshold
        plot!(shape, color=:purple, fillalpha=0.99, linewidth=1, linecolor=:black, linealpha=linealpha, foreground_color_border=:white, foreground_color_axis=:white, label="")
    else
        plot!(shape, color=:purple, fillalpha=0.50, linewidth=1, linecolor=:black, linealpha=linealpha, foreground_color_border=:white, foreground_color_axis=:white, label="")
    end
end

function plot_gamma_value(results_path, min_ResMat, max_ResMat; num_dfa_states=1)
    function plot_cell(extent, prob_value)
        x = [extent["x1"][1], extent["x1"][1], extent["x1"][2], extent["x1"][2]]
        y = [extent["x2"][1], extent["x2"][2], extent["x2"][2], extent["x2"][1]]
        shape = Plots.Shape(x, y)
        plot!(shape, color=:black, fillalpha=1-prob_value, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
    end

    extents, post_extents = deserialize_region_data("$results_path/regions.bson") 
    num_regions = extents.count - 1
    X = extents[extents.count]
    minx = minimum(X["x1"])
    maxx = maximum(X["x1"])
    miny = minimum(X["x2"])
    maxy = maximum(X["x2"])

    minPrs = min_ResMat[1:num_dfa_states:end, 3]
    maxPrs = max_ResMat[1:num_dfa_states:end, 4]

    base_str = "gamma-values" 
    plt = plot(#title=L"$k Step Minimum Safety",
                       aspect_ratio=1,
                       size=(300,300), dpi=300,
                       xlims=[minx, maxx], ylims=[miny, maxy],
                       #xlabel="X Units", ylabel="Y Units",
                       xtickfont=font(10),
                       ytickfont=font(10),
                       titlefont=font(10),
                       grid=false,
                       backgroundcolor=128)
    plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
    [plot_cell(extents[i], abs(minPrs[i]-maxPrs[i])) for i in 1:num_regions]
    savefig(plt, "$results_path/$base_str.png")
end

function plot_gp_fields(results_path, dyn_fn)
    files = readdir("$results_path/gps")
    gp_file = filter(x->occursin(".bin", x), files)[1]
    f = open("$results_path/gps/$gp_file")
    gp_save_info = deserialize(f)
    close(f)
    extents, post_extents = deserialize_region_data("$results_path/regions.bson")

    if length(gp_save_info) == 1
        plot_gp_1d(gp_set, region_data, results_path, dyn_fn)
        return
    elseif length(gp_save_info) > 2
        plot_gp_field_slice(results_path, dyn_fn, 0.0; nl_flag=false) 
        return
    end

    X = extents[extents.count]
    minx = minimum(X["x1"])
    maxx = maximum(X["x1"])
    miny = minimum(X["x2"])
    maxy = maximum(X["x2"])

    x_range = minx:0.25:maxx
    y_range = miny:0.25:maxy
    meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
    x, y = meshgrid(x_range, y_range)
    gp_x1 = gp_save_info[:gp_set]["x1"]
    gp_x2 = gp_save_info[:gp_set]["x2"]
    delta_x_gp = [x[i] - predict_y(gp_x1, hcat([x[i], y[i]]))[1][1] for i in 1:length(x)]
    delta_y_gp = [y[i] - predict_y(gp_x2, hcat([x[i], y[i]]))[1][1] for i in 1:length(x)]
    u_gp = delta_x_gp*0.5
    v_gp = delta_y_gp*0.5
    plt1 = plot(size=(300,300), dpi=600, aspect_ratio=1, grid=false, foreground_color_border=:white, foreground_color_axis=:white, 
                xlims = [minx, maxx], ylims = [miny, maxy],
                xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], 
                xtickfont=font(10), ytickfont=font(10), titlefont=font(10))
    plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="", linealpha=0.1)
    quiver!(x, y, quiver=(-u_gp, -v_gp), color=:black, linewidth=1.5)
    plt2 = plot(size=(300,300), dpi=600, aspect_ratio=1, grid=false, foreground_color_border=:white, foreground_color_axis=:white,
                xlims = [minx, maxx], ylims = [miny, maxy],
                xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], 
                xtickfont=font(10), ytickfont=font(10), titlefont=font(10))
    plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="", linealpha=0.1)

    u_true = 0.5*[x[i] - dyn_fn([x[i], y[i]])[1] for i in 1:length(x)] 
    v_true = 0.5*[y[i] - dyn_fn([x[i], y[i]])[2] for i in 1:length(x)] 
    quiver!(x, y, quiver=(-u_true, -v_true), color=:black, linewidth=1.5)
    plt = plot(plt1, plt2, layout=(1,2))
    savefig(plt1, "$results_path/gp-vector-fields.png")
    savefig(plt2, "$results_path/true-vector-fields.png")

    plt = plot(gp_x1, var=true, title="Variance of GP1", fill=true)
    savefig(plt, "$results_path/gp1-var.pdf")
    plt = plot(gp_x2, var=true, title="Variance of GP2", fill=true)
    savefig(plt, "$results_path/gp2-var.pdf")

    # Plot where the data was collected and a convex hull
    data_shape = VPolygon(convex_hull([gp_x1.x[:, i] for i=1:gp_x1.nobs]))
    plt = plot(data_shape, fillalpha=0.5)
    [scatter!([gp_x1.x[1, i]], [gp_x1.x[2,i]], label="", color=:black) for i=1:gp_x1.nobs]
    savefig(plt, "$results_path/data-hull.pdf")
end

function plot_gp_1d(gp_set, region_data, results_path, dyn_fn)
    extents = region_data[:extents]
    X = extents[extents.count]
    minx = minimum(X["x1"])
    maxx = maximum(X["x1"])
    x = collect(minx:0.25:maxx)
    gp_x1 = gp_set["x1"]
    delta_x_gp = [x[i] - predict_y(gp_x1,[x[i]])[1][1] for i in 1:length(x)]

    plt1 = plot(gp_x1)
    savefig(plt1, "$results_path/gp-vector-fields.png")
    plt = plot(gp_x1, var=true, title="Variance of GP1")
    savefig(plt, "$results_path/gp1-var.pdf")
end

function plot_gp_field_slice(results_path, dyn_fn, x3; nl_flag=false)
    f = open("$results_path/gaussian_processes_set.bin")
    gp_set = deserialize(f)
    close(f)
    #A_mats = collect(keys(gp_set))
    #A = A_mats[1]
    #A = [0.9 0.4; -0.4 0.5] # Rotation
    #gp_dict = gp_set[A_mats[1]]

    extents, post_extents = deserialize_region_data("$results_path/regions.bson")
    X = extents[extents.count]
    minx = minimum(X["x1"])
    maxx = maximum(X["x1"])
    miny = minimum(X["x2"])
    maxy = maximum(X["x2"])

    x_range = minx:0.25:maxx
    y_range = miny:0.25:maxy
    meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
    y, x = meshgrid(x_range, y_range)
    gp_x1 = gp_set["x1"]
    gp_x2 = gp_set["x2"]
    delta_x_gp = [x[i] - predict_y(gp_x1, hcat([x[i], y[i], x3]))[1][1] for i in 1:length(x)]
    delta_y_gp = [y[i] - predict_y(gp_x2, hcat([x[i], y[i], x3]))[1][1] for i in 1:length(x)]
    u_gp = delta_x_gp*0.5
    v_gp = delta_y_gp*0.5
    plt1 = plot(size=(300,300), dpi=600, aspect_ratio=1, grid=false, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], xtickfont=font(10), ytickfont=font(10), titlefont=font(10))
    plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
    quiver!(x, y, quiver=(-u_gp, -v_gp), color=:black, linewidth=1.5)
    plt2 = plot(size=(300,300), dpi=600, aspect_ratio=1, grid=false, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], xtickfont=font(10), ytickfont=font(10), titlefont=font(10))
    plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="", linealpha=0.1)

    u_true = 0.5*[x[i] - dyn_fn([x[i], y[i], x3])[1] for i in 1:length(x)] 
    v_true = 0.5*[y[i] - dyn_fn([x[i], y[i], x3])[2] for i in 1:length(x)] 
    quiver!(x, y, quiver=(-u_true, -v_true), color=:black, linewidth=1.5)
    plt = plot(plt1, plt2, layout=(1,2))
    savefig(plt1, "$results_path/gp-vector-fields-$x3.png")
    savefig(plt2, "$results_path/true-vector-fields-$x3.png")

    # plt = plot(gp_x1, var=true, title="Variance of GP1", fill=true)
    # savefig(plt, "$results_path/gp1-var.pdf")
    # plt = plot(gp_x2, var=true, title="Variance of GP2", fill=true)
    # savefig(plt, "$results_path/gp2-var.pdf")

    # Plot where the data was collected and a convex hull
    # data_shape = VPolygon(convex_hull([gp_x1.x[:, i] for i=1:gp_x1.nobs]))
    # plt = plot(data_shape, fillalpha=0.5)
    # [scatter!([gp_x1.x[1, i]], [gp_x1.x[2,i]], label="", color=:black) for i=1:gp_x1.nobs]
    # savefig(plt, "$results_path/data-hull.pdf")
end

function plot_synthesis_results(results_path, res_mat, imdp, dfa, pimdp; plot_gp=false, min_threshold=0.95, trajectories=nothing, filename=nothing, old_mins=nothing)
    extents, post_extents = deserialize_region_data("$results_path/regions.bson")
    num_regions = extents.count - 1
    X = extents[extents.count]

    if X.count == 2
        plot_results_from_file(results_path, num_dfa_states=length(dfa.states), trajectories=trajectories, filename=filename)
        return
    end

    minx = minimum(X["x1"])
    maxx = maximum(X["x1"])
    miny = minimum(X["x2"])
    maxy = maximum(X["x2"])

    # mat_files = filter(x -> occursin(r"globally-safe-.*.mat", x), readdir(results_path))

    # if plot_gp
    #     f = open("$results_path/gaussian_processes_set.bin")
    #     gp_set = deserialize(f)
    #     gp_x1 = gp_set["x1"]
    #     close(f)
    # end

    # for mat_file in mat_files

    # res = matread("$results_path/$mat_file")
    # k = occursin("--1", mat_file) ? -1 : parse(Int, split(mat_file, "-")[3])
    # base_str = mat_file[1:end-4]
    # Plot the maximum results
    # Probabilities for each cell from 0 to 1 corresponding to regions 1 to N
    # TODO: Need to somehow generalize this skipping step
    Qsize = length(dfa.states)
    maxPrs = res_mat[1:Qsize:end, 4]
    minPrs = res_mat[1:Qsize:end, 3]
    min_minPr = minimum(minPrs)
    min_maxPr = minimum(maxPrs)

    # Get the number of x3 extents
    x3_extents = []
    
    for region_id in keys(extents)
        if !in(extents[region_id]["x3"], x3_extents)
            push!(x3_extents, extents[region_id]["x3"])
        end
    end

    function plot_cell(extent, prob_value, norm_val)
        x = [extent["x1"][1], extent["x1"][1], extent["x1"][2], extent["x1"][2]]
        y = [extent["x2"][1], extent["x2"][2], extent["x2"][2], extent["x2"][1]]
        shape = Plots.Shape(x, y)
        plot!(shape, color=:black, fillalpha=1-prob_value, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
    end

    for (j, x3_extent) in enumerate(x3_extents)
        indeces = []
        for i=1:num_regions
            if extents[i]["x3"] == x3_extent
                push!(indeces, i)
            end
        end
        base_str = "synthesis-res-x3$j"

        plt_max = plot(#title="$k Step Maximum Safety",
                        aspect_ratio=1,
                        size=(300,300), dpi=300,
                        xlims=[minx, maxx], ylims=[miny, maxy],
                        #xlabel="X Units", ylabel="Y Units",
                        xtickfont=font(10),
                        ytickfont=font(10),
                        titlefont=font(10),
                        grid=false)
        plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
        [plot_cell(extents[i], maxPrs[i], min_maxPr) for i in indeces]
        savefig(plt_max, "$results_path/$base_str-max-heatmap.png")
        # Plot the minimum results
        plt_min = plot(#title=L"$k Step Minimum Safety",
                        aspect_ratio=1,
                        size=(300,300), dpi=300,
                        xlims=[minx, maxx], ylims=[miny, maxy],
                        #xlabel="X Units", ylabel="Y Units",
                        xtickfont=font(10),
                        ytickfont=font(10),
                        titlefont=font(10),
                        grid=false,
                        backgroundcolor=128)
        plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
        [plot_cell(extents[i], minPrs[i], min_minPr) for i in indeces]
        savefig(plt_min, "$results_path/$base_str-min-heatmap.png")

        plt_min = plot(#title=L"$k Step Minimum Safety",
                        aspect_ratio=1,
                        size=(300,300), dpi=300,
                        xlims=[minx, maxx], ylims=[miny, maxy],
                        #xlabel="X Units", ylabel="Y Units",
                        xtickfont=font(10),
                        ytickfont=font(10),
                        titlefont=font(10),
                        grid=false,
                        backgroundcolor=128)

        function plot_cell_verify(extent, min_prob_value, max_prob_value, threshold)
            x = [extent["x1"][1], extent["x1"][1], extent["x1"][2], extent["x1"][2]]
            y = [extent["x2"][1], extent["x2"][2], extent["x2"][2], extent["x2"][1]]
            shape = Plots.Shape(x, y)

            if min_prob_value >= threshold
                plot!(shape, color=:purple, fillalpha=0.10, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
            elseif max_prob_value < threshold
                plot!(shape, color=:purple, fillalpha=0.99, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
            else
                plot!(shape, color=:purple, fillalpha=0.50, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
            end
        end
        # Plot the cells
        plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
        [plot_cell_verify(extents[i], minPrs[i], maxPrs[i], min_threshold) for i in indeces]
        savefig(plt_min, "$results_path/$base_str-verification.png")

        if !isnothing(trajectories)
            for trajectory in trajectories
                [plot!([trajectory[i][1], trajectory[i+1][1]], [trajectory[i][2], trajectory[i+1][2]], color=:black, label="", linewith=2) for i=1:length(trajectory)-1]
                scatter!([trajectory[end][1]], [trajectory[end][2]], color=:black, markershape=:star5)
            end
            savefig(plt_min, "$results_path/$filename")
        end
    end
    # end
end

function extent_dict_to_shape(extent)
    xmin = minimum(extent["x1"])
    xmax = maximum(extent["x1"])
    ymin = minimum(extent["x2"])
    ymax= maximum(extent["x2"])

    x = [xmin, xmin, xmax, xmax]
    y = [ymin, ymax, ymax, ymin]

    return Plots.Shape(x, y)
end

"""
Plot the results from file using the legacy format.
"""
function plot_online_results(extents, res_mat, plots_basename, num_dfa_states; min_threshold=0.95, trajectories=nothing, trajectory_dict=nothing, extents_dict=nothing, online_data_dict=nothing)
    num_regions = extents.count - 1
    X = extents[extents.count]
    minx = minimum(X["x1"])
    maxx = maximum(X["x1"])
    miny = minimum(X["x2"])
    maxy = maximum(X["x2"])

    # if plot_gp
    #     files = readdir("$results_path/gps")
    #     gp_file = filter(x->occursin(".bin", x), files)[1]
    #     f = open("$results_path/gps/$gp_file")
    #     gp_set = deserialize(f)
    #     gp_x1 = gp_set["x1"]
    #     close(f)
    # end

    # Plot the maximum results
    maxPrs = res_mat[1:num_dfa_states:end,4]
    minPrs = res_mat[1:num_dfa_states:end,3] 
    min_minPr = minimum(minPrs)
    min_maxPr = minimum(maxPrs)

    function plot_cell(extent, prob_value, norm_val)
        x = [extent["x1"][1], extent["x1"][1], extent["x1"][2], extent["x1"][2]]
        y = [extent["x2"][1], extent["x2"][2], extent["x2"][2], extent["x2"][1]]
        shape = Plots.Shape(x, y)
        plot!(shape, color=:black, fillalpha=1-prob_value, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
    end

    plt_max = plot(#title="$k Step Maximum Safety",
                    aspect_ratio=1,
                    size=(300,300), dpi=300,
                    xlims=[minx, maxx], ylims=[miny, maxy],
                    #xlabel="X Units", ylabel="Y Units",
                    xtickfont=font(10),
                    ytickfont=font(10),
                    titlefont=font(10),
                    grid=false)
    plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
    [plot_cell(extents[i], maxPrs[i], min_maxPr) for i in 1:num_regions]
    savefig(plt_max, "$plots_basename-max-heatmap.png")
    # Plot the minimum results
    plt_min = plot(#title=L"$k Step Minimum Safety",
                    aspect_ratio=1,
                    size=(300,300), dpi=300,
                    xlims=[minx, maxx], ylims=[miny, maxy],
                    #xlabel="X Units", ylabel="Y Units",
                    xtickfont=font(10),
                    ytickfont=font(10),
                    titlefont=font(10),
                    grid=false,
                    backgroundcolor=128)
    plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
    [plot_cell(extents[i], minPrs[i], min_minPr) for i in 1:num_regions]

    if !isnothing(extents_dict)
        for key in keys(extents_dict)
            for extent in extents_dict[key]
                shape = extent_dict_to_shape(extent)
                if key == "!a∧!b∧c"
                    plot!(shape, color=:black, label="", fillalpha=1.0)
                else
                    plot!(shape, color=:black, label="", fillalpha=0.0) 
                end
            end
        end
    end

    savefig(plt_min, "$plots_basename-min-heatmap.png")

    plt_min = plot(#title=L"$k Step Minimum Safety",
                    aspect_ratio=1,
                    size=(300,300), dpi=300,
                    xlims=[minx, maxx], ylims=[miny, maxy],
                    #xlabel="X Units", ylabel="Y Units",
                    xtickfont=font(10),
                    ytickfont=font(10),
                    titlefont=font(10),
                    grid=false,
                    backgroundcolor=128)

    function plot_cell_verify(extent, min_prob_value, max_prob_value, threshold)
        x = [extent["x1"][1], extent["x1"][1], extent["x1"][2], extent["x1"][2]]
        y = [extent["x2"][1], extent["x2"][2], extent["x2"][2], extent["x2"][1]]
        shape = Plots.Shape(x, y)

        # if min_prob_value >= threshold
        plot!(shape, color=:black, fillalpha=min_prob_value/2, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
        # elseif max_prob_value < threshold
            # plot!(shape, color=:purple, fillalpha=0.99, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
        # else
            # plot!(shape, color=:purple, fillalpha=0.50, linealpha=0.05, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], label="")
        # end
    end
    # Plot the cells
    plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, linewidth=2, label="")
    [plot_cell_verify(extents[i], minPrs[i], maxPrs[i], min_threshold) for i in 1:num_regions]

    # savefig(plt_min, "$plots_basename-verification.png")

    if !isnothing(trajectories)
        for trajectory in trajectories
            if length(trajectory) > 1
                [plot!([trajectory[i][1], trajectory[i+1][1]], [trajectory[i][2], trajectory[i+1][2]], color=:black, label="", linewith=2) for i=1:length(trajectory)-1]
                scatter!([trajectory[end][1]], [trajectory[end][2]], color=:black, markershape=:star5, label="")
            end
        end

        savefig(plt_min, "$plots_basename-online-sim.png")
    end

    if !isnothing(extents_dict)
        for key in keys(extents_dict)
            for extent in extents_dict[key]
                shape = extent_dict_to_shape(extent)
                if key == "!a∧!b∧c"
                    plot!(shape, color=:black, label="", fillalpha=1.0)
                else
                    plot!(shape, color=:black, label="", fillalpha=0.0) 
                end
            end
        end
    end

    if !isnothing(trajectory_dict)
        for (jj, key) in enumerate(keys(trajectory_dict))
            trajectories = trajectory_dict[key]
            for trajectory in trajectories
                if length(trajectory) > 1
                    [plot!([trajectory[i][1], trajectory[i+1][1]], [trajectory[i][2], trajectory[i+1][2]],  color=jj, label="", linewith=2) for i=1:length(trajectory)-1]
                    scatter!([trajectory[end][1]], [trajectory[end][2]], color=jj, markershape=:star5, label="")
                end 
            end
            scatter!([trajectories[1][1][1]], [trajectories[1][1][2]], markershape=:circle, color=jj, label=key)
        end
        savefig(plt_min, "$plots_basename-online-sim-dict.png")
    end

    if !isnothing(online_data_dict) 
        # Plot where the data was collected and a convex hull
        x_data = []
        for i=1:online_data_dict.count
            for dp in online_data_dict[i]["x_online"]
                push!(x_data, dp)
            end
        end
        # data_shape = VPolygon(convex_hull(x_data))
        # plot!(data_shape, fillalpha=0.25)
        [scatter!([x_data[i][1]], [x_data[i][2]], label="", color=:black) for i=1:length(x_data)]
        savefig(plt_min, "$plots_basename-online-data-hull.png")
    end
end

function plot_extent_outlines(extents, extents_dict, plots_basename)
    num_regions = extents.count - 1
    X = extents[extents.count]
    minx = minimum(X["x1"])
    maxx = maximum(X["x1"])
    miny = minimum(X["x2"])
    maxy = maximum(X["x2"])

    plt_min = plot(#title=L"$k Step Minimum Safety",
                    aspect_ratio=1,
                    size=(300,300), dpi=300,
                    xlims=[minx, maxx], ylims=[miny, maxy],
                    #xlabel="X Units", ylabel="Y Units",
                    xtickfont=font(10),
                    ytickfont=font(10),
                    titlefont=font(10),
                    grid=false,
                    backgroundcolor=128)
    plot!(Plots.Shape([minx, minx, maxx, maxx], [miny, maxy, maxy, miny]), fillalpha=0, linecolor=:black, foreground_color_border=:white, foreground_color_axis=:white, xticks = [minx, 0, maxx], yticks = [miny, 0, maxy], linewidth=2, label="")

    for key in keys(extents_dict)
        for extent in extents_dict[key]
            shape = extent_dict_to_shape(extent)
            if key == "!a∧!b∧c"
                plot!(shape, color=:black, label="", fillalpha=1.0)
            else
                plot!(shape, color=:black, label="", fillalpha=0.0) 
            end
        end
    end
    savefig(plt_min, "$plots_basename-extent-dict.png")

end
