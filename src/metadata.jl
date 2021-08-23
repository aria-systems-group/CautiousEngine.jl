"
Creates the experiment directory form the parameter structure.
"
function create_experiment_directory(params)
    if !isnothing(params.data_params)
        data_tag = @sprintf("m%d-σ%1.3f-rs%d", params.data_params.data_num, params.data_params.noise_sigma, params.random_seed)
    else
        data_tag = "known-system"
    end
    for dim_key in keys(params.discretization_step)
        data_tag = @sprintf("%s-%0.3f", data_tag, params.discretization_step[dim_key])
    end
    exp_dir = @sprintf("%s/modes/%s/%s", params.experiment_directory, params.system_params.mode_tag, data_tag)
    !isdir(exp_dir) && mkpath(exp_dir)
    return exp_dir
end

function initialize_log(exp_dir; logging=Logging.Info)
    glogger = SimpleLogger(stdout, logging)
    logfile = open("$exp_dir/log.txt", "w+")
    text_logger = SimpleLogger(logfile, logging)
    demux_logger = TeeLogger(glogger, text_logger)
    global_logger(demux_logger)
    @info "Experiment directory: " exp_dir
    return logfile
end

function save_time_info(exp_dir, time_info)
    time_filename = "$exp_dir/time_info.bson"
    bson(time_filename, time_info)
end
