function set_time_limit!(model::Model, total_time_limit::Float64, starting_time::UInt)
    set_time_limit_sec(model, maximum([0.0, total_time_limit - (time_ns() - starting_time) / 1.0e9]))
end
