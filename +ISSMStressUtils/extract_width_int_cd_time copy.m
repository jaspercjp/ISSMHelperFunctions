function [rxx_t_arr, rxx_it_t_arr] = extract_width_int_stress_time(md, x_coords)
    num_t_outputs = length(md.results.TransientSolution);
    rxx_t_arr = nan*ones(num_t_outputs, length(x_coords));
    rxx_it_t_arr = nan*ones(num_t_outputs, length(x_coords));
    for t=1:10:length(md.results.TransientSolution)
        cd_t_arr(t,:) = ISSMStressUtils.extract_width_int_cd(md, x_coords, t);
    end
end