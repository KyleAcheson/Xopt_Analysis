function [Q, weights, important_trajs, traj_spin] = load_fields(fpath)

    load(fpath);
    
    Q = geometries;
    weights = weights;
    important_trajs = trajs_gt1;
    traj_spin = traj_spin_gt1;

end
