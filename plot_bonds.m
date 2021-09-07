clear all 
close all

fpath = '/Users/kyleacheson/MATLAB/SCATTERING/ROT_AVG/MeV_UED/Theory/OUTPUTS/Optimised_OUTPUT.mat';
FLAGweight = 1; % 0 = assume equal weighting, 1 = multiply by optimised weights
dt = 0.5; % time step in theory
diss_cutoff = 3.4;
diss_label = 2; % atom label in matrix of bonds
bound_label = 3;

[Q, weights, important_trajs, traj_spin] = load_fields(fpath);

[Natom, ~, Ntraj, Nts] = size(Q);
time = 0:dt:(Nts-1)*dt;

if FLAGweight == 1
    count = 0;
    for tr=1:Ntraj
        if ismember(tr, important_trajs)
            count = count + 1;
            Qi(:, :, count, :) = Q(:, :, tr, :);
        end
    end
    
    Q = Qi; clear Qi;
    Ntraj = length(important_trajs);
    weights = weights(important_trajs);
    
end

% distance matrix
for tr=1:Ntraj
    for ts=1:Nts
        for a=1:Natom
            for b=a+1:Natom
                dd = norm(Q(a,1:3,tr,ts) - Q(b,1:3,tr,ts));
                D(a,b,tr,ts) = dd;
                D(b,a,tr,ts) = dd;
            end
        end
    end
end

bond1 = zeros(Ntraj,Nts);
bond2 = zeros(Ntraj,Nts);
bond3 = zeros(Ntraj,Nts);

for tr=1:Ntraj % sort lengths into dissociative, bound and s-s
    for ts=1:Nts
      if D(1,diss_label,tr,end) > D(1,bound_label,tr,end)  
          bond1(tr,ts) = D(1,diss_label,tr,ts);   % bond 1 = dissociative bond
          bond2(tr,ts) = D(1,bound_label,tr,ts);  % bond 2 = remaining bound bond
      elseif D(1,diss_label,tr,end) < D(1,bound_label,tr,end)
          bond1(tr,ts) = D(1,bound_label,tr,ts);
          bond2(tr,ts) =  D(1,diss_label,tr,ts);
      end
      bond3(tr,ts) = D(diss_label,bound_label,tr,ts);
      costheta = (bond1(tr,ts)^2 + bond2(tr,ts)^2-bond3(tr,ts)^2)/(2*bond1(tr,ts)*bond2(tr,ts));
      angle(tr,ts) = acos(costheta)*(180/pi);
      
    end
end

for tr=1:Ntraj
    if traj_spin(tr) ~= 0
        tind = min(find(bond1(tr,:) > diss_cutoff));
        diss_times(tr) = time(tind); % time at which dissociation occurs
    else
        diss_times(tr) = time(end); % bound trajs
    end
       
end


singlet_times = diss_times(traj_spin == 1);
triplet_times = diss_times(traj_spin == 3);
singlet_weights = weights(traj_spin == 1);
triplet_weights = weights(traj_spin == 3);

norm_singlet_weights = singlet_weights./sum(singlet_weights);
norm_triplet_weights = triplet_weights./sum(triplet_weights);

singlet_time = 0;
for tr=1:length(singlet_times)
    singlet_time = singlet_time + (singlet_times(tr) * norm_singlet_weights(tr));
end
triplet_time = 0;
for tr=1:length(triplet_times)
    triplet_time = triplet_time + (triplet_times(tr) * norm_triplet_weights(tr));
end

avg_singlet_time = sum(singlet_times)./length(singlet_times);
avg_triplet_time = sum(triplet_times)./length(triplet_times);

disp(['Number Bound Trajs: ', num2str(sum(traj_spin == 0))])
disp(['Number Singlet Trajs: ', num2str(sum(traj_spin == 1))])
disp(['Number Triplet Trajs: ', num2str(sum(traj_spin == 3))])

disp(['Average Singlet Dissociation Time = ', num2str(avg_singlet_time), ' fs'])
disp(['Average Triplet Dissociation Time = ', num2str(avg_triplet_time), ' fs'])
disp(['Normalised Weighted Singlet Dissociation Time = ', num2str(singlet_time), ' fs'])
disp(['Normalised Weighted Dissociation Time = ', num2str(triplet_time), ' fs'])

T = zeros(Natom+1, Ntraj);

for tr=1:Ntraj
   
    tind = find(time == diss_times(tr));

    [peaks, locs1, width, prom] = findpeaks(bond1(tr,1:tind)); % find peaks of bound frag vibration
    T(1, tr) = sum(diff(locs1).*dt)./length(locs1); % period of vibration

    [peaks, locs2, width, prom] = findpeaks(bond2(tr,1:tind)); % find peaks of bound frag vibration
    T(2, tr) = sum(diff(locs2).*dt)./length(locs2); % period of vibration

    [peaks, locs3, width, prom] = findpeaks(bond3(tr,1:tind)); % find peaks of bound frag vibration
    T(3, tr) = sum(diff(locs3).*dt)./length(locs3); % period of vibration
    
    [peaks, locs4, width, prom] = findpeaks(angle(tr,1:tind)); % find peaks of bound frag vibration
    T(4, tr) = sum(diff(locs4).*dt)./length(locs4); % period of vibration
    
    
    figure
    subplot(4,1,1)
    plot(time(1:tind), bond1(tr, 1:tind), 'b');
    hold on
    for i=1:length(locs1)
        xline(locs1(i)*dt)
    end
    subplot(4,1,2)
    plot(time(1:tind), bond2(tr, 1:tind), 'r');
    hold on
    for i=1:length(locs2)
        xline(locs2(i)*dt)
    end
    subplot(4,1,3)
    plot(time(1:tind), bond3(tr, 1:tind), 'g');
    hold on
    for i=1:length(locs3)
        xline(locs3(i)*dt)
    end
    subplot(4,1,4)
    plot(time(1:tind), angle(tr, 1:tind), 'c');
    hold on
    for i=1:length(locs4)
        xline(locs4(i)*dt)
    end
    
    
    sgtitle(num2str(tr))

   
    
end



% for tr=1:Ntraj
% [peaks, locs, width, prom] = findpeaks(bond2(tr,:)); % find peaks of bound frag vibration
% P(tr) = sum(prom)./length(prom);
% T(tr) = sum(diff(locs).*dt)./length(locs); % period of vibration
% end

if FLAGweight == 1
    
    cs1 = zeros(1, Nts);
    cs2 = zeros(1, Nts);
    ss = zeros(1, Nts);
    ang = zeros(1, Nts);
    
    for tr=1:Ntraj
        cs1 = cs1 + (bond1(tr, :) * weights(tr));
        cs2 = cs2 + (bond2(tr, :) * weights(tr));
        ss = ss + (bond3(tr, :) * weights(tr));
        ang = ang + angle(tr, :) * weights(tr);
    end
    
else
    
    cs1 = sum(bond1, 1)/Ntraj; % average
    cs2 = sum(bond2, 1)/Ntraj;
    ss = sum(bond3, 1)/Ntraj;
    ang = sum(angle, 1)/Ntraj;
    
end


figure
subplot(4,1,1)
plot(time, cs1, 'r')
xlabel('Time (fs)')
ylim([0 4])
ylabel(['Distance (', char(197), ')'])
subplot(4,1,2)
plot(time, cs2, 'b')
xlabel('Time (fs)')
ylabel(['Distance (', char(197), ')'])
subplot(4,1,3)
plot(time, ss, 'g')
xlabel('Time (fs)')
ylim([0 8])
ylabel(['Distance (', char(197), ')'])
subplot(4,1,4)
plot(time, ang, 'c')
xlabel('Time (fs)')
ylim([0 8])
ylabel(['Distance (', char(197), ')'])

%legend('C-S$_a$ Distance', 'C-S$_b$ Distance', 'S-S Distance', 'interpreter','latex', 'Location', 'NorthWest')


