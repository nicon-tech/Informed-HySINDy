%% validation for hopping model
% run after Ex_Hopping.m

% =========================================================================
% run new data for validation step.
% =========================================================================
options_test = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_diverge);

% start with mass slightly below it's equilibrium position
% velocity downward
yinitvec_val = [0.84 -0.11];
%     0.77 -0.12
%     0.83 -0.13
%     0.79 -0.10
%     0.82 -0.11];
p.yinitvec = yinitvec_val;

[data]  = RUN_Hopper_training(p,0);

t_val_spring = data.tout_spring;
y_val_spring = data.yout_spring + rand(size(data.yout_spring))*p.eps;
a_val_spring = data.aout_spring;

t_val_flight = data.tout_flight;
y_val_flight = data.yout_flight + rand(size(data.yout_flight))*p.eps;
a_val_flight = data.aout_flight;

%% Find cluster in validation data that is nearest the centroid of the original clusters
% initialize a structure library
II_Xi_aic_val_spring =[]; ind_all_aic_val_spring = [];

Xicomb = Cluster_spring.Xistruct.Xicomb;
clear x2 y2
x2= y_val_spring(1:end-1,1); % height
y2= y_val_spring(1:end-1,2); % velocity

% % %     dimMan = 2; % dimension of the Manifold
% % %     numOfPts = length(x2); % number of points in the data
k_val_spring = 30; % number of points in validation spring cluster
idx_xy2 = knnsearch([x2 y2],...
        [min(x2) 0],'K',k_val_spring);
    idx_xy2 = idx_xy2(:,2:end);
% indstart_val_spring = 1; %test the model prediction comparing data_val with model, starting from a chosen position
% idx_xy2 = indstart_val_spring:k-1;

val_ntimes_spring = floor(0.8/dt); %length of comparision not including ICs % NOTE_NN: How do I choose the value that should be has val_ntimes? %% take out from for loop ii
data = [x2 y2]; %% take out from for loop ii
[xA, x0clust, tvectest] = buildinit_fromcluster(val_ntimes_spring, 1, idx_xy2, data, dt, k_val_spring);

val.tA =tvectest;
val.xA =xA;
val.x0 =x0clust;
val.options = options_test;

clear IC aic_c aic_cv tdivmax tdivmin nodynam_any tdivavg tdiv
clear nodynam abserror_bf_switch mean_e numterms
for kk = 1:length(Xicomb)
    disp('Cluster spring')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_spring = Xicomb{kk}
    
    clear  savetB savexB
    [savetB, savexB] = validation_sims(Xi_spring, Thetalib_spring ,val, plottag);
    
    for ll = 1:length(xA)
        clear tdiv1_val_spring nodynam1 abserror_val_spring ...
            abserror_avg1_val_spring RMSE_val_spring xAcomp xBcomp
        [tlength, nterms] = size(savexB{ll});
        [tdiv1_val_spring, ~, abserror_avg1_val_spring, RMSE_valspring] = calc_tlength(xA{ll}, savexB{ll}, val);
        tdiv_val_spring(ll,:) = tdiv1_val_spring;
        abserror_bf_switch_val_spring(ll,:) = abserror_avg1_val_spring; %abs error befor that the switching takes place
    end
    abserror_val_spring = [abserror_bf_switch_val_spring(:,2); abserror_bf_switch_val_spring(:,1)];
    IC_val_val_spring{kk} = ICcalculations(abserror_val_spring, nnz(Xi_spring),size(x0clust,2));
    aic_c_val_spring(kk,:) = IC_val_val_spring{kk}.aic_c;
    numterms_val_spring(kk) = nnz(Xi_spring);
    tdivmax_val_spring(kk,:) = max(tdiv_val_spring);
    tdivmin_val_spring(kk,:) = min(tdiv_val_spring);
    tdivavg_val_spring(kk,:) = mean(tdiv_val_spring);
    mean_e_val_spring(kk,:) = mean(abserror_bf_switch_val_spring);
    valsims_val_spring{kk}.savetB = savetB;
    valsims_val_spring{kk}.savexB = savexB;
    valsims_val_spring{kk}.val = val;
    
end

minIC_val_spring = min(aic_c_val_spring);
Xi_select_val_spring = find( abs(aic_c_val_spring-minIC_val_spring)<3)
Xicomb{Xi_select_val_spring}
% Make a library of models below AIC_c threshold
[II_Xitemp_val_spring, II_ind_val_spring] = build_Xi_Library(Xicomb(Xi_select_val_spring), II_Xi_aic_val_spring);
II_Xi_aic_val_spring = II_Xitemp_val_spring;
ind_all_aic_val_spring = [ind_all_aic_val_spring; II_ind_val_spring]

Cluster_val_spring.tdivmax = tdivmax_val_spring;
Cluster_val_spring.abserror_bf = mean_e_val_spring;
Cluster_val_spring.Xi_ind = Xi_select_val_spring;
Cluster_val_spring.IC = IC_val_val_spring;
Cluster_val_spring.aic_c =aic_c_val_spring;
Cluster_val_spring.numterms =numterms_val_spring;
Cluster_val_spring.valsims = valsims_val_spring;
% 
% dateformatout = 'mmddyyyy';
% save([datestr(now, dateformatout) 'Hopping_validation_eps1e_6.mat'])

