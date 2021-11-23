%% Hopping - supervised learning
% Run Hopping Model
% run this file first, then 
% Ex_Hopping_validation.m
% and
% Ex_Hopping_PLOTS.m

close all
clear all
clc

p.plottag = 1;
addpath('hopper'); addpath('SINDy'); addpath('validation');

%% Training Hopping
%The implementation of the training phase is shown below
%% data generation
% kappa = k*y0/m*g  % nondimensional parameter which represents the balance between the spring and gravity forces
p.kappa = 10;

thick_unddata = 0.35; % range undefined data position

p.eps = 1e-6; % noise level

rng(1)
% In some situations, setting the seed alone will not guarantee the same
% results. This is because the generator that the random number functions
% draw from might be different than you expect when your code executes.
% For long-term repeatability, specify the seed and the generator type together.
% For example, the following code sets the seed to 1 and the generator to Mersenne Twister.
% rng(1,'twister');

% variables are: time, t, vertical position, y, vertical velocity, yp

tend = 5;
dt = 0.033; %0.033
p.tspan = 0:dt:tend;

% start with mass slightly below it's equilibrium position
% velocity downward
p.yinitvec = [0.8 -0.1
    0.78 -0.1
    0.82 -0.1];
% p.yinitvec = [0.8 -0.1
%     0.78 -0.1
%     0.82 -0.1];

p.options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_hopper);
% The root finding mechanism employed by the ODE solver in conjunction with
% the event function has these limitations:
% If a terminal event occurs during the first step of the integration, then
% the solver registers the event as nonterminal and continues integrating.
% If more than one terminal event occurs during the first step, then only
% the first event registers and the solver continues integrating.
% Zeros are determined by sign crossings between steps. Therefore, zeros
% with an even number of crossings between steps can be missed.
% If the solver steps past events, try reducing RelTol and AbsTol to
% improve accuracy. Alternatively, set MaxStep to place an upper bound on
% the step size. Adjusting tspan does not change the steps taken by the solver.

[data]  = RUN_Hopper_training(p,thick_unddata);

t_out_und = data.tout_und;
y_out_und = data.yout_und + rand(size(data.yout_und))*p.eps;
a_out_und = data.aout_und;

t_out_spring = data.tout_spring;
y_out_spring = data.yout_spring + rand(size(data.yout_spring))*p.eps;
a_out_spring = data.aout_spring;

t_out_flight = data.tout_flight;
y_out_flight = data.yout_flight + rand(size(data.yout_flight))*p.eps;
a_out_flight = data.aout_flight;

%No clustered data
[data]  = RUN_Hopper(p);
y_out = data.yout;
a_out = data.aout;
t_out = data.tout;

%%
figure
plot(t_out_und,y_out_und(:,1),'r.')
hold on
plot(t_out_spring,y_out_spring(:,1),'k.')
plot(t_out_flight,y_out_flight(:,1),'b.')
% xlimits = [];
% ylimits = [];
% axis([xlimits ylimits])
l = get(gca, 'Children');
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$x$ [ ]','Interpreter','latex')
legend('und.area','springing','flying')
title('time-position data')
set(gca,'FontSize',9)
set(l, 'linewidth', 1.5)
%
figure
plot(t_out_und,y_out_und(:,2),'r.')
hold on
plot(t_out_spring,y_out_spring(:,2),'k.')
plot(t_out_flight,y_out_flight(:,2),'b.')
% xlimits = [];
% ylimits = [];
% axis([xlimits ylimits])
l = get(gca, 'Children');
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\dot{x}$ [ ]','Interpreter','latex')
legend('und.area','springing','flying')
title('time-velocity data')
set(gca,'FontSize',9)
set(l, 'linewidth', 1.5)
%
figure
plot(y_out_und(:,1),y_out_und(:,2),'r.')
hold on
plot(y_out_spring(:,1),y_out_spring(:,2),'k.')
plot(y_out_flight(:,1),y_out_flight(:,2),'b.')
plot([1 1],[min(y_out_und(:,2)) max(y_out_und(:,2))])
% xlimits = [];
% ylimits = [];
% axis([xlimits ylimits])
l = get(gca, 'Children');
xlabel('$x$ [s]','Interpreter','latex')
ylabel('$\dot{x}$ [ ]','Interpreter','latex')
legend('und.area','springing','flying')
title('phase space data')
set(gca,'FontSize',9)
set(l, 'linewidth', 1.5)



%% undefined data
tStart = tic;

plottag =0;
dimMan = 2; % dimension of the Manifold

% Run CCM on y_out
clear x y
x= y_out_und(1:end-1,1); % height
y= y_out_und(1:end-1,2); % velocity

numOfPts = length(x); % number of points in the data
%k = dimMan+2;
k = 6; % number of points in a cluster

idx_xy = knnsearch([x/max(x) y/max(y)],[x/max(x) y/max(y)],'K',k);
%dati dinamici simili => ogni cluster contribuisce con una sua libreria
idx_xy = idx_xy(:,2:end);
%si scrive idx_xy(:,2:end) per eliminare il primo elemento per eliminare il
%primo elemento dal cluster (cioè il più vicino) che rappresenta il
%campione stesso che viene utilizzato come centroid per l'applicazione di
%knnsearch

% run HySINDY on the undefined data with small (few data) stiff model
% initialize a structure library
II_Xi = [];
for ii = 1:floor(numOfPts)
    
    % stack all the "position" values in cluster ii into a vector
    X = [x(idx_xy(ii,:)') y(idx_xy(ii,:)')];
    % stack all the "velocity" values in cluster ii into a vector
    XP = [a_out_und(idx_xy(ii,:)',1) a_out_und(idx_xy(ii,:)',2)];
    
    
    % create function library
    polyorder_und = 0;  % search space up to i-th order polynomials
    polyorder_und_x = 1;  % search space up to i-th order polynomials
    polyorder_und_y = 0;  % search space up to i-th order polynomials
    usesine_und = 0;    % no trig functions
    laurent_und = 0;    % no laurent terms
    dyin_und = 0;        % just regular SINDy
    dyorder_und = 0;    % no derivatives in library
    Theta = poolDatady(X,2,polyorder_und,polyorder_und_x,polyorder_und_y,usesine_und, laurent_und...
        , dyin_und, dyorder_und);
    m = size(Theta,2);
    ThetaN = zeros(size(Theta,1),m);
    for rr = 1:m
%         normTheta(rr) = mean(Theta(:,rr));
        ThetaN(:,rr) = Theta(:,rr)/mean(Theta(:,rr));
    end
    Thetalib.Theta = ThetaN;
    Thetalib.normTheta = 1;
    Thetalib.dx = XP;
    Thetalib.polyorder = polyorder_und;
    Thetalib.polyorderx = polyorder_und_x;
    Thetalib.polyordery = polyorder_und_y;
    Thetalib.polyorder = polyorder_und;    
    Thetalib.usesine = usesine_und;
    
    lambdavals.numlambda = 20;
    lambdavals.lambdastart = -5;
    lambdavals.lambdaend = 2;
    
    Xistruct =multiD_Lambda(Thetalib,lambdavals);
    Xicomb = Xistruct.Xicomb;
    
    % add new models to library and find the indices for this set
    [II_Xitemp, II_ind] = build_Xi_Library(Xicomb, II_Xi);
    II_Xi = II_Xitemp;

    Cluster{ii}.centroidx =  mean(X(:,1));
    Cluster{ii}.centroidy = mean(X(:,2));
     % store results by cluster number
    Cluster{ii}.Xistruct = Xistruct;
    Cluster{ii}.II_ind = II_ind;
    Cluster{ii}.X = X;
    Cluster{ii}.XP = XP;
end

%find switch point
clear x2 y2 x y II_Xi
%x2 & y2 can be training or validation data
x2= y_out_und(1:end-1,1); % height
y2= y_out_und(1:end-1,2); % velocity

numOfPts = length(x2); % number of points in the data
% k = 16; % number of points in a cluster

val_ntimes = floor(0.2/dt); %length of comparision not including ICs
data2 = [x2 y2];

Cluster_xswitch_max = [];

options_test = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_diverge);
    
for ii = 1:length(Cluster)
    Xicomb = Cluster{ii}.Xistruct.Xicomb;
    
    idx_xy2 = knnsearch([x2 y2],...
        [Cluster{ii}.centroidx Cluster{ii}.centroidy],'K',k);
    idx_xy2 = idx_xy2(:,2:end);

    [xA, x0clust, tvectest] = buildinit_fromcluster(val_ntimes, 1, idx_xy2, data2, dt, k);
    
    val.tA =tvectest;
    val.xA =xA;
    val.x0 =x0clust;
    val.options = options_test;
    
    clear IC aic_c aic_cv xswitch_max xswitch_min nodynam_any xswitch_avg xswitch 
    clear nodynam abserror_bf_switch mean_e numterms idx_xy2
    for kk = 1:length(Xicomb)
        disp(['Cluster num:' num2str(ii) ' of' num2str(length(Cluster))])
        disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
        Xi = Xicomb{kk}
        
        clear  savetB savexB
        [savetB, savexB] = validation_sims(Xi, Thetalib, val, plottag);
        
        for ll = 1:length(xA)
            clear xswitch nodynam1 abserror abserror_avg1 RMSE xAcomp xBcomp
            [tlength, nterms] = size(savexB{ll});
            [xswitch, abserror, abserror_avg1, RMSE] = calc_tlength(xA{ll}, savexB{ll},val);
            xswitch(ll,:) = xswitch;
            abserror_bf_switch(ll,:) = abserror_avg1; %abs error befor that the switching takes place
         
        end
        abserror = [abserror_bf_switch(:,2); abserror_bf_switch(:,1)];
        IC{kk} = ICcalculations(abserror, nnz(Xi),size(x0clust,2));
        aic_c(kk,:) = IC{kk}.aic_c;
        xswitch_max(kk,:) = max(xswitch);
        xswitch_min(kk,:) = min(xswitch);
        xswitch_avg(kk,:) = mean(xswitch);
        mean_e(kk,:) = mean(abserror_bf_switch);
    
        
    end
    
    minIC = min(aic_c);
    
    Cluster_xswitch_max = [Cluster_xswitch_max ; mean(xswitch_max)];
    
%     Cluster{ii}.xswitch_max = min(xswitch_max);
%     Cluster{ii}.abserror_bf = mean_e;
%     Cluster{ii}.IC = IC;
%     Cluster{ii}.aic_c =aic_c; 

end
clear Cluster Theta normTheta ThetaN Thetalib lambdavals Xistruct Xicomb
clear II_ind val Xi x2 y2 x0clust xA X XP savetB savexB

%% Useful in some case to avoid the identification of discontinuity in the undefined region boundary layer
idx=find(Cluster_xswitch_max ~= max(y_out_und(:,1)) | Cluster_xswitch_max ~= min(y_out_und(:,1)));
Cluster_xswitch_max_correction = Cluster_xswitch_max(idx); clear idx
%% Use kmeans to divide the correct cluster of possible switch point from the other
% plot(Cluster_xswitch_max,'*') %pre-analysis of switch point detection
% The xswitch data set(Cluster_xswitch_max) is divided into multiple
% clusters to distinguish the clusters of the true xswitch from the false
% ones detected
optskmeans = statset('Display','final');
[idxkmeans,Ckmeans] = kmeans(Cluster_xswitch_max_correction,7,'Distance','cityblock',...
    'Replicates',5,'Options',optskmeans);
figure;
plot(Cluster_xswitch_max_correction(idxkmeans==1,1),'r.','MarkerSize',11)
hold on
plot(Cluster_xswitch_max_correction(idxkmeans==2,1),'b.','MarkerSize',11)
plot(Cluster_xswitch_max_correction(idxkmeans==3,1),'g.','MarkerSize',11)
plot(Cluster_xswitch_max_correction(idxkmeans==4,1),'k.','MarkerSize',11)
plot(Cluster_xswitch_max_correction(idxkmeans==5,1),'y.','MarkerSize',11)
plot(Cluster_xswitch_max_correction(idxkmeans==6,1),'m.','MarkerSize',11)
plot(Cluster_xswitch_max_correction(idxkmeans==7,1),'c.','MarkerSize',11)
plot(Ckmeans(:,1),'kx','MarkerSize',11,'LineWidth',1.5)

%%
% We find switch position by use xswitch_max because the algorithm
% findchangepts find the switch prematurely. Then we use mean to add
% statistic
num_idxkmeans = unique(idxkmeans);
histckmeans = [num_idxkmeans,histc(idxkmeans(:),num_idxkmeans)];
[num_idxkmeans, switch_idxkeams] = sort(histckmeans(:,2),'descend');
switch_idxkmeans = switch_idxkeams(1,1);
xswitch_median = median(Cluster_xswitch_max_correction(idxkmeans==switch_idxkmeans,1)); %switch position
% xswitch_mean = mean(Cluster_xswitch_max(idxkmeans==switch_idxkmeans,1)); %switch position

% We find switch position by use xswitch_max because the algorithm
% findchangepts find the switch prematurely. Then we use mean to add
% statistic
    
switch_index_spring = find(y_out(:,1)<=xswitch_median);
t_out_spring = t_out(switch_index_spring);
y_out_spring = y_out(switch_index_spring,:);
a_out_spring = a_out(switch_index_spring,:);
switch_index_flight = find(y_out(:,1)>xswitch_median);
t_out_flight = t_out(switch_index_flight);
y_out_flight = y_out(switch_index_flight,:);
a_out_flight = a_out(switch_index_flight,:);

%% spring mainfold
% initialize a structure library
ind_all_spring = []; II_Xi_spring = [];

% stack all the "position" values in cluster ii into a vector
[~,indtspring] = sort(t_out_spring);
indtspring = [indtspring(1:20);]% indtspring(40:60); indtspring(80:100)]; %modification 11/11/2021 to avoid the interference of x^2 function
y_out_spring = y_out_spring(indtspring,:);%we preserve the time sequence
X_spring = y_out_spring;
% stack all the "velocity" values in cluster ii into a vector
a_out_spring = a_out_spring(indtspring,:);%we preserve the time sequence
XP_spring = a_out_spring;
   
% create function library
% SINDy
    
% create function library
polyorder_spring = 2;  % search space up to 2nd order polynomials
polyorder_spring_x = 1;  % search space up to i-nd order polynomials
polyorder_spring_y = 1;  % search space up to i-nd order polynomials
usesine_spring = 0;    % no trig functions
laurent_spring = 0;    % no laurent terms
dyin_spring = 0;        % just regular SINDy
dyorder_spring = 0;    % no derivatives in library
Theta_spring = poolDatady(X_spring,2,polyorder_spring,polyorder_spring_x,polyorder_spring_y,usesine_spring...
    , laurent_spring, dyin_spring, dyorder_spring);
m = size(Theta_spring,2);
ThetaN_spring = zeros(size(Theta_spring,1),m);
for rr = 1:m
%     normTheta(rr) = mean(Theta_spring(:,rr));
    ThetaN_spring(:,rr) = Theta_spring(:,rr)/mean(Theta_spring(:,rr));
end
Thetalib_spring.Theta = Theta_spring;
Thetalib_spring.normTheta = 1;
Thetalib_spring.dx = XP_spring;
Thetalib_spring.polyorder = polyorder_spring;
Thetalib_spring.polyorderx = polyorder_spring_x;
Thetalib_spring.polyordery = polyorder_spring_y;
Thetalib_spring.usesine = usesine_spring;

lambdavals_spring.numlambda = 20;
lambdavals_spring.lambdastart = -2;
lambdavals_spring.lambdaend = 2;

Xistruct_spring = multiD_Lambda(Thetalib_spring,lambdavals_spring);
Xicomb_spring = Xistruct_spring.Xicomb;  

% add new models to library and find the indices for this set
[II_Xitemp_spring, II_ind_spring] = build_Xi_Library(Xicomb_spring, II_Xi_spring);
II_Xi_spring = II_Xitemp_spring;
ind_all_spring = [ind_all_spring; II_ind_spring];

% store results
Cluster_spring.Xistruct = Xistruct_spring;
Cluster_spring.II_ind = II_ind_spring;
Cluster_spring.X = X_spring;
Cluster_spring.XP = XP_spring;


%% flight manifold
% initialize a structure library
ind_all_flight = []; II_Xi_flight = [];

% stack all the "position" values in cluster ii into a vector
[~,indtflight] = sort(t_out_flight);
y_out_flight = y_out_flight(indtflight,:);%we preserve the time sequence
X_flight = y_out_flight;
% stack all the "velocity" values in cluster ii into a vector
a_out_flight = a_out_flight(indtflight,:);%we preserve the time sequence
XP_flight = a_out_flight;


% create function library
% SINDy
    
% create function library
polyorder_flight = 2;  % search space up to 2nd order polynomials
polyorder_flight_x = 1;  % search space up to i-nd order polynomials
polyorder_flight_y = 1;  % search space up to i-nd order polynomials
usesine_flight = 0;    % no trig functions
laurent_flight = 0;    % no laurent terms
dyin_flight = 0;        % just regular SINDy
dyorder_flight = 0;    % no derivatives in library
Theta_flight = poolDatady(X_flight,2,polyorder_flight,polyorder_flight_x,polyorder_flight_y,usesine_flight...
    , laurent_flight, dyin_flight, dyorder_flight);
m = size(Theta_flight,2);
ThetaN_flight = zeros(size(Theta_flight,1),m);
for rr = 1:m
%     normTheta(rr) = mean(Theta_flight(:,rr));
    ThetaN_flight(:,rr) = Theta_flight(:,rr)/mean(Theta_flight(:,rr));
end
Thetalib_flight.Theta = Theta_flight;
Thetalib_flight.normTheta = 1;
Thetalib_flight.dx = XP_flight;
Thetalib_flight.polyorder = polyorder_flight;
Thetalib_flight.polyorderx = polyorder_flight_x;
Thetalib_flight.polyordery = polyorder_flight_y;
Thetalib_flight.usesine = usesine_flight;

lambdavals_flight.numlambda = 20;
lambdavals_flight.lambdastart = -2;
lambdavals_flight.lambdaend = 2;

Xistruct_flight =multiD_Lambda(Thetalib_flight,lambdavals_flight);
Xicomb_flight = Xistruct_flight.Xicomb;  

% add new models to library and find the indices for this set
[II_Xitemp_flight, II_ind_flight] = build_Xi_Library(Xicomb_flight, II_Xi_flight);
II_Xi_flight = II_Xitemp_flight;
ind_all_flight = [ind_all_flight; II_ind_flight];

% store results
Cluster_flight.Xistruct = Xistruct_flight;
Cluster_flight.II_ind = II_ind_flight;
Cluster_flight.X = X_flight;
Cluster_flight.XP = XP_flight;

tEnd = toc(tStart);


%%
dateformatout = 'mmddyyyy';
% save([datestr(now, dateformatout) 'Hopping_data_eps_zero.mat'])
% save([datestr(now, dateformatout) 'Hopping_data_eps_1e-6.mat'])

%% define correct Xi for checking
% Xi_truespring_xdot = [xdot=1]
% Xi_truespring_xdotdot = [costant=11 ; x=-10 ; xdot=0]
% Xi_truehop(2,2) = -10;
% Xi_truefly_xdot = [xdot=1]
% Xi_truefly_xdotdot = [constant=-1]
%% plot validation
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

clear data yout aout tout % clear data to recreate the valition data set

%validation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.yinitvec = [0.8 -0.1];
%No clustered data
figure
[data]  = RUN_Hopper(p);
y_out = data.yout;
a_out = data.aout;
t_out = data.tout;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%training data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.yinitvec = [0.78 -0.1];
%No clustered data
figure
[data_training]  = RUN_Hopper(p);
y_out_training = data_training.yout;
a_out_training = data_training.aout;
t_out_training = data_training.tout;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot_index_switch = find(y_out(:,1) == xswitch_median);
y0_spring = [1  -0.9]; %for data set ...
y0_flight = [1   0.9]; %for data set ...
   
tspan_spring = 0;
tspan_flight = 0;
for mm=1:numel(t_out_spring)-71 %we cut 144 points to avoid that springing dynamic continue beyond its own manifold
    tspan_spring = [tspan_spring mm*dt];
end
for mm=1:numel(t_out_flight)-135 %we cut 268 points to avoid that flying dynamic continue beyond its own manifold
    tspan_flight = [tspan_flight mm*dt];
end
Xicomb = Cluster_spring.Xistruct.Xicomb;
for kk = 1%:length(Xicomb)
    disp('Cluster spring')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_spring = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_spring,x_mod_spring] = ode45(@(t,x)sparseGalerkin(t,x,Xi_spring,polyorder_spring,polyorder_spring_x,polyorder_spring_y,usesine_spring),tspan_spring,y0_spring,options);  % approximate
    a_mod_spring = diff(x_mod_spring(:,2));
%     a_mod_spring = sparseGalerkin(t_mod_spring,x_mod_spring',Xi_spring,polyorder_spring,polyorder_spring_x,polyorder_spring_y,usesine_spring);  % approximate
end

clear Xicomb
Xicomb = Cluster_flight.Xistruct.Xicomb;
for kk = 1%:length(Xicomb)
    disp('Cluster flight')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_flight = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_flight,x_mod_flight] = ode45(@(t,x)sparseGalerkin(t,x,Xi_flight,polyorder_flight,polyorder_flight_x,polyorder_flight_y,usesine_flight),tspan_flight,y0_flight,options);  % approximate
    a_mod_flight = diff(x_mod_flight(:,2));
%     a_mod_flight = sparseGalerkin(t_mod_flight,x_mod_flight',Xi_flight,polyorder_flight,polyorder_flight_x,polyorder_flight_y,usesine_flight);  % approximate
end

figure
plot(y_out(:,1),y_out(:,2),'k--','LineWidth',0.9) % measurement - training data
hold on
plot(y_out_training(:,1),y_out_training(:,2),'k-','LineWidth',0.9) % measurement - validation data
plot(x_mod_spring(:,1),x_mod_spring(:,2),'r','LineWidth',0.9) % approx model
plot(x_mod_flight(:,1),x_mod_flight(:,2),'b','LineWidth',0.9) % approx model
plot([xswitch_median xswitch_median], [min(y_out(:,2)) max(y_out(:,2))],'k','LineWidth',0.9)
xlabel('$x [ ]$','interpreter','latex')
ylabel('$\dot{x} [ ]$','interpreter','latex')
title('HySINDy - Phase space','interpreter','latex')
legend('validation data','training data','springing','flying','$x{switch}$','interpreter','latex')

figure
plot(t_mod_spring(2:end),a_mod_spring,'r','LineWidth',0.9) % approx model
hold on
t_mod_flight = t_mod_flight+t_mod_spring(end);
plot(t_mod_flight(2:end),a_mod_flight,'b','LineWidth',0.9) % approx model
xlabel('$t [ ]$','interpreter','latex')
ylabel('$\stackrel{..}{x} [ ]$','interpreter','latex');
title('HySINDy - Acceleration','interpreter','latex')
legend('spring','flight','interpreter','latex')

figure
plot(t_mod_spring,x_mod_spring(:,1),'r','LineWidth',0.9) % approx model
hold on
plot(t_mod_flight,x_mod_flight(:,1),'b','LineWidth',0.9) % approx model
xlabel('$t [ ]$','interpreter','latex')
ylabel('$x [ ]$','interpreter','latex');
title('HySINDy - Position','interpreter','latex')
legend('spring','flight','interpreter','latex')


figure
plot(t_mod_spring,x_mod_spring(:,2),'r','LineWidth',0.9) % approx model
hold on
plot(t_mod_flight,x_mod_flight(:,2),'b','LineWidth',0.9) % approx model
xlabel('$t [ ]$','interpreter','latex')
ylabel('$\dot{x} [ ]$','interpreter','latex');
title('HySINDy - velocity','interpreter','latex')
legend('spring','flight','interpreter','latex')










