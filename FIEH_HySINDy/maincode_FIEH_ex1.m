%% main code - FIEH - supervised learning 
% This script is autonomous.
% Run this code to train and validate the model from data.

%Study of electromagnetic vibration energy harvesting with free/impact
%motion for low frequency operation

%% control parameters
% Shown below the parameters on which to act to adapt and improve the
% performance of the code to the various example cases:
% cluster dimension k
% val_ntimes
% thick_unddata
% switching and model library
%%

close all
clear all
clc

plottag =1;
addpath('FIEH'); addpath('SINDy'); addpath('validation');

%% training data generation
M = 0.03; %[kg] mass
Ca = 0.02; %[Ns/m] air damping coefficient
C = 1; %[Ns/m] spring viscous damping coefficient
m = 50; %[Am^2] magnetic moment
mu0 = 4*pi*1e-7; %[N/A^2] vacum permeability
D = 0.021; %[m] coil mean diameter
L = 0.05; %[m] coil length
N = 75; %[] total number of turns
R = 0.948; %[Omega] coil resistance
kappa = 196; %[N/m] FIEH stiffness
stroke = 0.03; %[m] stroke
mu = 0.2; %friction coefficient
g = 9.81; %[m/s^2] gravity acceleration
Yi = [0.01; 0.05] ; %0.01; %[m] exitaction amplitude
omega = 2*pi*6; %[Hz] exitaction frequency

thick_unddata = 0.005; % range undefined data position

eps = 0; % noise level

rng(1)
% In some situations, setting the seed alone will not guarantee the same
% results. This is because the generator that the random number functions
% draw from might be different than you expect when your code executes.
% For long-term repeatability, specify the seed and the generator type together.
% For example, the following code sets the seed to 1 and the generator to Mersenne Twister.
% rng(1,'twister');

% variables are: time, t, vertical position, y, vertical velocity, yp

tend = 11;% 12 time for training; 17 time for predict
dt = 0.001;%0.0001;
tspan = 0:dt:tend;

%initial condition
yinit = [0.01 0.02 2*pi*6 %for training 13/11/2021
    %0.02 0.025 2*pi*6 %for training --> this i.c. is not use for stable
    %data training analysis, cause it doesn't make sense.
    ];

% yinit = [0.01 0.02 2*pi*6]; %for training - old data
% yinit = [0.02 0.025 2*pi*6]; %for validation - old data

% yoverFcn = @(t, y) events_FIEH(t,y,stroke,Yi);
% options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', yoverFcn);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% p.options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_FIEH);
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

% tend = tspan(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_out time vector
% y_out contain both vertical position and velocity
% ydot_out contain both vertical acceleration and velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation model:
% M*ypp = -(Ca+Ci+Ce)*(ydot-Yi*omega*cos(omega*t))...
%-ke*(abs(y-Yi*sin(omega*t))-stroke/2)...
%-mu*M*g*(ydot-Yi*omega*cos(omega*t))/(abs(y-Yi*sin(omega*t)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_out = [];
y_out = [];
z_out = [];
zdot_out = [];
y_ext = [];
ydot_ext = [];
for ii=1:size(yinit,1)
    for jj=1:size(Yi,1)
        [t_out_ode45,y_out_ode45] = ode45(@(t,y) fieh(t,y,Ca,C,M,Yi(jj,:),omega,kappa,stroke,mu,g,N,D,mu0,m,R,L),tspan,yinit(ii,:));
        % a = fieh(t_out,y_out',Ca,C,M,Yi,omega*ones(1,numel(y_out(:,3))),kappa,stroke,mu,g,N,D,mu0,m,R,L)';

        %
        z_out_ode45 = [(y_out_ode45(:,1)-Yi(jj,:)*sin(y_out_ode45(:,3))) , (y_out_ode45(:,2)-Yi(jj,:)*omega*cos(y_out_ode45(:,3)))]; %augmented data state
        zdot_out_ode45 = diff(z_out_ode45,1,1); %Compute the first-order difference
        t_out_ode45 = t_out_ode45(1:end-1,:);
        z_out_ode45 = z_out_ode45(1:end-1,:);
        y_out_ode45 = y_out_ode45(1:end-1,:);
        
        % store only stable dynamic data
        clear ind_stab
        t_out_ode45 = t_out_ode45(t_out_ode45>=4);
        ind_stab = size(z_out_ode45,1)-size(t_out_ode45,1)+1;
        z_out_ode45 = z_out_ode45(ind_stab:end,:);
        zdot_out_ode45 = zdot_out_ode45(ind_stab:end,:);
        y_ext_ode45 = Yi(jj,:)*sin(y_out_ode45(ind_stab:end,3));
        ydot_ext_ode45 = Yi(jj,:)*omega*cos(y_out_ode45(ind_stab:end,3));
        
        % t_out = t_out(3805:3971);
        % z_out = z_out(3805:3971,:);
        % zdot_out = zdot_out(3805:3971,:);
        % y_ext = Yi*sin(y_out(3805:3971,3));
        % ydot_ext = Yi*omega*cos(y_out(3805:3971,3));


        t_out = [t_out; t_out_ode45];
        y_out = [y_out; y_out_ode45];
        z_out = [z_out; z_out_ode45];
        zdot_out = [zdot_out; zdot_out_ode45];
        y_ext = [y_ext; y_ext_ode45];
        ydot_ext = [ydot_ext; ydot_ext_ode45];        
    end
end

figure(1)
plot(t_out,z_out(:,1),'b--')
hold on
plot(t_out,z_out(:,2),'r')
xlabel('$t$ [ ]','Interpreter','latex')
ylabel('$\dot{x}$ [ ]','Interpreter','latex')
title('FIEH System')

% % %% flux rate of change
% dphidz = []; %[W/m] flux rate of change
% for jj = 1:size(y_out,1)
%     yy = y_out(jj,:);
%     dphi_star = 0;
%     for n = -N/2:N/2
%         dphi_star = dphi_star + ((yy(1)-Yi*sin(yy(3)))+n*(L/(N-1)))*...
%             (((2*((yy(1)-Yi*sin(yy(3)))+n*(L/(N-1)))^2)/(((D^2)/4+((yy(1)-Yi*sin(yy(3)))+n*(L/(N-1)))^2)^(5/2)))...
%             -(3/(((D^2)/4+((yy(1)-Yi*sin(yy(3)))+n*(L/(N-1)))^2)^(3/2))));
%     end
%     dphidz = [dphidz , (mu0*m/2*dphi_star)];
% end
% 
% figure(3);clf;
% plot(y_out(:,1)-Yi*sin(y_out(:,3)),dphidz)

PercUndData = 100*2*thick_unddata/(max(z_out(:,1))-min(z_out(:,1)));

%%
%plot - check phase space
figure
plot(z_out(:,1),z_out(:,2),'k')
xlabel('$x$ [ ]','Interpreter','latex')
ylabel('$\dot{x}$ [ ]','Interpreter','latex')
title('phase space data')

% % method 1
% % index undefined data
% indx_und = find(z_out(:,1)<=(min(z_out(:,1))+thick_unddata) | z_out(:,1)>=(max(z_out(:,1))-thick_unddata));
% % indx_und = find(abs(z_out(:,1))<=(0.015+thick_unddata/2) & abs(z_out(:,1))<=(-0.015-thick_unddata/2));
% % index defined data
% indx_free = find(z_out(:,1)>(min(z_out(:,1))+thick_unddata) & z_out(:,1)<(max(z_out(:,1))-thick_unddata));

% method 2
% index undefined data
indx_und = find(abs(z_out(:,1))<=(0.015+thick_unddata/2) & abs(z_out(:,1))>=(0.015-thick_unddata/2));
% index defined data
indx_free = find(z_out(:,1)>(-0.015+thick_unddata/2) & z_out(:,1)<(0.015-thick_unddata/2));
% index defined data
indx_impact = find(z_out(:,1)<(-0.015-thick_unddata/2) | z_out(:,1)>(0.015+thick_unddata/2));

%% Clustered data
z_out_und = z_out(indx_und,:);
zdot_out_und = zdot_out(indx_und,:);
t_out_und = t_out(indx_und);

z_out_free = z_out(indx_free,:);
zdot_out_free = zdot_out(indx_free,:);
t_out_free = t_out(indx_free);

z_out_impact = z_out(indx_impact,:);
zdot_out_impact = zdot_out(indx_impact,:);
t_out_impact = t_out(indx_impact);

%% No clustered data
% data.z_out = z_out;
% data.zdot_out = zdot_out;
% data.t_out = t_out;

%save external "forcing"
% data.y_ext = y_ext;
% data.ydot_ext = ydot_ext;

%% check - plot
%position
figure
plot(t_out_und,z_out_und(:,1),'r.')
hold on
plot(t_out_free,z_out_free(:,1),'b.')
plot(t_out_impact,z_out_impact(:,1),'g.')
plot(t_out,y_ext,'k--')
% xlimits = [];
% ylimits = [];
% axis([xlimits ylimits])
l = get(gca, 'Children');
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$x$ [ ]','Interpreter','latex')
legend('und.area','free','external forcing')
title('time-position data')
set(gca,'FontSize',9)
set(l, 'linewidth', 1.5)
%velocity
figure
plot(t_out_und,z_out_und(:,2),'r.')
hold on
plot(t_out_free,z_out_free(:,2),'b.')
plot(t_out_impact,z_out_impact(:,2),'g.')
plot(t_out,ydot_ext,'k--')
l = get(gca, 'Children');
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\dot{x}$ [ ]','Interpreter','latex')
legend('und.area','free','external forcing')
title('time-velocity data')
set(gca,'FontSize',9)
set(l, 'linewidth', 1.5)
%acceleration
figure
plot(t_out_und,zdot_out_und(:,2),'r.')
hold on
plot(t_out_free,zdot_out_free(:,2),'b.')
plot(t_out_impact,zdot_out_impact(:,2),'g.')
l = get(gca, 'Children');
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\stackrel{..}{x}$ [ ]','Interpreter','latex')
legend('und.area','free')
title('time-acceleration data')
set(gca,'FontSize',9)
set(l, 'linewidth', 1.5)
%phase space
figure
plot(z_out_und(:,1),zdot_out_und(:,1),'r.')
hold on
plot(z_out_free(:,1),zdot_out_free(:,1),'b.')
plot(z_out_impact(:,1),zdot_out_impact(:,1),'g.')
% xlimits = [];
% ylimits = [];
% axis([xlimits ylimits])
l = get(gca, 'Children');
xlabel('$x$ [ ]','Interpreter','latex')
ylabel('$\dot{x}$ [ ]','Interpreter','latex')
legend('und.area','free')
title('phase space data')
set(gca,'FontSize',9)
set(l, 'linewidth', 1.5)

%% undefined data
tStart = tic;

plottag =0;
dimMan = 2; % dimension of the Manifold

% Run CCM on y_out
clear x y
x= z_out_und(1:end-1,1); % height
y= z_out_und(1:end-1,2); % velocity

numOfPts = length(x); % number of points in the data
%k = dimMan+2;
k = 7; % number of points in a cluster

idx_xy = knnsearch([x/max(x) y/max(y)],[x/max(x) y/max(y)],'K',k);
%similar dynamic data => each cluster contributes with its own library
idx_xy = idx_xy(:,2:end);

% run HySINDY on the undefined data with small (few data) stiff model
% initialize a structure library
II_Xi = [];
for ii = 1:floor(numOfPts)
     
%     %approach 1
%     % stack all the "position" values in cluster ii into a vector
%     X = [x(idx_xy(ii,:)') y(idx_xy(ii,:)')];
%     % stack all the "velocity" values in cluster ii into a vector
%     XP = [zdot_out_und(idx_xy(ii,:)',1) zdot_out_und(idx_xy(ii,:)',2)];
    
    %approach 2
    % stack all the velocity values in cluster ii into a vector
    X = y(idx_xy(ii,:)');
    % stack all the acceleration values in cluster ii into a vector
    XP = zdot_out_und(idx_xy(ii,:)',2);
    
    
    % create function library
    polyorder_und = 1;  % search space up to i-th order polynomials
    polyorder_und_x = 1;  % search space up to i-th order polynomials
    polyorder_und_y = 0;  % search space up to i-th order polynomials
    usesine_und = 0;    % no trig functions
    laurent_und = 0;    % no laurent terms
    dyin_und = 0;        % just regular SINDy
    dyorder_und = 0;    % no derivatives in library
    Theta = poolDatady(X,1,polyorder_und,polyorder_und_x,polyorder_und_y,usesine_und, laurent_und...
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
    lambdavals.lambdastart = -4;
    lambdavals.lambdaend = 0;
    
    Xistruct =multiD_Lambda(Thetalib,lambdavals);
    Xicomb = Xistruct.Xicomb;
    
    % add new models to library and find the indices for this set
    [II_Xitemp, II_ind] = build_Xi_Library(Xicomb, II_Xi);
    II_Xi = II_Xitemp;

%     % approach 1
%     Cluster{ii}.centroidx =  mean(X(:,1));
%     Cluster{ii}.centroidy = mean(X(:,2));

    % approach 2
    Cluster{ii}.centroidx =  mean(X(:,1));
    
     % store results by cluster number
    Cluster{ii}.Xistruct = Xistruct;
    Cluster{ii}.II_ind = II_ind;
    Cluster{ii}.X = X;
    Cluster{ii}.XP = XP;
end

%find switch point
clear x2 y2 x y II_Xi
%x2 & y2 can be training or validation data
x2= z_out_und(1:end-1,1); % height
y2= z_out_und(1:end-1,2); % velocity

numOfPts = length(x2); % number of points in the data
% k = 16; % number of points in a cluster

val_ntimes = floor(0.01/dt); %length of comparision not including ICs
% %approach 1
% data2 = [x2 y2];

%approach 2
data2 = y2;

Cluster_zswitch_max = [];

options_test = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_diverge);
    
for ii = 1:length(Cluster)
    Xicomb = Cluster{ii}.Xistruct.Xicomb;
    
%     %approach 1
%     idx_xy2 = knnsearch([x2 y2],...
%         [Cluster{ii}.centroidx Cluster{ii}.centroidy],'K',k);
%     idx_xy2 = idx_xy2(:,2:end);

    %approach 2
    idx_xy2 = knnsearch(y2,...
        [Cluster{ii}.centroidx],'K',k);
    idx_xy2 = idx_xy2(:,2:end);

    [xA, x0clust, tvectest] = buildinit_fromcluster(val_ntimes, 1, idx_xy2, data2, dt, k);
    
    val.tA =tvectest;
    val.xA =xA;
    val.x0 =x0clust;
    val.options = options_test;
    
    clear IC aic_c aic_cv xswitch_max xswitch_min nodynam_any xswitch_avg xswitch 
    clear nodynam abserror_bf_switch mean_e numterms idx_xy2
    zswitch_avg = [];
    for kk = 1:length(Xicomb)
        disp(['Cluster num:' num2str(ii) ' of' num2str(length(Cluster))])
        disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
        Xi = Xicomb{kk}
        
        clear  savetB savexB
        [savetB, savexB] = validation_sims(Xi, Thetalib, val, plottag);
        
        zswitch = [];
        abserror_bf_switch = [];
        for ll = 1:length(xA)
            clear xswitch nodynam1 abserror abserror_avg1 RMSE xAcomp xBcomp
            [tlength, nterms] = size(savexB{ll});
            [z_switch, abserror, abserror_avg1, RMSE] = calc_tlength(xA{ll}, savexB{ll}, val);
%             [z_switch, abserror] = calc_tlength_for_switching_identification(xA{ll}, savexB{ll});
%             zswitch(ll,:) = z_switch;
            zswitch = [zswitch; z_switch];
%             abserror_bf_switch(ll,:) = abserror_avg1; %abs error befor that the switching takes place
            abserror_bf_switch  = [abserror_bf_switch; abserror_avg1];       
        end
%         %approach 1
%         abserror = [abserror_bf_switch(:,2); abserror_bf_switch(:,1)];
        
        %approach 2
        abserror = abserror_bf_switch(:,1);
            
        IC{kk} = ICcalculations(abserror, nnz(Xi),size(x0clust,2));
        aic_c(kk,:) = IC{kk}.aic_c;
%         zswitch_max(kk,:) = max(zswitch);
%         zswitch_min(kk,:) = min(zswitch);
%         zswitch_avg(kk,:) = median(zswitch);
%         mean_e(kk,:) = mean(abserror_bf_switch);
        zswitch_avg = [zswitch_avg; median(zswitch)];
    
        
    end
    
    minIC = min(aic_c);
    
    Cluster_zswitch_max = [Cluster_zswitch_max ; median(zswitch_avg)];
    
%     Cluster{ii}.xswitch_max = min(xswitch_max);
%     Cluster{ii}.abserror_bf = mean_e;
%     Cluster{ii}.IC = IC;
%     Cluster{ii}.aic_c =aic_c; 

end
clear Cluster Theta normTheta ThetaN Thetalib lambdavals Xistruct Xicomb
clear II_ind val Xi x2 y2 x0clust xA X XP savetB savexB

%% Use kmeans to divide the correct cluster of possible switch point from the other
% plot(Cluster_zswitch_max,'*') %pre-analysis of switch point detection
% The zswitch data set(Cluster_zswitch_max) is divided into multiple
% clusters to distinguish the clusters of the true zswitch from the false
% ones detected
optskmeans = statset('Display','final');
[idxkmeans,Ckmeans] = kmeans(Cluster_zswitch_max,13,'Distance','cityblock',...
    'Replicates',17,'Options',optskmeans);
figure;
plot(Cluster_zswitch_max(idxkmeans==1,1),'r.','MarkerSize',11)
hold on
plot(Cluster_zswitch_max(idxkmeans==2,1),'b.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==3,1),'g.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==4,1),'k.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==5,1),'y.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==6,1),'m.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==7,1),'c.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==8,1),'.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==9,1),'.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==10,1),'.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==11,1),'.','MarkerSize',11)
plot(Cluster_zswitch_max(idxkmeans==12,1),'.','MarkerSize',11)
% plot(Cluster_zswitch_max(idxkmeans==13,1),'.','MarkerSize',11)
% plot(Cluster_zswitch_max(idxkmeans==14,1),'.','MarkerSize',11)
% plot(Cluster_zswitch_max(idxkmeans==15,1),'.','MarkerSize',11)
% plot(Cluster_zswitch_max(idxkmeans==16,1),'.','MarkerSize',11)
% plot(Cluster_zswitch_max(idxkmeans==17,1),'.','MarkerSize',11)
% plot(Cluster_zswitch_max(idxkmeans==18,1),'.','MarkerSize',11)
% plot(Cluster_zswitch_max(idxkmeans==19,1),'.','MarkerSize',11)

plot(Ckmeans(:,1),'kx','MarkerSize',13,'LineWidth',1.5)

%%
% We find switch position by use xswitch_max because the algorithm
% findchangepts find the switch prematurely. Then we use mean to add
% statistic
num_idxkmeans = unique(idxkmeans);
histckmeans = [num_idxkmeans,histc(idxkmeans(:),num_idxkmeans)];
[num_idxkmeans, switch_idxkeams] = sort(histckmeans(:,2),'descend');
switch_idxkmeans1 = switch_idxkeams(1,1);
switch_idxkmeans2 = switch_idxkeams(2,1);
zswitch_median = sort([median(Cluster_zswitch_max(idxkmeans==switch_idxkmeans1,1)) , median(Cluster_zswitch_max(idxkmeans==switch_idxkmeans2,1))]);
%
ind_negative = find(z_out_und(:,1)<0);
[~,idxdown]=min(abs(z_out_und(ind_negative,2)-zswitch_median(:,1)));
idx_pos1 = ind_negative(idxdown);
zswitch_median_down = z_out_und(idx_pos1,1); %switch position 1
%
ind_positive = find(z_out_und(:,1)>0);
[~,idxtop]=min(abs(z_out_und(ind_positive,2)-zswitch_median(:,2)));
idx_pos2 = ind_positive(idxtop);
zswitch_median_top = z_out_und(idx_pos2,1); %switch position 2
% zswitch_median_down = zswitch_median(1,1); %switch position 1
% zswitch_median_top = zswitch_median(1,2); %switch position 2

% true switches
% zswitch_median_down = -0.015; %switch position 1
% zswitch_median_top = 0.015; %switch position 2

% We find switch position by use xswitch_max because the algorithm
% findchangepts find the switch prematurely. Then we use mean to add
% statistic
    
switch_index_impactDown = find(z_out(:,1)<=zswitch_median_down);
t_out_impactDown = t_out(switch_index_impactDown);
z_out_impactDown = z_out(switch_index_impactDown,:);
zdot_out_impactDown = zdot_out(switch_index_impactDown,:);
%
switch_index_impactTop = find(z_out(:,1)>=zswitch_median_top);
t_out_impactTop = t_out(switch_index_impactTop);
z_out_impactTop = z_out(switch_index_impactTop,:);
zdot_out_impactTop = zdot_out(switch_index_impactTop,:);
%
switch_index_free = find(z_out(:,1)>zswitch_median_down & z_out(:,1)<zswitch_median_top);
t_out_free = t_out(switch_index_free);
z_out_free = z_out(switch_index_free,:);
zdot_out_free = zdot_out(switch_index_free,:);

%% Down impact mainfold
% initialize a structure library
ind_all_impactDown = []; II_Xi_impactDown = [];

% stack all the "position" values in cluster ii into a vector
[~,indtimpactDown] = sort(t_out_impactDown);
z_out_impactDown = z_out_impactDown(indtimpactDown,:);%we preserve the time sequence
X_impactDown = z_out_impactDown;
% stack all the "velocity" values in cluster ii into a vector
zdot_out_impactDown = zdot_out_impactDown(indtimpactDown,:);%we preserve the time sequence
XP_impactDown = zdot_out_impactDown;
   
% create function library
% SINDy
    
% create function library
polyorder_impactDown = 2;  % search space up to 2nd order polynomials
polyorder_impactDown_x = 1;  % search space up to i-nd order polynomials
polyorder_impactDown_y = 1;  % search space up to i-nd order polynomials
usesine_impactDown = 0;    % no trig functions
laurent_impactDown = 0;    % no laurent terms
dyin_impactDown = 0;        % just regular SINDy
dyorder_impactDown = 0;    % no derivatives in library
Theta_impactDown = poolDatady(X_impactDown,2,polyorder_impactDown,polyorder_impactDown_x,polyorder_impactDown_y,usesine_impactDown...
    , laurent_impactDown, dyin_impactDown, dyorder_impactDown);
m = size(Theta_impactDown,2);
ThetaN_impactDown = zeros(size(Theta_impactDown,1),m);
for rr = 1:m
%     normTheta(rr) = mean(Theta_impact(:,rr));
    ThetaN_impactDown(:,rr) = Theta_impactDown(:,rr)/mean(Theta_impactDown(:,rr));
end
Thetalib_impactDown.Theta = Theta_impactDown;
Thetalib_impactDown.normTheta = 1;
Thetalib_impactDown.dx = XP_impactDown;
Thetalib_impactDown.polyorder = polyorder_impactDown;
Thetalib_impactDown.polyorderx = polyorder_impactDown_x;
Thetalib_impactDown.polyordery = polyorder_impactDown_y;
Thetalib_impactDown.usesine = usesine_impactDown;

lambdavals_impact.numlambda = 20;
lambdavals_impact.lambdastart = -6;
lambdavals_impact.lambdaend = 4;

Xistruct_impactDown = multiD_Lambda(Thetalib_impactDown,lambdavals_impact);
Xicomb_impactDown = Xistruct_impactDown.Xicomb;  

% add new models to library and find the indices for this set
[II_Xitemp_impactDown, II_ind_impactDown] = build_Xi_Library(Xicomb_impactDown, II_Xi_impactDown);
II_Xi_impactDown = II_Xitemp_impactDown;
ind_all_impactDown = [ind_all_impactDown; II_ind_impactDown];

% store results
Cluster_impactDown.Xistruct = Xistruct_impactDown;
Cluster_impactDown.II_ind = II_ind_impactDown;
Cluster_impactDown.X = X_impactDown;
Cluster_impactDown.XP = XP_impactDown;


%% Top impact mainfold
% initialize a structure library
ind_all_impactTop = []; II_Xi_impactTop = [];

% stack all the "position" values in cluster ii into a vector
[~,indtimpactTop] = sort(t_out_impactTop);
z_out_impactTop = z_out_impactTop(indtimpactTop,:);%we preserve the time sequence
X_impactTop = z_out_impactTop;
% stack all the "velocity" values in cluster ii into a vector
zdot_out_impactTop = zdot_out_impactTop(indtimpactTop,:);%we preserve the time sequence
XP_impactTop = zdot_out_impactTop;
   
% create function library
% SINDy
    
% create function library
polyorder_impactTop = 3;  % search space up to 2nd order polynomials
polyorder_impactTop_x = 1;  % search space up to i-nd order polynomials
polyorder_impactTop_y = 1;  % search space up to i-nd order polynomials
usesine_impactTop = 0;    % no trig functions
laurent_impactTop = 0;    % no laurent terms
dyin_impactTop = 0;        % just regular SINDy
dyorder_impactTop = 0;    % no derivatives in library
Theta_impactTop = poolDatady(X_impactTop,2,polyorder_impactTop,polyorder_impactTop_x,polyorder_impactTop_y,usesine_impactTop...
    , laurent_impactTop, dyin_impactTop, dyorder_impactTop);
m = size(Theta_impactTop,2);
ThetaN_impactTop = zeros(size(Theta_impactTop,1),m);
for rr = 1:m
%     normTheta(rr) = mean(Theta_impact(:,rr));
    ThetaN_impactTop(:,rr) = Theta_impactTop(:,rr)/mean(Theta_impactTop(:,rr));
end
Thetalib_impactTop.Theta = Theta_impactTop;
Thetalib_impactTop.normTheta = 1;
Thetalib_impactTop.dx = XP_impactTop;
Thetalib_impactTop.polyorder = polyorder_impactTop;
Thetalib_impactTop.polyorderx = polyorder_impactTop_x;
Thetalib_impactTop.polyordery = polyorder_impactTop_y;
Thetalib_impactTop.usesine = usesine_impactTop;

lambdavals_impact.numlambda = 20;
lambdavals_impact.lambdastart = -6;
lambdavals_impact.lambdaend = 4;

Xistruct_impactTop = multiD_Lambda(Thetalib_impactTop,lambdavals_impact);
Xicomb_impactTop = Xistruct_impactTop.Xicomb;  

% add new models to library and find the indices for this set
[II_Xitemp_impactTop, II_ind_impactTop] = build_Xi_Library(Xicomb_impactTop, II_Xi_impactTop);
II_Xi_impactTop = II_Xitemp_impactTop;
ind_all_impactTop = [ind_all_impactTop; II_ind_impactTop];

% store results
Cluster_impactTop.Xistruct = Xistruct_impactTop;
Cluster_impactTop.II_ind = II_ind_impactTop;
Cluster_impactTop.X = X_impactTop;
Cluster_impactTop.XP = XP_impactTop;


%% free manifold
% initialize a structure library
ind_all_free = []; II_Xi_free = [];

% stack all the "position" values in cluster ii into a vector
[~,indtfree] = sort(t_out_free);
z_out_free = z_out_free(indtfree,:);%we preserve the time sequence
X_free = z_out_free;
% stack all the "velocity" values in cluster ii into a vector
zdot_out_free = zdot_out_free(indtfree,:);%we preserve the time sequence
XP_free = zdot_out_free;


% create function library
% SINDy
    
% create function library
polyorder_free = 2;  % search space up to 2nd order polynomials
polyorder_free_x = 1;  % search space up to i-nd order polynomials
polyorder_free_y = 1;  % search space up to i-nd order polynomials
usesine_free = 0;    % no trig functions
laurent_free = 0;    % no laurent terms
dyin_free = 0;        % just regular SINDy
dyorder_free = 0;    % no derivatives in library
Theta_free = poolDatady(X_free,2,polyorder_free,polyorder_free_x,polyorder_free_y,usesine_free...
    , laurent_free, dyin_free, dyorder_free);
m = size(Theta_free,2);
ThetaN_free = zeros(size(Theta_free,1),m);
for rr = 1:m
%     normTheta(rr) = mean(Theta_free(:,rr));
    ThetaN_free(:,rr) = Theta_free(:,rr)/mean(Theta_free(:,rr));
end
Thetalib_free.Theta = Theta_free;
Thetalib_free.normTheta = 1;
Thetalib_free.dx = XP_free;
Thetalib_free.polyorder = polyorder_free;
Thetalib_free.polyorderx = polyorder_free_x;
Thetalib_free.polyordery = polyorder_free_y;
Thetalib_free.usesine = usesine_free;

lambdavals_free.numlambda = 20;
lambdavals_free.lambdastart = -8;
lambdavals_free.lambdaend = 4;

Xistruct_free =multiD_Lambda(Thetalib_free,lambdavals_free);
Xicomb_free = Xistruct_free.Xicomb;  

% add new models to library and find the indices for this set
[II_Xitemp_free, II_ind_free] = build_Xi_Library(Xicomb_free, II_Xi_free);
II_Xi_free = II_Xitemp_free;
ind_all_free = [ind_all_free; II_ind_free];

% store results
Cluster_free.Xistruct = Xistruct_free;
Cluster_free.II_ind = II_ind_free;
Cluster_free.X = X_free;
Cluster_free.XP = XP_free;
%%
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
%% plot
% options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_FIEH);
clear Yi
%% validation data generation
M = 0.03; %[kg] mass --> nell'articolo è 0.03 %TEST NICO 0.09
Ca = 0.02; %[Ns/m] air damping coefficient
C = 1; %[Ns/m] spring viscous damping coefficient
m = 50; %[Am^2] magnetic moment
mu0 = 4*pi*1e-7; %[N/A^2] vacum permeability
D = 0.021; %[m] coil mean diameter
L = 0.05; %[m] coil length
N = 75; %[] total number of turns %con 52 si ottiene circa la stessa variazione di flusso magnetico
R = 0.948; %[Omega] coil resistance
kappa = 196; %[N/m] FIEH stiffness --> nell'articolo è 196 %TEST NICO 70
stroke = 0.03; %[m] stroke
mu = 0.2; %friction coefficient --> nell'articolo è 0.2 %TEST NICO 0
g = 9.81; %[m/s^2] gravity acceleration
Yi = [0.03] ; %0.03; %[m] exitaction amplitude
omega = 2*pi*6; %[Hz] exitaction frequency

thick_unddata = 0.005;%0.013; % range undefined data position

eps = 0; % noise level

rng(1)
% In some situations, setting the seed alone will not guarantee the same
% results. This is because the generator that the random number functions
% draw from might be different than you expect when your code executes.
% For long-term repeatability, specify the seed and the generator type together.
% For example, the following code sets the seed to 1 and the generator to Mersenne Twister.
% rng(1,'twister');

% variables are: time, t, vertical position, y, vertical velocity, yp

tend = 9.17;
dt = 0.001;%0.0001; %TEST NICO 0.0001
tspan = 0:dt:tend;

%initial condition
yinit = [0.01 0.02 2*pi*6 %for training 13/11/2021
    %0.02 0.025 2*pi*6 %for training --> this i.c. is not use for stable
    %data training analysis, cause it doesn't make sense.
    ];

% yinit = [0.01 0.02 2*pi*6]; %for training - old data
% yinit = [0.02 0.025 2*pi*6]; %for validation - old data

% yoverFcn = @(t, y) events_FIEH(t,y,stroke,Yi);
% options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', yoverFcn);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% p.options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_FIEH);
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

% tend = tspan(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_out time vector
% y_out contain both vertical position and velocity
% ydot_out contain both vertical acceleration and velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation model:
% M*ypp = -(Ca+Ci+Ce)*(ydot-Yi*omega*cos(omega*t))...
%-ke*(abs(y-Yi*sin(omega*t))-stroke/2)...
%-mu*M*g*(ydot-Yi*omega*cos(omega*t))/(abs(y-Yi*sin(omega*t)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_out = [];
y_out = [];
z_out = [];
zdot_out = [];
y_ext = [];
ydot_ext = [];
for ii=1:size(yinit,1)
    for jj=1:size(Yi,1)
        [t_out_ode45,y_out_ode45] = ode45(@(t,y) fieh(t,y,Ca,C,M,Yi(jj,:),omega,kappa,stroke,mu,g,N,D,mu0,m,R,L),tspan,yinit(ii,:));
        % a = fieh(t_out,y_out',Ca,C,M,Yi,omega*ones(1,numel(y_out(:,3))),kappa,stroke,mu,g,N,D,mu0,m,R,L)';

        %
        z_out_ode45 = [(y_out_ode45(:,1)-Yi(jj,:)*sin(y_out_ode45(:,3))) , (y_out_ode45(:,2)-Yi(jj,:)*omega*cos(y_out_ode45(:,3)))]; %augmented data state
        zdot_out_ode45 = diff(z_out_ode45,1,1); %Compute the first-order difference
        t_out_ode45 = t_out_ode45(1:end-1,:);
        z_out_ode45 = z_out_ode45(1:end-1,:);
        y_out_ode45 = y_out_ode45(1:end-1,:);
        
        % store only stable dynamic data
        clear ind_stab
        t_out_ode45 = t_out_ode45(t_out_ode45>=9);
        ind_stab = size(z_out_ode45,1)-size(t_out_ode45,1)+1;
        z_out_ode45 = z_out_ode45(ind_stab:end,:);
        zdot_out_ode45 = zdot_out_ode45(ind_stab:end,:);
        y_ext_ode45 = Yi(jj,:)*sin(y_out_ode45(ind_stab:end,3));
        ydot_ext_ode45 = Yi(jj,:)*omega*cos(y_out_ode45(ind_stab:end,3));
        
        % t_out = t_out(3805:3971);
        % z_out = z_out(3805:3971,:);
        % zdot_out = zdot_out(3805:3971,:);
        % y_ext = Yi*sin(y_out(3805:3971,3));
        % ydot_ext = Yi*omega*cos(y_out(3805:3971,3));


        t_out = [t_out; t_out_ode45];
        y_out = [y_out; y_out_ode45];
        z_out = [z_out; z_out_ode45];
        zdot_out = [zdot_out; zdot_out_ode45];
        y_ext = [y_ext; y_ext_ode45];
        ydot_ext = [ydot_ext; ydot_ext_ode45];        
    end
end

figure(1)
plot(t_out,z_out(:,1),'b--')
hold on
plot(t_out,z_out(:,2),'r')
xlabel('$t$ [ ]','Interpreter','latex')
ylabel('$\dot{x}$ [ ]','Interpreter','latex')
title('FIEH System')

%plot - check phase space
figure
plot(z_out(:,1),z_out(:,2),'k')
xlabel('$x$ [ ]','Interpreter','latex')
ylabel('$\dot{x}$ [ ]','Interpreter','latex')
title('phase space data')

%% plot
plot_index_switch = find(z_out(:,1) == zswitch_median_down);
%setting the initial condition of the validation data, from which to run
%the model discovered
y0_impactDown = [-0.0159874  -1.2847]; %for data set ...
y0_impactTop = [0.0156102   1.32878]; %for data set ...
y0_free1 = [-0.0149886   0.675417]; %for data set ...
y0_free2 = [0.0147908   -0.700633]; %for data set ...
   
tspan_impact = 0;
tspan_free = 0;
% for mm=1:numel(t_out_impactDown)-1+24000 %we cut 144 points to avoid that springing dynamic continue beyond its own manifold
%     tspan_impact = [tspan_impact mm*dt];
% end
% for mm=1:numel(t_out_free)-1+110000%we cut 268 points to avoid that flying dynamic continue beyond its own manifold
%     tspan_free = [tspan_free mm*dt];
% end
Xicomb = Cluster_impactDown.Xistruct.Xicomb;
for kk = 10%:length(Xicomb)
    disp('Cluster down impact')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_impact_down = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_impactDown,x_mod_impactDown] = ode45(@(t,x)sparseGalerkin(t,x,Xi_impact_down,polyorder_impactDown,polyorder_impactDown_x,polyorder_impactDown_y,usesine_impactDown),[0 48],y0_impactDown,options);  % approximate
    a_mod_impactDown = diff(x_mod_impactDown(:,2));
%     a_mod_springDown = sparseGalerkin(t_mod_spring,x_mod_spring',Xi_spring,polyorder_spring,polyorder_spring_x,polyorder_spring_y,usesine_spring);  % approximate
end
%
clear Xicomb
Xicomb = Cluster_impactTop.Xistruct.Xicomb;
for kk = 42%:length(Xicomb)
    disp('Cluster top impact')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_impact_top = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_impactTop,x_mod_impactTop] = ode45(@(t,x)sparseGalerkin(t,x,Xi_impact_top,polyorder_impactTop,polyorder_impactTop_x,polyorder_impactTop_y,usesine_impactTop),[0 48],y0_impactTop,options);  % approximate
    a_mod_impactTop = diff(x_mod_impactTop(:,2));
%     a_mod_spring = sparseGalerkin(t_mod_spring,x_mod_spring',Xi_spring,polyorder_spring,polyorder_spring_x,polyorder_spring_y,usesine_spring);  % approximate
end
%
clear Xicomb
Xicomb = Cluster_free.Xistruct.Xicomb;
LL=6
for kk = LL%:length(Xicomb)
    disp('Cluster free')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_free = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_free1,x_mod_free1] = ode45(@(t,x)sparseGalerkin(t,x,Xi_free,polyorder_free,polyorder_free_x,polyorder_free_y,usesine_free),[0 32],y0_free1,options);  % approximate
    a_mod_free1 = diff(x_mod_free1(:,2));
%     a_mod_flight = sparseGalerkin(t_mod_flight,x_mod_flight',Xi_flight,polyorder_flight,polyorder_flight_x,polyorder_flight_y,usesine_flight);  % approximate
end
%
clear Xicomb
Xicomb = Cluster_free.Xistruct.Xicomb;
for kk = LL%:length(Xicomb)
    disp('Cluster free')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_free = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_free2,x_mod_free2] = ode45(@(t,x)sparseGalerkin(t,x,Xi_free,polyorder_free,polyorder_free_x,polyorder_free_y,usesine_free),[0 30],y0_free2,options);  % approximate
%     [t_mod_free2,x_mod_free2] = ode45(@(t,x)sparseGalerkin(t,x,Xi_free,polyorder_free,polyorder_free_x,polyorder_free_y,usesine_free),[],y0_free2,options);  % approximate
    a_mod_free2 = diff(x_mod_free2(:,2));
%     a_mod_flight = sparseGalerkin(t_mod_flight,x_mod_flight',Xi_flight,polyorder_flight,polyorder_flight_x,polyorder_flight_y,usesine_flight);  % approximate
end

figure
plot(z_out(:,1),z_out(:,2),'k--','LineWidth',0.9) % measurements
hold on
plot(x_mod_impactDown(:,1),x_mod_impactDown(:,2),'r','LineWidth',0.9) % approx model
plot(x_mod_free1(:,1),x_mod_free1(:,2),'b','LineWidth',0.9) % approx model
plot([zswitch_median_down zswitch_median_down], [1.3*min(z_out(:,2)) 1.3*max(z_out(:,2))],'k','LineWidth',0.9)
plot([zswitch_median_top zswitch_median_top], [1.3*min(z_out(:,2)) 1.3*max(z_out(:,2))],'k','LineWidth',0.9)
plot(x_mod_impactTop(:,1),x_mod_impactTop(:,2),'r','LineWidth',0.9) % approx model
plot(x_mod_free2(:,1),x_mod_free2(:,2),'b','LineWidth',0.9) % approx model
xlabel('$x [ ]$','interpreter','latex')
ylabel('$\dot{x} [ ]$','interpreter','latex')
title('HySINDy - Phase space','interpreter','latex')
legend('training data','impact','free','$x{switch}$','interpreter','latex')


figure
plot(t_mod_impactDown(2:end),a_mod_impactDown,'r','LineWidth',0.9) % approx model
hold on
t1plot = t_mod_free1+t_mod_impactDown(end);
plot(t1plot(2:end),a_mod_free1,'b','LineWidth',0.9) % approx model
t2plot = t1plot(end)+t_mod_impactTop;
plot(t2plot(2:end),a_mod_impactTop,'r','LineWidth',0.9) % approx model
t3plot = t2plot(end)+t_mod_free2;
plot(t3plot(2:end),a_mod_free2,'b','LineWidth',0.9) % approx model
% plot(t_out,zdot_out(:,2),'k--','LineWidth',0.9) % measurements
xlabel('$t [ ]$','interpreter','latex')
ylabel('$\stackrel{..}{x} [ ]$','interpreter','latex');
title('HySINDy - Acceleration','interpreter','latex')
legend('impact','free','interpreter','latex')

figure
subplot(2,1,1)
plot(t_out./dt,z_out(:,1),'k--','LineWidth',0.9) % measurements
hold on
plot(t_mod_impactDown,x_mod_impactDown(:,1),'r','LineWidth',0.9) % approx model
plot(t1plot,x_mod_free1(:,1),'b','LineWidth',0.9) % approx model
plot(t2plot,x_mod_impactTop(:,1),'r','LineWidth',0.9) % approx model
plot(t3plot,x_mod_free2(:,1),'b','LineWidth',0.9) % approx model
xlabel('$t [ ]$','interpreter','latex')
ylabel('$x [ ]$','interpreter','latex');
title('HySINDy - Position','interpreter','latex')
legend('training data','impact','free','interpreter','latex')

subplot(2,1,2)
plot(t_out./dt,z_out(:,2),'k--','LineWidth',0.9) % measurements
hold on
plot(t_mod_impactDown,x_mod_impactDown(:,2),'r','LineWidth',0.9) % approx model
plot(t1plot,x_mod_free1(:,2),'b','LineWidth',0.9) % approx model
plot(t2plot,x_mod_impactTop(:,2),'r','LineWidth',0.9) % approx model
plot(t3plot,x_mod_free2(:,2),'b','LineWidth',0.9) % approx model
xlabel('$t [ ]$','interpreter','latex')
ylabel('$\dot{x} [ ]$','interpreter','latex');
title('HySINDy - velocity','interpreter','latex')
legend('training data','impact','free','interpreter','latex')
