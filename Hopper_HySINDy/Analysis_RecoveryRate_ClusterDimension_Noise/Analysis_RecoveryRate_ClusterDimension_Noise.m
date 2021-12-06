%% Hopper Noise and Cluster dimension analysis
% Here is reported (with slight changes) the code developed by Mangan et.al
% to establish the influence of noise and cluster dimension on the recovery
% rate of SINDy algorthm.

% This script is autonomous.
% Run to carry out an analysis to establish the influence of noise and
% cluster dimension on the recovery rate.
%% create data sets for noise and cluster size plots

close all
clear all
clc
addpath('hopper'); addpath('SINDy'); addpath('validation');
changeplot
p.plottag =0;


% kappa = k*y0/m*g  % nondimensional param
p.kappa = 10;

% variables are: time, t, vertical position, y, vertical velocity, yp

tend = 5;
dt = 0.01;
p.tspan = [0:dt:tend];
p.options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_hopper);

% start with mass slightly below it's equilibrium position
% velocity downward
rng(2)
p.yinitvec = [abs(1+0.5*randn(100,1)), 0.5*randn(100,1)];
% noise conditions
% epsvec = [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1];
epsvec = logspace(-4, 1,50);

polyorder = 2;  % search space up to this order polynomials
usesine = 0; laurent = 0;   dyin =0;  dyorder = 0;
Thetalib.polyorder = polyorder;
Thetalib.usesine = usesine;
Thetalib.normTheta = 1;

% define correct Xi for checking

if polyorder == 2
    Xi_truefly = zeros(6,2);
    Xi_truehop = zeros(6,2);
elseif polyorder == 3
    Xi_truefly = zeros(10,2);
    Xi_truehop = zeros(10,2);
end
Xi_truehop(3,1) = 1;
Xi_truehop(1,2) = 11;
Xi_truehop(2,2) = -10;
Xi_truefly(3,1) = 1;
Xi_truefly(1,2) = -1;

%parameters for lambda sweep
lambdavals.numlambda = 20;
lambdavals.lambdastart = -2;
lambdavals.lambdaend = 2;
lambdavals.method = 'reg';

% run clean data
[data]  = RUN_Hopper(p);
% divide clean data
dtol = 1e-5;
[datafly, datahop] = dividehopper(data,dtol);

% run validation data
p.yinitvec = [abs(1+0.5*randn(10,1)), 0.5*randn(10,1)];
[valdata]  = RUN_Hopper(p);
[valfly, valhop] = dividehopper(valdata,dtol);

%%
for kk = 1:length(epsvec)
    rng(1)
    eps = epsvec(kk) % noise level
    for mm = 1:20 % do the whole thing this many times
        datafly2 = datafly;
        datahop2 = datahop;
        datafly2.y = datafly.y + eps*randn(size(datafly.y));
        datahop2.y = datahop.y + eps*randn(size(datahop.y));
        
        %     prep validation data
        valfly2 = valfly; valhop2 = valhop;
        valfly2.y = valfly.y + eps*randn(size(valfly.y));
        valhop2.y = valhop.y + eps*randn(size(valhop.y));
        % prep validation data
        Thetahopval = poolDatady_for_Ex_Jumping_Noise_Datasets(valhop2.y,2,polyorder,usesine, laurent, dyin, dyorder);
        Thetaflyval = poolDatady_for_Ex_Jumping_Noise_Datasets(valfly2.y,2,polyorder,usesine, laurent, dyin, dyorder);
        
        % find clusters
        %     clustersize = 10:20:10000;
        clustersize = floor(logspace(1, log10(min(length(datafly.y), length(datahop.y))),40));
        [hop_xy, fly_xy] = ClusterFlyHop(datafly2, datahop2, clustersize, p.plottag);
        % add cluster which is all data in hop and fly
        hop_xy{length(hop_xy)+1} =1:length(datahop2.t);
        fly_xy{length(fly_xy)+1} =1:length(datafly2.t);
        
        % make data matrices for each size cluster
        % create function library
        
        for ll =  1:length(clustersize)
            clustersize(ll)
            % concatenate data for hopping
            Xhop = [datahop2.y(hop_xy{ll}',:)];
            XPhop = [datahop2.a(hop_xy{ll}',:)];

            % conatenate data for flying
            Xfly = [datafly2.y(fly_xy{ll}',:)];
            XPfly = [datafly2.a(fly_xy{ll}',:)];
            
            Thetahop = poolDatady_for_Ex_Jumping_Noise_Datasets(Xhop,2,polyorder,usesine, laurent, dyin, dyorder);
            Thetafly = poolDatady_for_Ex_Jumping_Noise_Datasets(Xfly,2,polyorder,usesine, laurent, dyin, dyorder);
            
            %Calc condition number
            condTheta(kk,ll,mm,1) = cond(Thetahop);
            condTheta(kk,ll,mm,2) = cond(Thetafly);

            
            % run SINDy on cluster
            Thetalib.Theta = Thetahop;
            Thetalib.dx = XPhop;
            Xistructhop =multiD_Lambda_for_Ex_Jumping_Noise_Datasets(Thetalib,lambdavals);
            
            Thetalib.Theta = Thetafly;
            Thetalib.dx = XPfly;
            Xistructfly =multiD_Lambda_for_Ex_Jumping_Noise_Datasets(Thetalib,lambdavals);
            
            % check if correct model was recovered
            clear Xihop Xifly numXihop numXifly
            Xihop = reshape(cell2mat(Xistructhop.Xicomb), [6],[2],[])
            Xifly = reshape(cell2mat(Xistructfly.Xicomb), [6],[2],[])
            numXihop = size(Xihop,3);
            numXifly = size(Xifly,3);
            
            
            clear errorhop1_a errorhop1_v errorfly1_a errorfly1_v
            
            foundhop = 0;
            for rr = 1:numXihop
                foundhop = foundhop + ~any(any((abs(Xi_truehop)>0) - (abs(Xihop(:,:,rr))>0)))
                if ~any(any((abs(Xi_truehop)>0) - (abs(Xihop(:,:,rr))>0))) % if correct
                    numnonzero = nnz(Xihop(:,:,rr));
                end
                % calc error with validation data
                errorhop1_v(rr) = sqrt(sum((valhop2.a(:,1)-Thetahopval*Xihop(:,1,rr)).^2));
                errorhop1_a(rr) = sqrt(sum((valhop2.a(:,2)-Thetahopval*Xihop(:,2,rr)).^2));
            end
            foundfly = 0;
            for rr =1:numXifly
                foundfly = foundfly + ~any(any((abs(Xi_truefly)>0) - (abs(Xifly(:,:,rr))>0)))
                % calc error with validation data
                errorfly1_v(rr) = sqrt(sum((valfly2.a(:,1)-Thetaflyval*Xifly(:,1,rr)).^2));
                errorfly1_a(rr) = sqrt(sum((valfly2.a(:,2)-Thetaflyval*Xifly(:,2,rr)).^2));
            end
            hopyes(kk,ll,mm) = foundhop;
            flyyes(kk,ll,mm) = foundfly;
            errorhop_v(kk,ll,mm) = min(errorhop1_v);
            errorhop_a(kk,ll,mm) = min(errorhop1_a);
            errorfly_v(kk,ll,mm) = min(errorfly1_v);
            errorfly_a(kk,ll,mm) = min(errorfly1_a);
           
            % find the number of nonzero entries for model with min error
            clear minerror_modelnum_hop minerror_modelnum_fly coeff_nonzero_hop coeff_nonzero_fly
            clear infXi_hop infXi_fly
            minerror_modelnum_hop =intersect(find(min(errorhop1_a)==errorhop1_a),find(min(errorhop1_v)==errorhop1_v));
            coeff_nonzero_hop = nnz(Xihop(:,:,minerror_modelnum_hop))
            infXi_hop = norm(Xihop(:,:,minerror_modelnum_hop),'inf')
            minerror_modelnum_fly =intersect(find(min(errorfly1_a)==errorfly1_a),find(min(errorfly1_v)==errorfly1_v));
            coeff_nonzero_fly = nnz(Xifly(:,:,minerror_modelnum_fly))
            infXi_fly = norm(Xifly(:,:,minerror_modelnum_fly),'inf')
            
            % find the smallest coefficient of the same model
            clear smallest_coeff_hop smallest_coeff_fly
            smallest_coeff_hop = min(min(abs(nonzeros(Xihop(:,:,minerror_modelnum_hop)))))
            smallest_coeff_fly = min(min(abs(nonzeros(Xifly(:,:,minerror_modelnum_fly)))))
            
            % calculate alternative bound from appendix
            boundhop(kk,ll,mm) = sqrt(coeff_nonzero_hop)*infXi_hop*condTheta(kk,ll,mm,1)*eps/(1-condTheta(kk,ll,mm,1)*eps);
            boundfly(kk,ll,mm) = sqrt(coeff_nonzero_fly)*infXi_fly*condTheta(kk,ll,mm,2)*eps/(1-condTheta(kk,ll,mm,2)*eps);
            boundhop2(kk,ll,mm) = sqrt(coeff_nonzero_hop)*infXi_hop*condTheta(kk,ll,mm,1)*eps;
            boundfly2(kk,ll,mm) = sqrt(coeff_nonzero_fly)*infXi_fly*condTheta(kk,ll,mm,2)*eps;
        end        
        
    end
end
%%
dateformatout = 'mmddyyyy';
save([datestr(now, dateformatout) 'Jumping_Noise_Datasets.mat'])

%% take averages across the instances  (3rd variable)
pathname = '';
mkdir(pathname);
changeplot
%%
hopyesavg = mean(hopyes, 3);
flyyesavg = mean(flyyes, 3);
condThetaavg = mean(condTheta,3);


%% Plots


dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'SpringDynamics.fig'];
f = figure(14);
[C,h]=contourf(epsvec, clustersize,  hopyesavg');
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Model found for Spring Dynamics')
f.Children.XScale = 'log'; f.Children.YScale = 'log';
colormapnew
c = colorbar; 
% colorbar('off')
h.LevelStep = 0.2;
axis([1e-4 10 10 1.5e4])
f.CurrentAxes.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
f.CurrentAxes.YTick = [1e1 1e2 1e3 1e4];

%%
savefig(plotname);
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('');
print([plotname(1:end-4) '.svg'],'-dsvg')
%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'FlyDynamics.fig'];
f = figure(20);
[C,h]=contourf(epsvec, clustersize,  flyyesavg');
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Model found for Flying Dynamics')
f.Children.XScale = 'log'; f.Children.YScale = 'log';
colormapnew
c = colorbar; 
% colorbar('off')
h.LevelStep = 0.2;
axis([1e-4 10 10 1.5e4])
f.CurrentAxes.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
f.CurrentAxes.YTick = [1e1 1e2 1e3 1e4];
savefig(plotname);
%%
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('');
print([plotname(1:end-4) '.svg'],'-dsvg')
%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'CondSpringDynamics.fig'];
f = figure(16);
contourf(epsvec, clustersize, log10(condThetaavg(:,:,:,1)'))
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Spring Dynamics')
colormapnew
c = colorbar;
c.Label.String = 'k = log_{10}(condition #)';
hold on
caxis([0 6])
[C,h]=contour(epsvec, clustersize,  hopyesavg');
h.LineColor = [0.5 0.5 0.5];
h.LineStyle = '-';h.LineWidth = 1; h.LevelStep = 0.5;
h.Parent.XScale = 'log'; h.Parent.YScale = 'log';
axis([1e-4 10 10 1.5e4])
h.Parent.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
h.Parent.YTick = [1e1 1e2 1e3 1e4];
savefig(plotname);
%%
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('');
print([plotname(1:end-4) '.svg'],'-dsvg')
%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'CondFlyDynamics.fig'];
f = figure(17);
contourf(epsvec, clustersize, log10(condThetaavg(:,:,:,2)'))
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Flying Dynamics')
colormapnew
c = colorbar;
c.Label.String = 'k = log_{10}(condition #)';
hold on
caxis([0 6])
[C,h]=contour(epsvec, clustersize,  flyyesavg');
h.Parent.XScale = 'log'; h.Parent.YScale = 'log';
h.LineColor = [0.5 0.5 0.5];
h.LineStyle = '-';h.LineWidth = 1; h.LevelStep = 0.5;
axis([1e-4 10 10 1.5e4])
h.Parent.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
h.Parent.YTick = [1e1 1e2 1e3 1e4];
%%
savefig(plotname);
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('')
print([plotname(1:end-4) '.svg'],'-dsvg')
%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'CondEpsSpringDynamics.fig'];
f = figure(18);
contourf(epsvec, clustersize, log10(repmat(epsvec', [1 size(condThetaavg,2)])'.*condThetaavg(:,:,:,1)'))
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Spring Dynamics')
colormapnew
c = colorbar;
c.Label.String = 'k = log_{10}(\epsilon \times condition #)';
hold on
caxis([-2 6])
[C,h]=contour(epsvec, clustersize,  hopyesavg');
h.LineColor = [0.5 0.5 0.5];
h.LineStyle = '-';h.LineWidth = 1; h.LevelStep = 0.5;
h.Parent.XScale = 'log'; h.Parent.YScale = 'log';
axis([1e-4 10 10 1.5e4])
h.Parent.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
h.Parent.YTick = [1e1 1e2 1e3 1e4];
savefig(plotname);
%%
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('');
print([plotname(1:end-4) '.svg'],'-dsvg')

%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'CondEpsFlyDynamics.fig'];
f = figure(19);
eps_cond =repmat(epsvec', [1 size(condThetaavg,2)])'.*condThetaavg(:,:,:,2)';
% contourf(epsvec, clustersize, log10(repmat(epsvec', [1 size(condThetaavg,2)])'.*condThetaavg(:,:,:,2)'))
contourf(epsvec, clustersize, eps_cond)
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Flying Dynamics')
colormapnew
c = colorbar;
c.Label.String = 'k = log_{10}(\epsilon \times condition #)';
hold on
caxis([-2 6])
[C,h]=contour(epsvec, clustersize,  flyyesavg');
h.Parent.XScale = 'log'; h.Parent.YScale = 'log';
h.LineColor = [0.5 0.5 0.5];
h.LineStyle = '-';h.LineWidth = 1; h.LevelStep = 0.5;
axis([1e-4 10 10 1.5e4])
h.Parent.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
h.Parent.YTick = [1e1 1e2 1e3 1e4];
savefig(plotname);
%%
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel('')
ylabel('')
title('')
print([plotname(1:end-4) '.svg'],'-dsvg')


