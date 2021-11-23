%% FIEH_PLOT - plots for paper figure FIEH
% input a save location in the line below:
pathname = 'plot';
% mkdir(pathname);

dateformatout = 'mmddyyyy';

addpath('FIEH'); addpath('SINDy'); addpath('validation');
%%
%% plot
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

plot_index_switch = find(z_out(:,1) == zswitch_median_down);
y0_impactDown = [-0.015  -0.125519]; %for data set ...
y0_impactTop = [0.015   0.125653]; %for data set ...
y0_free1 = [-0.015   0.0621174]; %for data set ...
y0_free2 = [0.015   -0.061806]; %for data set ...
   
tspan_impact = 0;
tspan_free = 0;
for mm=1:numel(t_out_impactDown)-1+2400 %we cut 144 points to avoid that springing dynamic continue beyond its own manifold
    tspan_impact = [tspan_impact mm*dt];
end
for mm=1:numel(t_out_free)-1 +6700%we cut 268 points to avoid that flying dynamic continue beyond its own manifold
    tspan_free = [tspan_free mm*dt];
end
Xicomb = Cluster_impactDown.Xistruct.Xicomb;
for kk = 1%:length(Xicomb)
    disp('Cluster down impact')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_impact_down = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_impactDown,x_mod_impactDown] = ode45(@(t,x)sparseGalerkin(t,x,Xi_impact_down,polyorder_impact,polyorder_impact_x,polyorder_impact_y,usesine_impact),tspan_impact,y0_impactDown,options);  % approximate
    a_mod_impactDown = diff(x_mod_impactDown(:,2));
%     a_mod_springDown = sparseGalerkin(t_mod_spring,x_mod_spring',Xi_spring,polyorder_spring,polyorder_spring_x,polyorder_spring_y,usesine_spring);  % approximate
end
%
clear Xicomb
Xicomb = Cluster_impactTop.Xistruct.Xicomb;
for kk = 1%:length(Xicomb)
    disp('Cluster top impact')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_impact_top = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_impactTop,x_mod_impactTop] = ode45(@(t,x)sparseGalerkin(t,x,Xi_impact_top,polyorder_impact,polyorder_impact_x,polyorder_impact_y,usesine_impact),tspan_impact,y0_impactTop,options);  % approximate
    a_mod_impactTop = diff(x_mod_impactTop(:,2));
%     a_mod_spring = sparseGalerkin(t_mod_spring,x_mod_spring',Xi_spring,polyorder_spring,polyorder_spring_x,polyorder_spring_y,usesine_spring);  % approximate
end
%
clear Xicomb
Xicomb = Cluster_free.Xistruct.Xicomb;
for kk = 35%:length(Xicomb)
    disp('Cluster free')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_free = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_free1,x_mod_free1] = ode45(@(t,x)sparseGalerkin(t,x,Xi_free,polyorder_free,polyorder_free_x,polyorder_free_y,usesine_free),tspan_free,y0_free1,options);  % approximate
    a_mod_free1 = diff(x_mod_free1(:,2));
%     a_mod_flight = sparseGalerkin(t_mod_flight,x_mod_flight',Xi_flight,polyorder_flight,polyorder_flight_x,polyorder_flight_y,usesine_flight);  % approximate
end
%
clear Xicomb
Xicomb = Cluster_free.Xistruct.Xicomb;
for kk = 35%:length(Xicomb)
    disp('Cluster free')
    disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
    Xi_free = [Xicomb{kk}(:,1), Xicomb{kk}(:,2)]
    [t_mod_free2,x_mod_free2] = ode45(@(t,x)sparseGalerkin(t,x,Xi_free,polyorder_free,polyorder_free_x,polyorder_free_y,usesine_free),tspan_free,y0_free2,options);  % approximate
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
plot(t_out./dt-1171,z_out(:,1),'k--','LineWidth',0.9) % measurements
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
plot(t_out./dt-1171,z_out(:,2),'k--','LineWidth',0.9) % measurements
hold on
plot(t_mod_impactDown,x_mod_impactDown(:,2),'r','LineWidth',0.9) % approx model
plot(t1plot,x_mod_free1(:,2),'b','LineWidth',0.9) % approx model
plot(t2plot,x_mod_impactTop(:,2),'r','LineWidth',0.9) % approx model
plot(t3plot,x_mod_free2(:,2),'b','LineWidth',0.9) % approx model
xlabel('$t [ ]$','interpreter','latex')
ylabel('$\dot{x} [ ]$','interpreter','latex');
title('HySINDy - velocity','interpreter','latex')
legend('training data','impact','free','interpreter','latex')

%% Manifolds Library
% plot coefficients

% delta_aic_c = Cluster_impact1.aic_c-min(Cluster_impact1.aic_c)
nLib = length(Cluster_free.Xistruct.Xicomb); %number of libraries in the chosen manifold 'free'
clustnum = 1; %Manifold impact1
figure(clustnum)
ii =1;
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'Coeff' num2str(clustnum) '.fig']
for  kk = 1:nLib
    
    vc_ITop = Cluster_impactTop.Xistruct.Xicomb{(6*kk-5)}(:,1) %velocity coefficient
    ac_ITop = Cluster_impactTop.Xistruct.Xicomb{(6*kk-5)}(:,2) %acceleration coefficient
    %     omegac = Cluster_impactTop.Xistruct.Xicomb{kk}(:,3) %omega coefficient
    
    vc_IDown = Cluster_impactDown.Xistruct.Xicomb{(2*kk-1)}(:,1) %velocity coefficient
    ac_IDown = Cluster_impactDown.Xistruct.Xicomb{(2*kk-1)}(:,2) %acceleration coefficient
    %     omegac = Cluster_impactDown.Xistruct.Xicomb{kk}(:,3) %omega coefficient
    
    vc_free = Cluster_free.Xistruct.Xicomb{kk}(:,1) %velocity coefficient
    ac_free = Cluster_free.Xistruct.Xicomb{kk}(:,2) %acceleration coefficient
    %     omegac = Cluster_impactTop.Xistruct.Xicomb{kk}(:,3) %omega coefficient
    
    subplot(nLib,2,ii)
    %     subplot(2,nLib,ii)
    %Fourier base
    %     baseLibrary = categorical({'1','sin(x)','sin(2x)','sin(3x)','sin(xdot)','sin(2xdot)'...
    %         ,'cos(x)','cos(2x)','cos(3x)','cos(xdot)','cos(2xdot)','cos(\tau)','sin(\tau)'});
    %     baseLibrary = reordercats(baseLibrary,{'1','sin(x)','sin(2x)','sin(3x)','sin(xdot)','sin(2xdot)'...
    %         ,'cos(x)','cos(2x)','cos(3x)','cos(xdot)','cos(2xdot)','cos(\tau)','sin(\tau)'});
    %Polynomial base
    baseLibrary = categorical({'1','x','xdot','x^2','x*xdot','xdot^2'});
    baseLibrary = reordercats(baseLibrary,{'1','x','xdot','x^2','x*xdot','xdot^2'});
    b1 = bar(baseLibrary,[abs(vc_free)'; abs(vc_IDown)'; abs(vc_ITop)']);
%     b1.BarWidth = 0.2;
%     box off
%     b1.Parent.YScale = 'log';
%     b1.Parent.YColor = [1 1 1];
    %         b1.Parent.XAxis.TickLength = [0 0];
    %         b1.Parent.XAxis.Limits = [0.5 6.5];
    %         b1.Parent.XAxis.LineWidth=1;
%     b1.Parent.FontSize = 9;
    %         b1.Parent.Parent.PaperPosition = [0 0 3.5 2.5];
    %         b1.Parent.Parent.PaperSize = [3.5 2.5];
    xlabel('num. coefficient')
    title(['library velocity coefficients \Delta{aic_c}' ])
%     if  kk ~= nLib
%         b1.Parent.XAxis.TickValues = [];
%     end
    ii = ii+1;
    %
    subplot(nLib,2,ii)
    %     subplot(2,nLib,ii)
    b2 = bar(baseLibrary,[abs(ac_free)'; abs(ac_IDown)'; abs(ac_ITop)']);
%     b2.BarWidth = 0.2;
%     box off
%     b2.Parent.YScale = 'log';
%     b2.Parent.YColor = [1 1 1];
    %         b2.Parent.XAxis.TickLength = [0 0];
    %         b2.Parent.XAxis.Limits = [0.5 6.5];
    %         b2.Parent.XAxis.LineWidth=1;
%     b2.Parent.FontSize = 9;
    %         b2.Parent.Parent.PaperPosition = [0 0 3.5 2.5];
    %         b2.Parent.Parent.PaperSize = [3.5 2.5];
    xlabel('num. coefficient')
    title(['library acceleration \Delta{aic_c}'])
%     if  kk ~= nLib
%         b2.Parent.XAxis.TickValues = [];
%     end
    ii= ii+1;
    %
    %         subplot(nLib,3,ii)
% %     subplot(3,nLib,ii)
% %     b3 = bar(baseLibrary,abs(omegac));
% %     b3.BarWidth = 0.2;
% %     box off
% %     b3.Parent.YScale = 'log';
% %     b3.Parent.YColor = [1 1 1];
% %     %         b3.Parent.XAxis.TickLength = [0 0];
% %     %         b3.Parent.XAxis.Limits = [0.5 6.5];
% %     %         b3.Parent.XAxis.LineWidth=1;
% %     b3.Parent.FontSize = 9;
% %     %         b3.Parent.Parent.PaperPosition = [0 0 3.5 2.5];
% %     %         b3.Parent.Parent.PaperSize = [3.5 2.5];
% %     xlabel('num. coefficient')
% %     title(['library omega coefficients Cluster' num2str(clustnum)])
% %     if  kk ~= nLib
% %         b2.Parent.XAxis.TickValues = [];
% %     end
% %     ii= ii+1;
    
end
%     print([plotname(1:end-4) '.svg'],'-dsvg')


%% Figure of validation phase space
%% Figure of validation time series
%% relative AIC_c plot
close all
modelnum = 1; %  this is the modelnumber in the HighFreq_II library
 % use this code for models 1 and 2, need to run for 3 and 4 each
 % before running next cell
plotname = [pathname datestr(now, dateformatout) 'AIC_rel1_' num2str(modelnum) '.fig'];
% model1ind = find(min_aic_modelnum == modelnum);

AIC_rel = Cluster_impact1.aic_c - min(Cluster_impact1.aic_c);
numcoeff = Cluster_impact1.numterms(1:length(AIC_rel));
figure(ii)
semilogy(numcoeff, AIC_rel, 'ok')
ylabel('relative AIC_c')
xlabel('numcoeff')
title('relative AIC_c')
hold on

%%
% plot coefficients
% C2nz = C2; C2nz(C2==0) = NaN;
% figure
% plot(t_out(1:end-1),C2nz','.')
% xlabel('t [ ]','Interpreter','latex')
% ylabel('Coeff$\stackrel{..}{x}$ [ ]','Interpreter','latex')
% set(gca,'FontSize',9)
% title('Coeff$\stackrel{..}{x}$ vs time','Interpreter','latex')
% xlimits = [min(t_out) max(t_out)];
% ylimits = [min(min(C2))*1.05  max(max(C2))*1.05];
% axis([xlimits ylimits])
% % yticks = round(min(min(C2)),-1):100:round(max(max(C2)),-1)+50
% % xticks = [0 1 2 3 4 5 6];
% % PP = [0 0 3.5 1.5]*0.8;
% % PS = [3.5 1.5]*0.8;
% plotname = [pathname datestr(now, dateformatout) 'C2vtime.fig'];
% % savefigfile
% savefig(plotname);
    
%      ylimits = [-11 11]
%      yticks = -12:2:12
% plotname = [pathname datestr(now, dateformatout) 'C2_' num2str(ii) 'vtimezoom.fig'];
%      savefigfile  
%%
     
% for ii=1:size(C2,1)
%     plotname = [pathname datestr(now, dateformatout) 'C2_' num2str(ii) 'vtime.fig'];
% 
%     % plot(t_out(1:end-1),C2','*')
% %     if any(ii == [1 2])
% %         plot(t_out(1:end-1),abs(C2(ii,:))','*')
% %     else
%         plot(t_out(1:end-1),C2(ii,:)','.')
%         xlabel('t [ ]','Interpreter','latex')
%         ylabel('Coeff$\stackrel{..}{x}$ [ ]','Interpreter','latex')
%         set(gca,'FontSize',9)
%         title(['Coeff.' num2str(ii) ' $\stackrel{..}{x}$ vs time'],'Interpreter','latex')
%         xlimits = [min(t_out) max(t_out)];
%         ylimits = [min(min(C2(ii,:)))*1.05  max(max(C2(ii,:)))*1.05];
%         axis([xlimits ylimits])
%         savefig(plotname)
% %         yticks = round(min(C2(ii,:)),-2):100:round(max(C2(ii,:)),-2)
% %     end
% %     xlimits = [0 6];
% %     ylimits = [min([-10 C2(ii,:)])  max([C2(ii,:) 10] )];
% %     xticks = [0 1 2 3 4 5 6];
% %     
% % %     yticks = [5e0 1e1 1e2];
% %     PP = [0 0 3.5 1]*0.8;
% %     PS = [3.5 1]*0.8;
% %     savefigfile
% end
%%
% % plot coefficients
% C1nz = C1; C1nz(C1==0) = NaN;
% figure
% plot(t_out(1:end-1),C1nz','.')
% xlabel('t [ ]','Interpreter','latex')
% ylabel('Coeff$\dot{x}$ [ ]','Interpreter','latex')
% set(gca,'FontSize',9)
% title('Coeff$\dot{x}$ vs time','Interpreter','latex')
% xlimits = [min(t_out) max(t_out)];
% ylimits = [-8*10^6  8.5*10^6];
% axis([xlimits ylimits])
% % yticks = round(min(min(C1)),-1):100:round(max(max(C1)),-1)+50
% % xticks = [0 1 2 3 4 5 6];
% % PP = [0 0 3.5 1.5]*0.8;
% % PS = [3.5 1.5]*0.8;
% plotname = [pathname datestr(now, dateformatout) 'C1vtime.fig'];
% % savefigfile
% savefig(plotname);
%%
     
% for ii=1:size(C1,1)
%     plotname = [pathname datestr(now, dateformatout) 'C1_' num2str(ii) 'vtime.fig'];
% 
%     % plot(t_out(1:end-1),C1','*')
% %     if any(ii == [1 2])
% %         plot(t_out(1:end-1),abs(C1(ii,:))','*')
% %     else
%         plot(t_out(1:end-1),C1(ii,:)','.')
%         xlabel('t [ ]','Interpreter','latex')
%         ylabel('Coeff$\dot{x}$ [ ]','Interpreter','latex')
%         set(gca,'FontSize',9)
%         title(['Coeff.' num2str(ii) ' $\dot{x}$ vs time'],'Interpreter','latex')
% %         xlimits = [min(t_out) max(t_out)];
% %         ylimits = [min(min(C1(ii,:)))*1.05  max(max(C1(ii,:)))*1.05];
% %         axis([xlimits ylimits])
%         savefig(plotname)
% 
% %         yticks = round(min(C1(ii,:)),-2):100:round(max(C1(ii,:)),-2)
% %     end
% %     xlimits = [0 6];
% %     ylimits = [min([-10 C1(ii,:)])  max([C1(ii,:) 10] )];
% %     xticks = [0 1 2 3 4 5 6];
% %     
% % %     yticks = [5e0 1e1 1e2];
% %     PP = [0 0 3.5 1]*0.8;
% %     PS = [3.5 1]*0.8;
% %     savefigfile
% end
%% plot aicmin vs time
t1ind = floor(length(HMt_min)/3);
plotname = [pathname datestr(now, dateformatout) 'minaicval.fig'];
figure
semilogy(HMt(1,1:t1ind), HMerror_minaic(1,1:t1ind),'.')
hold on
plot(HMt(2,1:t1ind), HMerror_minaic(2,1:t1ind),'.')
plot(HMt(3,1:t1ind), HMerror_minaic(3,1:t1ind),'.')
plot(HMt(4,1:t1ind), HMerror_minaic(4,1:t1ind),'.')
plot(HMt(5,1:t1ind), HMerror_minaic(5,1:t1ind),'.')
xlabel('t [ ]','Interpreter','latex')
ylabel('AICmin [ ]','Interpreter','latex')
legend('HighFreqModel1','HighFreqModel2','HighFreqModel3','HighFreqModel4','HighFreqModel5')
set(gca,'FontSize',9)
title('AICmin vs time','Interpreter','latex')
% axis([0 6  0.7*min(min(HMerror_minaic(HMerror_minaic(:,1:t1ind)>0))) 1.5*max(max(HMerror_minaic(HMerror_minaic(:,1:t1ind)>0)))])
xlimits = [0 2.131];
ylimits = [0.7*min(min(HMerror_minaic(HMerror_minaic(:,1:t1ind)>0))) 1.5*max(max(HMerror_minaic(HMerror_minaic(:,1:t1ind)>0)))];
axis([xlimits ylimits]);


%% Pareto front
%Da fare
