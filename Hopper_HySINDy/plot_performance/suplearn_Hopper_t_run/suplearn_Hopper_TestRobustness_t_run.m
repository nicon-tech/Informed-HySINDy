%% t_run Hopping performance - supervised learning
% This script is autonomous.
%Run to carry out an analysis on the computational performance of the algorithm

close all
clear all
clc

p.plottag =1;
addpath('hopper'); addpath('SINDy'); addpath('validation');

%% controlled parameters
% noise
% cluster dimension k
% val_ntimes
% thick_unddata
% model stiffness

thick_unddata_vect = 0.01:0.01:0.5;
eps_vect = [1e-6 1e-4 1e-2];
numleveldataprocessed = 3; %number of level data processed
% tStart = zeros(1,size(thick_unddata_vect,2));
tEnd = zeros(numel(thick_unddata_vect),numleveldataprocessed,numel(eps_vect));
% a0_mean = zeros(numel(thick_unddata_vect),numleveldataprocessed,numel(eps_vect));
a0_median = zeros(numel(thick_unddata_vect),numleveldataprocessed,numel(eps_vect));

%%
for leveldataprocessed=1:3 %level of data processed
    
    clear p.yinitvec
    
    if leveldataprocessed == 1
        % start with mass slightly below it's equilibrium position
        % velocity downward
        p.yinitvec = [0.8 -0.1
            0.78 -0.1
            0.82 -0.1];
        
    elseif leveldataprocessed == 2
        % start with mass slightly below it's equilibrium position
        % velocity downward
        p.yinitvec = [0.8 -0.1
            0.78 -0.1
            0.82 -0.1
            0.81 -0.12];
        
    elseif leveldataprocessed == 3
        % start with mass slightly below it's equilibrium position
        % velocity downward
        p.yinitvec = [0.8 -0.1
            0.78 -0.1
            0.82 -0.1
            0.81 -0.12
            0.79 -0.13];
    else
    end
    
    ind_eps = 1; %initialization of the index associated with the noise
    for eps = eps_vect
        ind_thick_unddata = 1; %initialization of the index associated with thick unddata
        for thick_unddata = thick_unddata_vect % range undefined data position
            %% data generation
            % kappa = k*y0/m*g  % nondimensional parameter which represents the balance between the spring and gravity forces
            p.kappa = 10;
            
            p.eps = eps; % noise level
            
            rng(1)
            % In some situations, setting the seed alone will not guarantee the same
            % results. This is because the generator that the random number functions
            % draw from might be different than you expect when your code executes.
            % For long-term repeatability, specify the seed and the generator type together.
            % For example, the following code sets the seed to 1 and the generator to Mersenne Twister.
            % rng(1,'twister');
            
            % variables are: time, t, vertical position, y, vertical velocity, yp
            
            tend = 5;
            dt = 0.033;
            p.tspan = 0:dt:tend;
            
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
            
            %% undefined data
            tStart = tic;           % pair 1: tic
            
            plottag =0;
            dimMan = 2; % dimension of the Manifold
            
            % Run CCM on y_out
            clear x y
%             x= y_out_und(1:end-1,1); % height
            y= y_out_und(1:end-1,2); % velocity
            
            numOfPts = length(y); % number of points in the data
            %k = dimMan+2;
            k = 6; % number of points in a cluster
            
            idx_xy = knnsearch(y/max(y),y/max(y),'K',k);
            % similar dynamic data => each cluster contributes with its own library
            idx_xy = idx_xy(:,2:end);
            
            % run HySINDY on the undefined data with small (few data) stiff model
            % initialize a structure library
            II_Xi = [];
            for ii = 1:floor(numOfPts)
                
                % stack all the "position" values in cluster ii into a vector
                X = y(idx_xy(ii,:)');
                % stack all the "velocity" values in cluster ii into a vector
                XP = a_out_und(idx_xy(ii,:)',2);
                
                
                % create function library
                polyorder_und = 0;  % search space up to i-th order polynomials
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
                lambdavals.lambdastart = -5;
                lambdavals.lambdaend = 2;
                
                Xistruct =multiD_Lambda(Thetalib,lambdavals);
                Xicomb = Xistruct.Xicomb;
                
                % add new models to library and find the indices for this set
                [II_Xitemp, II_ind] = build_Xi_Library(Xicomb, II_Xi);
                II_Xi = II_Xitemp;
                
                Cluster{ii}.centroidy =  mean(X(:,1));
%                 Cluster{ii}.centroidy = mean(X(:,2));
                % store results by cluster number
                Cluster{ii}.Xistruct = Xistruct;
                Cluster{ii}.II_ind = II_ind;
                Cluster{ii}.X = X;
                Cluster{ii}.XP = XP;
            end
            
            %find switch point
            clear x2 y2 x y II_Xi
%             x2= y_out_und(1:end-1,1); % height
            y2= y_out_und(1:end-1,2); % velocity
            
            numOfPts = length(y2); % number of points in the data
            %     k = 16; % number of points in a cluster
            
            val_ntimes = floor(0.2/dt); %length of comparision not including ICs
            data2 = y2;
            
            Cluster_xswitch_max = [];
            
            options_test = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_diverge);
            
            for ii = 1:length(Cluster)
                Xicomb = Cluster{ii}.Xistruct.Xicomb;
                
                idx_xy2 = knnsearch(y2,...
                    Cluster{ii}.centroidy,'K',k);
                idx_xy2 = idx_xy2(:,2:end);
                
                [xA, x0clust, tvectest] = buildinit_fromcluster(val_ntimes, 1, idx_xy2, data2, dt, k);
                
                val.tA =tvectest;
                val.xA =xA;
                val.x0 =x0clust;
                val.options = options_test;
                
                clear IC aic_c aic_cv xswitch_max xswitch_min nodynam_any xswitch_avg xswitch
                clear nodynam abserror_bf_switch mean_e numterms idx_xy2
                for kk = 1:length(Xicomb)
                    %                 disp(['Cluster num:' num2str(ii) ' of' num2str(length(Cluster))])
                    %                 disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
                    Xi = Xicomb{kk};
                    
                    clear  savetB savexB
                    [savetB, savexB] = validation_sims(Xi, Thetalib, val, plottag);
                    
                    abserror_bf_switch = zeros(length(xA),2);
                    for ll = 1:length(xA)
                        clear xswitch nodynam1 abserror abserror_avg1 RMSE xAcomp xBcomp
                        [tlength, nterms] = size(savexB{ll});
                        [xswitch, ~, abserror_avg1, RMSE] = calc_tlength(xA{ll}, savexB{ll},val);
                        xswitch(ll,:) = xswitch;
                        abserror_bf_switch(ll,:) = abserror_avg1; %abs error befor that the switching takes place
                    end
                    abserror = abserror_bf_switch(:,1);
                    IC{kk} = ICcalculations(abserror, nnz(Xi),size(x0clust,2));
                    aic_c(kk,:) = IC{kk}.aic_c;
                    xswitch_max(kk,:) = max(xswitch);
                    xswitch_min(kk,:) = min(xswitch);
                    xswitch_avg(kk,:) = mean(xswitch);
                    mean_e(kk,:) = mean(abserror_bf_switch);
                    
                    
                end
                
                minIC = min(aic_c);
                
                Cluster_xswitch_max = [Cluster_xswitch_max ; mean(xswitch_max)]; %calculated with mean
                %             Cluster_xswitch_max_cmedian = [Cluster_xswitch_max_cmedian ; median(xswitch_max)]; %calculated with median
                %             Cluster_xswitch_max_ctrimmean = [Cluster_xswitch_max ; trimmean(xswitch_max,10,'round')]; %calculated with trimmean
                %             %trim
                %             percentlow=10;
                %             trimlow = floor((numel(xswitch_max))*(percent/100)/2);
                %             xswitch_max=sort(xswitch_max);
                %             xswitch_max=xswitch_max(trimlow:numel(xswitch_max)-trimhigh);
                %             Cluster_xswitch_max_ctrimmean = [Cluster_xswitch_max_ctrimmean ; mean(xswitch_max)]; %calculated with trimmean
                
                %     Cluster{ii}.xswitch_max = min(xswitch_max);
                %     Cluster{ii}.abserror_bf = mean_e;
                %     Cluster{ii}.IC = IC;
                %     Cluster{ii}.aic_c =aic_c;
                
            end
            clear Cluster Theta normTheta ThetaN Thetalib lambdavals Xistruct Xicomb
            clear II_ind val Xi
            
            %% Use kmeans to divide the correct cluster of possible switch point from the other
            % plot(Cluster_xswitch_max,'*') %pre-analysis of switch point detection
            % The xswitch data set(Cluster_xswitch_max) is divided into multiple
            % clusters to distinguish the clusters of the true xswitch from the false
            % ones detected
            optskmeans = statset('Display','final');
            [idxkmeans,Ckmeans] = kmeans(Cluster_xswitch_max,7,'Distance','cityblock',...
                'Replicates',7,'Options',optskmeans);
            %         figure(2); %clustering of possible switches
            %         plot(Cluster_xswitch_max(idxkmeans==1,1),'r.','MarkerSize',11)
            %         hold on
            %         plot(Cluster_xswitch_max(idxkmeans==2,1),'b.','MarkerSize',11)
            %         plot(Cluster_xswitch_max(idxkmeans==3,1),'g.','MarkerSize',11)
            %         plot(Cluster_xswitch_max(idxkmeans==4,1),'k.','MarkerSize',11)
            %         plot(Cluster_xswitch_max(idxkmeans==5,1),'y.','MarkerSize',11)
            %         plot(Cluster_xswitch_max(idxkmeans==6,1),'m.','MarkerSize',11)
            %         plot(Cluster_xswitch_max(idxkmeans==7,1),'c.','MarkerSize',11)
            %         plot(Ckmeans(:,1),'kx','MarkerSize',11,'LineWidth',1.5)
            
            %%
            % We find switch position by use xswitch_max because the algorithm
            % findchangepts find the switch prematurely. Then we use mean to add
            % statistic
            num_idxkmeans = unique(idxkmeans);
            histckmeans = [num_idxkmeans,histc(idxkmeans(:),num_idxkmeans)];
            [num_idxkmeans, switch_idxkeams] = sort(histckmeans(:,2),'descend');
            switch_idxkmeans = switch_idxkeams(1,1);
            xswitch_median = median(Cluster_xswitch_max(idxkmeans==switch_idxkmeans,1)); %switch position
%             xswitch_mean = mean(Cluster_xswitch_max(idxkmeans==switch_idxkmeans,1)); %switch position

%             a0_mean(ind_thick_unddata,leveldataprocessed,ind_eps) = xswitch_mean; %switch position calculate with mean
            a0_median(ind_thick_unddata,leveldataprocessed,ind_eps) = xswitch_median; %switch position calculate with median
            
            switch_index_spring = find(y_out_und(:,1)<=xswitch_median);
            new_data_y_spring = y_out_und(switch_index_spring,:);
            new_data_a_spring = a_out_und(switch_index_spring,:);
            switch_index_flight = find(y_out_und(:,1)>xswitch_median);
            new_data_y_flight = y_out_und(switch_index_flight,:);
            new_data_a_flight = a_out_und(switch_index_flight,:);
            
            y_out_spring = [y_out_spring ; new_data_y_spring];
            a_out_spring = [a_out_spring ; new_data_a_spring];
            
            y_out_flight = [y_out_flight ; new_data_y_flight];
            a_out_flight = [a_out_flight ; new_data_a_flight];
            
            %% spring mainfold
            % initialize a structure library
            ind_all_spring = []; II_Xi_spring = [];
            
            % stack all the "position" values in cluster ii into a vector
            X_spring = y_out_spring;
            % stack all the "velocity" values in cluster ii into a vector
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
            Thetalib_spring.Theta = ThetaN_spring;
            Thetalib_spring.normTheta = 1;
            Thetalib_spring.dx = XP_spring;
            Thetalib_spring.polyorder = polyorder_spring;
            Thetalib_spring.polyorderx = polyorder_spring_x;
            Thetalib_spring.polyordery = polyorder_spring_y;
            Thetalib_spring.usesine = usesine_spring;
            
            lambdavals_spring.numlambda = 20;
            lambdavals_spring.lambdastart = -5;
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
            X_flight = y_out_flight;
            % stack all the "velocity" values in cluster ii into a vector
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
            Thetalib_flight.Theta = ThetaN_flight;
            Thetalib_flight.normTheta = 1;
            Thetalib_flight.dx = XP_flight;
            Thetalib_flight.polyorder = polyorder_flight;
            Thetalib_flight.polyorderx = polyorder_flight_x;
            Thetalib_flight.polyordery = polyorder_flight_y;
            Thetalib_flight.usesine = usesine_flight;
            
            lambdavals_flight.numlambda = 20;
            lambdavals_flight.lambdastart = -5;
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
            
            tEnd(ind_thick_unddata,leveldataprocessed,ind_eps) = toc(tStart);      % pair 1: toc
            
            ind_thick_unddata = ind_thick_unddata+1;
        end
        ind_eps = ind_eps+1;
    end
end
%%
dateformatout = 'mmddyyyy';
%     save([datestr(now, dateformatout) 'Hopping_data_eps_zero.mat'])
save([datestr(now, dateformatout) 'Hopping_data_eps_thickunddata.mat'])
            

%% plot trun vs thick_undefined_data
%
% input a save location in the line below:
pathname = 'plot_suplearn_Hopper_trun';
% mkdir(pathname);

dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'trun_vs_thickunddata.fig'];
figure(1234)
for jj=1:3
    for kk=1:numleveldataprocessed %plot trun for each level data processed
        PercUndData = 100.*thick_unddata_vect./(max(y_out_flight(:,1))-min(y_out_spring(:,1)));
        plot(PercUndData,(tEnd(:,kk,jj))','LineWidth',1.5);
        hold on
        xlimits = [min(PercUndData) max(PercUndData)];
        ylimits = [0.9*min(min(min(min(tEnd)))) 1.1*max(max(max(max((tEnd)))))];
        axis([xlimits ylimits])
        l = get(gca, 'Children');
        xlabel('% Und. Data [ ]','Interpreter','latex')
        ylabel('$t_{run}$ [s]','Interpreter','latex')
        legend('num-i.c.eps1e-6 = 3','num-i.c.eps1e-6 = 4','num-i.c.eps1e-6 = 5',...
           'num-i.c.eps1e-4 = 3','num-i.c.eps1e-4 = 4','num-i.c.eps1e-4 = 5',...
           'num-i.c.eps1e-2 = 3','num-i.c.eps1e-2 = 4','num-i.c.eps1e-2 = 5')
        % set(gca,'FontSize',13)
        title('t_{run}  vs  thick_{undefined-data}')
        set(gca,'FontSize',9)
        set(l, 'linewidth', 1.5)
    end
end
savefig(plotname)

%% plot switch position vs thick_undefined_data
% plotname = [pathname datestr(now, dateformatout) 'xswitchMean_vs_thickunddata.fig'];
% figure(2345)
% for jj=1:3
%     for kk=1:numleveldataprocessed %plot trun for each level data processed
%         plot(thick_unddata_vect,nonzeros(a0_mean(:,kk,jj))','LineWidth',1.5);
%         hold on
%         xlabel('thick_{und-data} [ ]')
%         ylabel('x_{switch-mean} [L]')
%         %     axis tight
%         legend('num-i.c.eps1e-6 = 3','num-i.c.eps1e-6 = 4','num-i.c.eps1e-6 = 5',...
%            'num-i.c.eps1e-4 = 3','num-i.c.eps1e-4 = 4','num-i.c.eps1e-4 = 5',...
%            'num-i.c.eps1e-2 = 3','num-i.c.eps1e-2 = 4','num-i.c.eps1e-2 = 5')
%         % set(gca,'FontSize',13)
%         title('x_{switch-mean}  vs  thick_{undefined-data}')
%     end
% end
% savefig(plotname)
%%
plotname = [pathname datestr(now, dateformatout) 'xswitchMedian_vs_thickunddata.fig'];
figure(2346)
for jj=1:3
    for kk=1:numleveldataprocessed %plot trun for each level data processed
        PercUndData = 100.*thick_unddata_vect./(max(y_out_flight(:,1))-min(y_out_spring(:,1)));
%         plot(PercUndData,nonzeros(a0_median(:,kk,jj))','LineWidth',1.5);
        plot(PercUndData,(a0_median(:,kk,jj))'./thick_unddata_vect);%type plot2
        hold on
        xlimits = [min(PercUndData) max(PercUndData)];
%         ylimits = [0.9*min(min(min(min(abs(a0_mean-1))))) 1.02*max(max(max(max(abs(a0_mean-1)))))];
        ylimits = [0 1];
        axis([xlimits ylimits])
        l = get(gca, 'Children');
        xlabel('% Und. Data [ ]','Interpreter','latex')
        ylabel('$\epsilon_{xswitch-mean}$ [ ]','Interpreter','latex')
        %     axis tight
        legend('num-i.c.eps1e-6 = 3','num-i.c.eps1e-6 = 4','num-i.c.eps1e-6 = 5',...
           'num-i.c.eps1e-4 = 3','num-i.c.eps1e-4 = 4','num-i.c.eps1e-4 = 5',...
           'num-i.c.eps1e-2 = 3','num-i.c.eps1e-2 = 4','num-i.c.eps1e-2 = 5')
        % set(gca,'FontSize',13)
        title('x_{switch-median}  vs  thick_{undefined-data}')
        set(gca,'FontSize',9)
        set(l, 'linewidth', 1.5)
    end
end
% plot([0,max(PercUndData)],[1 1],'k--','LineWidth',1.5);
savefig(plotname)
