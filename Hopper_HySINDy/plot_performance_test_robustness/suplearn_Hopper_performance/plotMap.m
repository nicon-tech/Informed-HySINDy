%% take averages across the instances  (3rd variable)
pathname = 'plotMap';
% mkdir(pathname);
% changeplot
%% plot trun vs thick_undefined_data & k_knn
% % input a save location in the line below:
pathname = 'plot_suplearn_Hopper_performance';
% % mkdir(pathname);
% 
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'trun_vs_PercUndData.fig'];
figure(1234)
for ii=1:3 %run for each noise level
    for jj=1:numleveldataprocessed %plot trun for each level data processed
        subplot(3,3,jj+3*(ii-1))
        for kk=1:numel(k_vect)
            PercUndData = 100.*thick_unddata_vect./(max(y_out_flight(:,1))-min(y_out_spring(:,1)));
            plot(PercUndData,(tEnd(:,kk,jj,ii)));
            hold on
            xlimits = [min(PercUndData) max(PercUndData)];
            ylimits = [0.9*min(min(min(min(tEnd)))) 1.1*max(max(max(max((tEnd)))))];
            axis([xlimits ylimits])
            l = get(gca, 'Children');
            xlabel('% Und. Data [ ]','Interpreter','latex')
            ylabel('$t_{run}$ [s]','Interpreter','latex')
            legend('k = 6','k = 7','k = 8','k = 9','k = 10')
            txt = ['i.c. ' int2str(jj) ' eps1e- ' int2str(6-2*(ii-1))];
            %         subtitle(txt)
            %         title('t_{run}  vs  thick_{undefined-data}')
            title(txt)
%             subtitle(['error' int2str(ii)])
            set(gca,'FontSize',9)
            set(l, 'linewidth', 1.5)
        end
    end
end
savefig(plotname)

%% Plots for paper

dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'plotMap_trunUnddataEps.fig'];
f = figure(1);
for jj=1:numleveldataprocessed
    subplot(3,1,jj)
    PercUndData = 100.*thick_unddata_vect./(max(y_out_flight(:,1))-min(y_out_spring(:,1)));
%             plot(PercUndData,(tEnd(:,kk,jj,ii)));
%     [C,h]=contourf(PercUndData, [1e-6 1e-4 1e-2], tEnd(jj+23*(jj-1):jj+23*jj,:,jj)');
    tEndCotourf = reshape(tEnd(:,5,jj,:),[24,3]);
    [C,h]=contourf(PercUndData, [1e-6 1e-4 1e-2], tEndCotourf');
    xlabel('% Und. Data [ ]','Interpreter','latex')
    ylabel('$noise level$ [ ]','Interpreter','latex')
    title(['i.c.' int2str(jj) ' n°10Knn'],'Interpreter','latex')
%     f.Children.XScale = 'log'; f.Children.YScale = 'log';
%     colormapnew
    c = colorbar;
    c.Label.String = ('t_{run} [s]');
    c.Label.FontSize = 11;
    f.CurrentAxes.FontSize = 9;
    f.CurrentAxes.LineWidth = 0.8;
    f.PaperUnits = 'centimeters';
    f.CurrentAxes.YTick = [1e-6 1e-4 1e-2];
    colormap(hot);
    % colorbar('off')
    % h.LevelStep = 0.2;
    % axis([1e-4 10 10 1.5e4])
    % f.CurrentAxes.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
    
    
end
savefig(plotname);
%% Plots for paper

dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'plotMap_xswitch_medianUnddataEps.fig'];
f = figure(1);
for jj=1:numleveldataprocessed
    subplot(3,1,jj)
    PercUndData = 100.*thick_unddata_vect./(max(y_out_flight(:,1))-min(y_out_spring(:,1)));
%             plot(PercUndData,(tEnd(:,kk,jj,ii)));
%     [C,h]=contourf(PercUndData, [1e-6 1e-4 1e-2], tEnd(jj+23*(jj-1):jj+23*jj,:,jj)');
    a0_medianCotourf = reshape((abs(a0_median(:,5,jj,:)-1)),[24,3]);
    a0_medianCotourf = a0_medianCotourf'./thick_unddata_vect;
    [C,h]=contourf(PercUndData, [1e-6 1e-4 1e-2], a0_medianCotourf);
    xlabel('% Und. Data [ ]','Interpreter','latex')
    ylabel('$noise level$ [ ]','Interpreter','latex')
    title(['i.c.' int2str(jj) ' n°10Knn'],'Interpreter','latex')
%     f.Children.XScale = 'log'; f.Children.YScale = 'log';
%     colormapnew
    c = colorbar;
    c.Label.String = ('\epsilon_{xswitch-median} [ ]');
    c.Label.FontSize = 11;
    f.CurrentAxes.FontSize = 9;
    f.CurrentAxes.LineWidth = 0.8;
    f.PaperUnits = 'centimeters';
    f.CurrentAxes.YTick = [1e-6 1e-4 1e-2];
    colormap(hot);
    % colorbar('off')
    % h.LevelStep = 0.2;
    % axis([1e-4 10 10 1.5e4])
    % f.CurrentAxes.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
    
    
end
savefig(plotname);
%% Plots for paper

dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'plotMap_xswitch_meanUnddataEps.fig'];
f = figure(1);
for jj=1:numleveldataprocessed
    subplot(3,1,jj)
    PercUndData = 100.*thick_unddata_vect./(max(y_out_flight(:,1))-min(y_out_spring(:,1)));
%             plot(PercUndData,(tEnd(:,kk,jj,ii)));
%     [C,h]=contourf(PercUndData, [1e-6 1e-4 1e-2], tEnd(jj+23*(jj-1):jj+23*jj,:,jj)');
    a0_meanCotourf = reshape((abs(a0_mean(:,5,jj,:)-1)),[24,3]);
    a0_meanCotourf = a0_meanCotourf'./thick_unddata_vect;
    [C,h]=contourf(PercUndData, [1e-6 1e-4 1e-2], a0_meanCotourf);
    xlabel('% Und. Data [ ]','Interpreter','latex')
    ylabel('$noise level$ [ ]','Interpreter','latex')
    title(['i.c.' int2str(jj) ' n°10Knn'],'Interpreter','latex')
%     f.Children.XScale = 'log'; f.Children.YScale = 'log';
%     colormapnew
    c = colorbar;
    c.Label.String = ('\epsilon_{xswitch-mean} [ ]');
    c.Label.FontSize = 11;
    f.CurrentAxes.FontSize = 9;
    f.CurrentAxes.LineWidth = 0.8;
    f.PaperUnits = 'centimeters';
    f.CurrentAxes.YTick = [1e-6 1e-4 1e-2];
    colormap(hot);
    % colorbar('off')
    % h.LevelStep = 0.2;
    % axis([1e-4 10 10 1.5e4])
    % f.CurrentAxes.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
    
    
end
savefig(plotname);



