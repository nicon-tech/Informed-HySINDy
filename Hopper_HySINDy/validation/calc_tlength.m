function [xswitch, abserror, abserror_avg, RMSE] = calc_tlength(xA, xB, val)

        % define vectors to compare, making sure they are the same
        % length, and calculate length
        [tlength, nterms] = size(xB);
        xAcomp =xA(1:tlength, :);
        xBcomp = xB;
        
        % calculate errors
        if length(xAcomp) == length(xB)
            error1 = xAcomp-xBcomp;
        else
            error('model and validation time-series are not the same length')
        end
        
        abserror = abs(error1);
        tvec = val.tA;
        
        for ll=1:nterms % for each term
            % find the time at which the error changes significantly
            % (switching point)
            if size(error1,1)>1
                switch_ind = findchangepts(error1(:,ll), 'Statistic','rms'); %Find the index where the root-mean-square level of the error1(:,ll) changes most significantly.
                %NOTE: ipt = findchangepts(x) returns the index at which
                %the mean of x changes most significantly.
                %If x is a vector with N elements, then findchangepts
                %partitions x into two regions, x(1:ipt-1) and x(ipt:N),
                %that minimize the sum of the residual (squared) error of
                %each region from its local mean.
                %If x is an M-by-N matrix, then findchangepts partitions x
                %into two regions, x(1:M,1:ipt-1) and x(1:M,ipt:N),
                %returning the column index that minimizes the sum of the
                %residual error of each region from its local M-dimensional mean.
                % => Altenatively: switch_ind = findchangepts(error1, 'Statistic','rms');
                % => In this way we can reduced the computational cost
                %(delete the for loop) and improve the code robustness
                % even if we find the switching point later than the other case.
            else
                switch_ind = [];
            end
            
            if any(size(switch_ind)<1)% check if it's empty
                if any(size(switch_ind) == 0)
                 switch_ind = length(abserror(:,ll));
                end
            end
%             xswitch(ll) = xAcomp(switch_ind,1);
            xswitch = xAcomp(switch_ind,1);

            if length(xAcomp) == length(xBcomp)
                abserror_avg(ll) = sum(abs(xAcomp(1:switch_ind,ll)-xBcomp(1:switch_ind,ll)))/switch_ind;
                RMSE(ll) = sqrt(sum(abs(xAcomp(1:switch_ind,ll)-xBcomp(1:switch_ind,ll)).^2)/switch_ind);
            else
                error('model and validation time-series are not the same length')
            end
        end
        
        
       