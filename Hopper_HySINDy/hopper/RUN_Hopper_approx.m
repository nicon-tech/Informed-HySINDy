%% Run Hopping HySINDy Model PLOT
function[data] = RUN_Hopper_approx(p,Xicomb_spring,Xicomb_flight,xswitch_max)
% variables in data are: time, t, vertical position, y, vertical velocity, yp and
% acceleration

Xicomb_spring = Xicomb_spring{1};%each time we analysed the result we need to search for better Xicomb
Xicomb_flight = Xicomb_flight{1};%each time we analysed the result we need to search for better Xicomb

kappa = p.kappa;
tspan_in = p.tspan;
plottag = p.plottag;
yinitvec = p.yinitvec;
options = p.options;

dt = tspan_in(2)-tspan_in(1);
tend = tspan_in(end);

current_t = 0; % starting value of the current time
t_out = []; y_out = []; a_out = [];

for mm = 1:size(yinitvec,1) %data simulation for each initial conditions
    mm
    yinit  = yinitvec(mm,:)
    current_t = 0; % starting value of the current time
    t_out = [t_out; 0];
    y_out = [y_out; yinit]; % y_out contain both vertical position and velocity
    a_out = [a_out; hopperspring(0,yinit',kappa)']; % a_out contain both vertical acceleration and velocity
    %NOTE_NN: a_out = [a_out; hopperspring(0,yinit',kappa)'] implica che si
    %parta in regime di spring
    t= 0;
    tspan = tspan_in;
    
    % m*xpp = -m*g - k(x-x0), x<=x0
    % m*xpp = -m*g,            x>x0
    while current_t < tend
        if length(tspan)>1 % check that we haven't reached the end within error
            if (yinit(1) <= xswitch_max) % if we are below xswitch_max we are springing
                disp('hopper springing')
                [t,y]=ode45(@(t,x)sparseGalerkin(t,x,Xicomb_spring,polyorder,usesine),...
                    tspan,yinit(:,kk),options);  % approximate
                
                figure(1)
                plot(t,y(:,1))
                xlabel('time [T]')
                ylabel('position [L]')
                % set(gca,'FontSize',13)
                title('position - time series')
                hold on
                drawnow
            elseif (yinit(1) > xswitch_max) % if we are over xswitch_max we are flying
                disp('hooper flying')
                [t,y]=ode45(@(t,x)sparseGalerkin(t,x,Xicomb_flight,polyorder,usesine),...
                    tspan,yinit(:,kk),options);  % approximate
                
                figure(1)
                plot(t,y(:,1))
                xlabel('time [T]')
                ylabel('position [L]')
                % set(gca,'FontSize',13)
                title('position - time series')
                hold on
                drawnow
            else
                disp('simulation ended early')
                current_t
                return
            end
            t_out = [t_out; t(2:end)];
            y_out = [y_out; y(2:end,:)];
            a_out = [a_out; a(2:end,:)];
            yinit = y(end,:);
            tspan = t(end):dt:tend+1;
            current_t = t_out(end);
            
        end
        
    end
   
    
end
 
% save data
data.yout = y_out;
data.aout = a_out;
data.tout = t_out;
