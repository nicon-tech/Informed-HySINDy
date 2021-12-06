% Run Hopping Model
function[data] = RUN_Hopper_training(p,thick_unddata)
% variables in data are: time, t, vertical position, y, vertical velocity, yp and
% acceleration

kappa = p.kappa;
tspan_in = p.tspan;
plottag = p.plottag;
yinitvec = p.yinitvec;
options = p.options;

dt = tspan_in(2)-tspan_in(1);
tend = tspan_in(end);

current_t = 0; % starting value of the current time
%Manifold_undefined
% t_out_und = []; y_out_und = []; a_out_und = [];
%Manifold_spring
t_out_spring = [];
y_out_spring = [];% y_out contain both vertical position and velocity
a_out_spring = [];% a_out contain both vertical acceleration and velocity
%Manifold_flight
t_out_flight = [];
y_out_flight = [];% y_out contain both vertical position and velocity
a_out_flight = [];% a_out contain both vertical acceleration and velocity

for mm = 1:size(yinitvec,1) %data creation(/<=>measurements) for each initial conditions
    mm
    yinit  = yinitvec(mm,:)
    current_t = 0; % starting value of the current time
    t = 0;
    tspan = tspan_in;
    
    % m*xpp = -m*g - k(x-x0), x<=x0
    % m*xpp = -m*g,            x>x0
    while current_t < tend
        if length(tspan)>1 % check that we haven't reached the end within error
            if abs(yinit(1)-1)<1e-12 % if we are within error of transition %Manifold_undefined
                if yinit(2)<0 % springing if velocity is negative
                   disp('hopper springing')
                    yinit(1) = 1-1e-12;
                    [t,y]=ode45(@(t,y) hopperspring(t,y,kappa),tspan,yinit,options);
                    tspring = t;
                    yspring = y;
                    a= hopperspring(t,y',kappa)';
                    aspring = a;
                    
                    figure(1)
                    plot(t,y(:,1))
                    hold on
                    drawnow
                    %Manifold_spring
                    t_out_spring = [t_out_spring; tspring(1:end)];
                    y_out_spring = [y_out_spring; yspring(1:end,:)];
                    a_out_spring = [a_out_spring; aspring(1:end,:)];
                    
                    current_t = t_out_spring(end);
                    
                else % flying
                    yinit(1) = 1+1e-12;
                    disp('hooper flying')
                    [t,y]=ode45(@(t,y) hopperflight(t,y),tspan,yinit,options);
                    tflight = t;
                    yflight = y;
                    a = hopperflight(t,y')';
                    aflight = a;
                 
                    plot(t,y(:,1))
                    hold on
                    drawnow
                    %Manifold_flight
                    t_out_flight = [t_out_flight; tflight(1:end)];
                    y_out_flight = [y_out_flight; yflight(1:end,:)];
                    a_out_flight = [a_out_flight; aflight(1:end,:)];
                    
                    current_t = t_out_flight(end);
                    
                end
            elseif (yinit(1)<1) % if we are below 1 we are springing
               disp('hopper springing')
                [t,y]=ode45(@(t,y) hopperspring(t,y,kappa),tspan,yinit,options);
                tspring = t;
                yspring = y;
                a= hopperspring(t,y',kappa)';
                aspring = a;
                
                plot(t,y(:,1))
                hold on
                drawnow
                %Manifold_spring
                t_out_spring = [t_out_spring; tspring(1:end)];
                y_out_spring = [y_out_spring; yspring(1:end,:)];
                a_out_spring = [a_out_spring; aspring(1:end,:)];
                
                current_t = t_out_spring(end);
                
            elseif (yinit(1) > 1) % if we are over 1 we are flying
                disp('hooper flying')
                [t,y]=ode45(@(t,y) hopperflight(t,y),tspan,yinit,options);
                tflight = t;
                yflight = y;
                a = hopperflight(t,y')';
                aflight = a;
                
                plot(t,y(:,1))
                hold on
                drawnow
                %Manifold_flight
                t_out_flight = [t_out_flight; tflight(1:end)];
                y_out_flight = [y_out_flight; yflight(1:end,:)];
                a_out_flight = [a_out_flight; aflight(1:end,:)];
                
                current_t = t_out_flight(end);
                
            else
                disp('simulation ended early')
                current_t
                return
            end                                
                yinit = y(end,:);
                tspan = t(end):dt:tend+1;

         end
        
    end
   
    
end

% index undefined data
indx_und_spring = find((y_out_spring(:,1)-1)<=0 & (1-y_out_spring(:,1))<(thick_unddata/2));
indx_und_flight = find((y_out_flight(:,1)-1)>=0 & (y_out_flight(:,1)-1)<(thick_unddata/2));
% index defined data
indx_spring = find((y_out_spring(:,1)-1)<=0 & (1-y_out_spring(:,1))>(thick_unddata/2));
indx_flight = find((y_out_flight(:,1)-1)>=0 & (y_out_flight(:,1)-1)>(thick_unddata/2));

% save data

data.yout_und = [y_out_spring(indx_und_spring,:) ; y_out_flight(indx_und_flight,:)];
data.aout_und = [a_out_spring(indx_und_spring,:) ; a_out_flight(indx_und_flight,:)];
data.tout_und = [t_out_spring(indx_und_spring,:) ; t_out_flight(indx_und_flight,:)];

data.yout_spring = y_out_spring(indx_spring,:);
data.aout_spring = a_out_spring(indx_spring,:);
data.tout_spring = t_out_spring(indx_spring,:);

data.yout_flight = y_out_flight(indx_flight,:);
data.aout_flight = a_out_flight(indx_flight,:);
data.tout_flight = t_out_flight(indx_flight,:);
% 
% figure(2); clf
% plot(y_out_spring(:,1),y_out_spring(:,2),'*')
% hold on
% plot(y_out_flight(:,1),y_out_flight(:,2),'*')
% plot(yout_und(:,1),yout_und(:,2))



