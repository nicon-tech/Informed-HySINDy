%% training data generation
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
Yi = [0.01; 0.05] ; %0.01; %[m] exitaction amplitude
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

tend = 11;
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
%
%plot - check phase space
% figure
plot(z_out(:,1),z_out(:,2),'-*b')
% xlabel('$x$ [ ]','Interpreter','latex')
% ylabel('$\dot{x}$ [ ]','Interpreter','latex')
% title('phase space data')