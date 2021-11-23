% Run Hopping Model
function[data] = RUN_FIEH(p,thick_unddata)
% variables in data are: time, t, vertical position, y, vertical velocity, yp and
% acceleration
M = p.M; %[kg] mass
Ca = p.Ca; %[Ns/m] air damping coefficient
C = p.C; %[Ns/m] spring viscous damping coefficient
m = p.m; %[Am^2] magnetic moment
mu0 = p.mu0; %[N/A^2] vacum permeability
D = p.D; %[m] coil mean diameter
L = p.L; %[m] coil length
N = p.N; %total number of turns
R = p.R; %[Omega] coil resistance
kappa = p.kappa; %[N/m] FIEH stiffness
stroke = p.stroke; %[m] stroke
mu = p.mu; %[ ] friction coefficient
g = p.g; %[m/s^2] gravity acceleration
Yi = p.Yi; %[m] exitaction amplitude
omega = p.omega; %[Hz] exitaction frequency

tspan_in = p.tspan;
plottag = p.plottag;
yinitvec = p.yinitvec;
options = p.options;

dt = tspan_in(2)-tspan_in(1);
tend = tspan_in(end);

current_t = 0; % starting value of the current time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_out = [];
y_out = [];% y_out contain both vertical position and velocity
ydot_out = [];% a_out contain both vertical acceleration and velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Manifold_undefined
t_out_und = []; y_out_und = []; a_out_und = [];

%Manifold_1
t_out_1 = [];
y_out_1 = [];% y_out contain both vertical position and velocity
a_out_1 = [];% a_out contain both vertical acceleration and velocity
%Manifold_2
t_out_2 = [];
y_out_2 = [];% y_out contain both vertical position and velocity
a_out_2 = [];% a_out contain both vertical acceleration and velocity
%Manifold_3
t_out_3 = [];
y_out_3 = [];% y_out contain both vertical position and velocity
a_out_3 = [];% a_out contain both vertical acceleration and velocity

for mm = 1:size(yinitvec,1) %data creation(/<=>measurements) for each initial conditions
    mm
    yinit  = yinitvec(mm,:)
    current_t = 0; % starting value of the current time
    t = 0;
    tspan = tspan_in;
    
    % M*ypp = -(Ca+Ci+Ce)*(ydot-Yi*omega*cos(omega*t))...
    %-ke*(abs(y-Yi*sin(omega*t))-stroke/2)...
    %-mu*M*g*(ydot-Yi*omega*cos(omega*t))/(abs(y-Yi*sin(omega*t)))

    while current_t < tend
        if length(tspan)>1 % check that we haven't reached the end within error
            [t,y] = ode45(@(t,y) fieh(t,y,Ca,C,M,Yi,omega,kappa,stroke,mu,g,N,D,mu0,m,R,L),tspan,yinit,options);

%             a = fieh_impact1(t,y',Ca,C,M,Yi,omega*ones(1,numel(y(:,3))),kappa,stroke,mu,g,N,D,mu0,m,R)';
            
            figure(1)
            plot(t,y(:,1))
            hold on
            drawnow
            
            figure(2)
            plot(t,(y(:,1)-Yi*sin(y(:,3))))
            hold on
            drawnow
            %Manifold_1
            t_out = [t_out; t(1:end)];
            y_out = [y_out; y(1:end,:)];
            
            current_t = t_out(end);
            
            yinit = y(end,:);
            tspan = t(end):dt:tend+1;
        end
            
        
    end
   
    
end
z_out = [(y_out(:,1)-Yi.*sin(y_out(:,3))) , (y_out(:,2)-Yi.*omega.*cos(y_out(:,3)))];
zdot_out = diff(z_out,1,1); %Compute the first-order difference
z_out = z_out(1:end-1,:);
y_out = y_out(1:end-1,:);

%store only stable dynamic data
t_out = t_out(t_out>=8);
ind_stab = size(z_out,1)-size(t_out,1)+1;
z_out = z_out(ind_stab:end,:);
zdot_out = zdot_out(ind_stab:end,:);
y_ext = Yi.*sin(y_out(ind_stab:end,3));
ydot_ext = Yi.*omega.*cos(y_out(ind_stab:end,3));


%plot - check phase space
figure
plot(z_out(:,1),z_out(:,2),'k')
xlabel('$x$ [ ]','Interpreter','latex')
ylabel('$\dot{x}$ [ ]','Interpreter','latex')
title('phase space data')

% index undefined data
indx_und = find(z_out(:,1)<=(min(z_out(:,1))+thick_unddata) | z_out(:,1)>=(max(z_out(:,1))-thick_unddata));
% index defined data
indx_free = find(z_out(:,1)>(min(z_out(:,1))+thick_unddata) & z_out(:,1)<(max(z_out(:,1))-thick_unddata));

% save data clustered data

data.z_out_und = z_out(indx_und,:);
data.zdot_out_und = zdot_out(indx_und,:);
data.t_out_und = t_out(indx_und);

% data.z_out_impact1 = z_out(indx_impact1,:);
% data.zdot_out_impact1 = zdot_out(indx_impact1,:);
% data.tout_impact1 = t_out(indx_impact1);

data.z_out_free = z_out(indx_free,:);
data.zdot_out_free = zdot_out(indx_free,:);
data.t_out_free = t_out(indx_free);

% data.z_out_impact2 = z_out(indx_impact2,:);
% data.zdot_out_impact2 = zdot_out(indx_impact2,:);
% data.t_out_impact2 = t_out(indx_impact2);

%save No clustered data
data.z_out = z_out;
data.zdot_out = zdot_out;
data.t_out = t_out;

%save external "forcing"
data.y_ext = y_ext;
data.ydot_ext = ydot_ext;