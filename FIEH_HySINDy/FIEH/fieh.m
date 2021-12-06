function dy = fieh(t,y,Ca,C,M,Yi,omega,kappa,stroke,mu,g,N,D,mu0,m,R,L)
%Ce (electrical damping coefficient) calculation
% syms n
% dphi_star = symsum(((y(1)-Yi*sin(y(3)))+n*(1/(N-1)))*...
%     (((2*((y(1)-Yi*sin(y(3)))+n*(1/(N-1)))^2)/((D^2/4+((y(1)-Yi*sin(y(3)))+n*(1/(N-1)))^2)^(5/2)))...
%     -(3/((D^2/4+((y(1)-Yi*sin(y(3)))+n*(1/(N-1)))^2)^(3/2))))...
%     ,n,-N/2,N/2);
% dphi_star = double(dphi_star);

dphi_star = 0;
for n = -N/2:N/2
  dphi_star = dphi_star + ((y(1)-Yi*sin(y(3)))+n*(L/(N-1)))*...
    (((2*((y(1)-Yi*sin(y(3)))+n*(L/(N-1)))^2)/(((D^2)/4+((y(1)-Yi*sin(y(3)))+n*(L/(N-1)))^2)^(5/2)))...
    -(3/(((D^2)/4+((y(1)-Yi*sin(y(3)))+n*(L/(N-1)))^2)^(3/2))));
end
% Ce = ((mu0*m/2*dphi_star)^2)/R;
Ce = 0; %TEST NICO


if (y(1)-Yi*sin(y(3)))>(stroke/2)
    % derivatives when FIEH is touching stop end 1
    dy=[
        y(2,:);
        -((Ca+C+Ce)/M).*(y(2,:)-Yi.*omega.*cos(y(3,:)))-kappa./M.*(abs(y(1,:)-Yi.*sin(y(3,:)))-stroke./2)-mu.*g.*(y(2,:)-Yi.*omega.*cos(y(3,:)))./(abs(y(2,:)-Yi.*omega.*cos(y(3,:))));
        omega;
        ];
elseif (y(1)-Yi*sin(y(3)))<(-stroke/2)
    % derivatives when FIEH is touching stop end 2
    dy=[
        y(2,:);
        -((Ca+C+Ce)/M).*(y(2,:)-Yi.*omega.*cos(y(3,:)))+kappa./M.*(abs(y(1,:)-Yi.*sin(y(3,:)))-stroke./2)-mu.*g.*(y(2,:)-Yi.*omega.*cos(y(3,:)))./(abs(y(2,:)-Yi.*omega.*cos(y(3,:))));
        omega;
        ];
else
    % derivatives when FIEH is free
    dy=[
        y(2,:);
        -((Ca+Ce)/M).*(y(2,:)-Yi.*omega.*cos(y(3,:)))-mu.*g.*(y(2,:)-Yi.*omega.*cos(y(3,:)))./(abs(y(2,:)-Yi.*omega.*cos(y(3,:))));
        omega;
        ];
end