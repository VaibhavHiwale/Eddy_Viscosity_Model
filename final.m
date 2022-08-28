clc; close all; clear all;

% Read DNS data [half-channel is given (till centerline)]
load y_dns.dat
load u_dns.dat
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat
load uv_dns.dat
load dns_data.dat
load d_y.dat
%------------------------------- PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu=1/395;
ustar=1;
rho = 1;
kappa=0.41;

% 1-equation model
c_mu=0.09;
sigma_k=1;

% Residual error limit
residue_limit = 10^(-6);
% Under relaxation factor
urf = 0.5;

%Initial conditions for old & new variables (U, dUdy, k, eps, nu_t, residue, ...)
U = u_dns;
U1 = zeros(97,1);
dUdy = zeros(97,1);
k = ones(97,1);
k1 = ones(97,1);
eps = ones(97,1);
nu_t = ones(97,1);
residue = 1;
r1=1;
r2=1;
pk = zeros(97,1);
ur= ones(97,1);
kr= zeros(97,1);
v_d = ones(97,1);
t_d = zeros(97,1);

%Boundary condition (U, k, eps, nu_t, ...)
U(97,1) = U(96,1); %mid
k(97,1) = k(96,1); %mid
U(1,1) = 0; %wall
k(1,1) = 0; %wall

%dy is created using Excel from DNS grid data
%------------------------------- CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate
while residue > residue_limit
 
    %Boundary condition (U, k, eps, nu_t, ...)
    U(97,1) = U(96,1); %mid
    k(97,1) = k(96,1); %mid
    U(1,1) = 0; %wall
    k(1,1) = 0; %wall


    %Compute eddy viscosity

for i=2:96
   
    nu_t(i,1) = 1/(k(i,1)^0.5);

end

       %Compute U
 for i=2:96
  
    U1(i,1) = 0.5*U(i+1,1)+0.5*U(i-1,1)+0.5*((d_y(i-1,1)^2 + d_y(i,1)^2)/(nu+nu_t(i,1)));

 end

    %Compute dUdy
 for i=2:96

    dUdy(i,1) = (U(i+1,1)-U(i-1,1))/(2*d_y(i-1,1));

 end

    %Compute Pk
for i=2:96

    pk(i,1) = nu_t(i,1)*(dUdy(i,1)^2);

end

    %Compute k
for i=2:96

    k1(i,1) = 0.5*(((d_y(i-1,1)^2 + d_y(i,1)^2)*(pk(i,1)-eps(i,1))/(nu+nu_t(i,1)))+2*urf*k(i+1,1)+2*urf*k(i-1,1));

end

    %Compute epsilon
for i=2:96

    eps(i,1)=0.09*k(i,1)^1.5;

end

    %Compute residue
 
for i=2:96
    
        ur(i,1) = abs(U(i,1)-U1(i,1));
        kr(i,1) = abs(k(i,1)-k1(i,1));
end

%maximum value of velocity difference
    r1 = ur(2,1);

for i=2:96
if r1<= ur(i,1)
    r1=ur(i,1);
end
end

%maximum value of kinetic energy difference
     r2 = kr(2,1);

for i=2:96
if r2<= kr(i,1)
    r2=kr(i,1);
end
end

%maximum value of residue from u and k
residue = max(r1,r2);

  
   %swap
for i=2:96
  
       U(i,1)=U1(i,1);
       k(i,1)=k1(i,1);

end
    
    %print k value
    disp(k(50));

end

 %Compute viscous diffusion
for i=2:96

      v_d (i,1) = nu*(2*(k(i+1,1)+2*k(i-1,1)-4*k(i,1))/(d_y(i-1,1)^2 + d_y(i,1)^2));

end    

%-------------------------------- PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_dns=0.5*(u2_dns+v2_dns+w2_dns);
eps_dns=dns_data(:,2)*ustar^4/nu; % eps is normalized by ustar^4/nu

figure(1)
hold on
plot(u_dns,y_dns,'bo');
plot(U,y_dns,'ro');
xlabel('U'); ylabel('y/h'); title('U-velocity');
legend('DNS','Prandtls 1-Eqn Model','Best'); legend boxoff;

figure(2)
hold on
plot(y_dns,k_dns,'bo');
plot(y_dns,k,'ro');
xlabel('y/h'); ylabel('k'); title('Turbulence kinetic energy');
legend('DNS','Prandtls 1-Eqn Model'); legend boxoff;

figure(3)
hold on
plot(y_dns,eps_dns,'bo')
plot(y_dns,eps,'ro')
xlabel('y/h'); ylabel('\epsilon'); title('Dissipation rate of k');
legend('DNS','Prandtls 1-Eqn Model'); legend boxoff;

figure(4)
plot(y_dns,-uv_dns,'bo')
xlabel('y/h'); ylabel('-<uv>'); title('Turbulence shear stress');
legend('DNS'); legend boxoff

figure(5)
plot(y_dns,nu_t,'ro')
xlabel('y/h'); ylabel('nu_t'); title('Turbulent Viscosity');
legend('Prandtls 1-Eqn Model'); legend boxoff

figure(6)
hold on
plot(y_dns,dns_data(:,3)/nu,'bo');
plot(y_dns,pk/nu,'ro')
xlabel('y/h'); ylabel('P_k'); title('Production rate of k');
legend('DNS','1-Eqn Model'); legend boxoff;

figure(7)
plot(y_dns,dns_data(:,5)/nu,'bo');
xlabel('y/h'); title('Turbulent diffusion of k');
legend('DNS'); legend boxoff

figure(8)
hold on
plot(y_dns,dns_data(:,6)/nu,'bo');
plot(y_dns, v_d/nu,'ro');
xlabel('y/h'); ylabel('V_k');title('Viscous diffusion of k');
legend('DNS','1-Eqn Model'); legend boxoff

figure(9)
hold on
plot(y_dns, v_d/nu,'ro');
plot(y_dns, eps,'bo')
plot(y_dns, pk/nu,'go')
plot(y_dns, dns_data(:,6)/nu,'rx');
plot(y_dns, eps_dns,'bx')
plot(y_dns, dns_data(:,3)/nu,'gx')
xlabel('y/h'); title('Budget of k');
legend('viscous diffusion of k','dissipation rate of k', 'production rate of k','DNS viscous diffusion of k','DNS dissipation rate of k', 'DNS production rate of k'); legend boxoff

figure(10)
hold on
plot(y_dns, v_d,'ro');
plot(y_dns, eps,'bo')
plot(y_dns, pk,'go')
xlabel('y/h'); title('Budget of k Prandtl 1-Eqn Model');
legend('viscous diffusion of k','dissipation rate of k', 'production rate of k'); legend boxoff

figure(11)
hold on
plot(y_dns, dns_data(:,6)/nu,'rx');
plot(y_dns, eps_dns,'bx')
plot(y_dns, dns_data(:,3)/nu,'gx')
xlabel('y/h'); title('Budget of k DNS');
legend('DNS viscous diffusion of k','DNS dissipation rate of k', 'DNS production rate of k'); legend boxoff