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

figure(9)
plot(y,nu_t,'ro')
xlabel('y'); ylabel('nu_t_RSM'); title('Turbulent Viscosity');