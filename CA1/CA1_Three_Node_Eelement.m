close all; clear; clc
% Physical Parameters
R_i = 10; R_o = 20;
P = 10; E = 1000; nu = 0.3;
lambda = E*nu/((1 + nu)*(1 - 2*nu));
mu = E/(2*(1 + nu));
% Spices (2 Node Linear)
C = [lambda + 2*mu, mu, mu; mu, lambda + 2*mu, mu; mu, mu, lambda + 2*mu];
x_iso = @(x_e,xi) [0.5*xi*(xi-1),(1-xi^2),0.5*xi*(xi+1)]*x_e;  % Isoparametric Mapping
B = @(h,xi,x_iso) [2*(xi-0.5)/h, -4*xi/h, 2*(xi+0.5)/h; xi*(xi-1)/(2*x_iso), (1-xi^2)/x_iso, xi*(xi+1)/(2*x_iso); ...
    xi*(xi-1)/(2*x_iso), (1-xi^2)/x_iso, xi*(xi+1)/(2*x_iso)];
% Three Point Gauss Guadrature
xi_1 = -sqrt(3/5); xi_2 = 0; xi_3 = sqrt(3/5);
w1 = 5/9; w2 = 8/9; w3 = 5/9;
% Finite Element Methods
d = zeros(42+41,3); 
sigma_rr = zeros(1,3);
sigma_tt = zeros(1,3);
nof_element = [11, 21, 41];
h = (R_o - R_i)./nof_element;
str = ["ko","c*","r."];
figure; hold on; title('3 Node')
for m = 1 : 3
    % Finite Element Discretization
    N_el = nof_element(m);
    x = linspace(R_i, R_o, N_el+1);  % End points of elements
    x_plot = linspace(R_i, R_o, N_el+1+N_el);
    % Prelocation for Global Assembly
    K = zeros(N_el+1+N_el, N_el+1+N_el);
    f = zeros(N_el+1+N_el, 1);

    for i = 1 : N_el
        x1 = x(i); x3 = x(i+1); x2 = 0.5*(x1+x3);
        x_e = [x1;x2;x3];
        x_xi1 = x_iso(x_e,xi_1); x_xi2 = x_iso(x_e,xi_2); x_xi3 = x_iso(x_e,xi_3);
        K_e = 4*pi*h(m)/2 * (w1*B(h(m),xi_1,x_xi1).'*C*B(h(m),xi_1,x_xi1)*x_xi1^2 ...
            + w2*B(h(m),xi_2,x_xi2).'*C*B(h(m),xi_2,x_xi2)*x_xi2^2 ...
            + w3*B(h(m),xi_3,x_xi3).'*C*B(h(m),xi_3,x_xi3)*x_xi3^2);
        K(2*(i-1)+1 : 2*i+1, 2*(i-1)+1 : 2*i+1) = K(2*(i-1)+1:2*i+1,2*(i-1)+1:2*i+1) + K_e;
    end
    f(1:3,1) = 4*pi*R_i^2*P*[1;0;0];
    d(1:N_el+1+N_el, m) = K \ f;
    plot(x_plot, d(1:N_el+1+N_el,m), str(m))
end

% Exact Solution
x_exact = linspace(R_i,R_o,42+41);
sigma_rr_exact = P*R_i^3*(R_o^3-x_exact.^3)./(x_exact.^3*(R_i^3-R_o^3));
sigma_tt_exact = -P*R_i^3*(2*x_exact.^3 + R_o^3)./(2*x_exact.^3*(R_i^3-R_o^3));
eps = C \ [sigma_rr_exact; sigma_tt_exact; sigma_tt_exact];
eps_tt = eps(2,:);
u_exact = x_exact.*eps_tt;
plot(x_exact, u_exact)
lgd = legend('$$u^{h}|N_{el} = 11$$','$$u^{h}|N_{el} = 21$$','$$u^{h}|N_{el} = 41$$','u');
set(lgd, 'Interpreter','latex');

% Error of Solution at Inner Surface and Outer Surfaec
error_inner = abs(d(1,:) - u_exact(1));
error_outer = abs([d(12+11,1),d(22+21,2),d(42+41,3)] - u_exact(end));
h = [(R_o-R_i)/11,(R_o-R_i)/21,(R_o-R_i)/41 ];
figure; loglog(h, error_inner,'*-'); 
tit = title('Error of $$u_{r}$$ at Inner Surface by 3-node elements'); set(tit,'Interpreter','latex');
xlabel('mesh size'); ylab = ylabel('$$|u_{r}-u_{r}^{h}|$$'); set(ylab, "Interpreter",'latex')
slope = (log(error_inner(3))-log(error_inner(1)))/(log(h(3))-log(h(1)))

figure; loglog(h,error_outer,'o-');
tit = title('Error of $$u_{r}$$ at Outer Surface by 3-node elements'); set(tit,'Interpreter','latex');
xlabel('mesh size'); ylab = ylabel('$$|u_{r}-u_{r}^{h}|$$'); set(ylab, "Interpreter",'latex')

% Error of Sigma_rr, Sigma_thetatheta, and their exact value at middle point
% Middle Point in element corresponds to xi = 0, then by isoparametric mapping, x = x_2 at mid-point
sigma = zeros(3,3);
R_mid = (R_i+R_o)/2;
for i = 1 : 3
    N_el = nof_element(i);
    d_e = d(N_el : N_el+2, i);
    B_mid = B(h(i),0,R_mid);
    sigma(:,i) = C*B_mid*d_e;
end
sigma_rr_ex = P*R_i^3*(R_o^3-R_mid^3)/(R_mid^3*(R_i^3-R_o^3));
sigma_tt_ex = -P*R_i^3*(2*R_mid^3+R_o^3)/(2*R_mid^3*(R_i^3-R_o^3));
error_sigma_rr = abs(sigma_rr_ex - sigma(1,:));
error_sigma_tt = abs(sigma_tt_ex - sigma(2,:));

figure; loglog(h,error_sigma_rr,'*-');  
tit = title('Error of $$\sigma_{rr}$$ At Middle Point by 3-node elements'); set(tit,'Interpreter','latex');
xlabel('mesh size'); ylab = ylabel('$$|\sigma_{rr}-\sigma_{rr}^{h}|$$'); set(ylab, "Interpreter",'latex')

figure; loglog(h,error_sigma_tt,'o-'); 
tit = title('Error of $$\sigma_{\theta\theta}$$ At Middle Point by 3-node elements'); set(tit,'Interpreter','latex');
xlabel('mesh size'); ylab = ylabel('$$|u_{r}-u_{r}^{h}|$$'); set(ylab, "Interpreter",'latex');

slope = (log(error_sigma_tt(3))-log(error_sigma_tt(1)))/(log(h(3))-log(h(1)))
slope = (log(error_sigma_rr(3))-log(error_sigma_rr(1)))/(log(h(3))-log(h(1)))