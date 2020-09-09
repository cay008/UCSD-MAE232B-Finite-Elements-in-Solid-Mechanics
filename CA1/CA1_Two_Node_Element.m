close all; clear; clc
% Physical Parameters
R_i = 10; R_o = 20;
P = 10; E = 1000; nu = 0.3;
lambda = E*nu/((1 + nu)*(1 - 2*nu));
mu = E/(2*(1 + nu));
% Spices (2 Node Linear)
x_iso = @(x_e,xi) [0.5*(1-xi), 0.5*(1+xi)]*x_e;  % Isoparametric Mapping, plug this as the last argument in B
C = [lambda + 2*mu, mu, mu; mu, lambda + 2*mu, mu; mu, mu, lambda + 2*mu];
B = @(h,xi,x_iso) [-1/h, 1/h; (1-xi)/(2*x_iso), (1+xi)/(2*x_iso); (1-xi)/(2*x_iso), (1+xi)/(2*x_iso)];
% Two Point Gauss Quadrature
xi_1 = -sqrt(1/3); xi_2 = sqrt(1/3);
w1 = 1; w2 = 1;
% Finite Element Methods
d = zeros(42,3);
sigma_rr = zeros(1,3);
sigma_tt = zeros(1,3);
nof_element = [11, 21, 41];
h = (R_o - R_i)./nof_element;
figure; hold on; title('2 Node');
str = ["ko","c*","r."];
for m = 1 : 3
    % Finite Element Discretization
    N_el = nof_element(m);
    x = linspace(R_i, R_o, N_el+1);
    % Prelocation for Global Assembly
    K = zeros(N_el+1, N_el+1);
    f = zeros(N_el+1, 1);
    
    for i = 1 : N_el
        x1 = x(i); x2 = x(i+1);
        x_e = [x1; x2];
        x_xi1 = x_iso(x_e,xi_1); x_xi2 = x_iso(x_e,xi_2);  % x at gauss integration point
        K_e = 4*pi*h(m)/2 * (w1*B(h(m),xi_1,x_xi1).'*C*B(h(m),xi_1,x_xi1)*x_xi1^2 ...
            + w2*B(h(m),xi_2,x_xi2).'*C*B(h(m),xi_2,x_xi2)*x_xi2^2);
        K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + K_e;
    end
    f(1:2,1) = 4*pi*R_i^2*P*[1;0];
    d(1:N_el+1, m) = K \ f;
    plot(x, d(1:N_el+1,m), str(m))
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

% Error of Solution at Inner Surface and Outer Suface
error_inner = abs(d(1,:) - u_exact(1));
error_outer = abs([d(12,1),d(22,2),d(42,3)] - u_exact(end));
figure; loglog(h, error_inner,'*-'); 
tit = title('Error of $$u_{r}$$ at Inner Surface by 2-node elements'); set(tit,'Interpreter','latex');
xlabel('mesh size'); ylab = ylabel('$$|u_{r}-u_{r}^{h}|$$'); set(ylab, "Interpreter",'latex')
slope = (log(error_inner(3))-log(error_inner(1)))/(log(h(3))-log(h(1)))

figure; loglog(h,error_outer,'o-');
tit = title('Error of $$u_{r}$$ at Outer Surface by 2-node elements'); set(tit,'Interpreter','latex');
xlabel('mesh size'); ylab = ylabel('$$|u_{r}-u_{r}^{h}|$$'); set(ylab, "Interpreter",'latex')
slope = (log(error_outer(3))-log(error_outer(1)))/(log(h(3))-log(h(1)))

% Error of Sigma_rr, Sigma_thetatheta, and their exact value at middle point
sigma = zeros(3,3);
R_mid = (R_i+R_o)/2;
for i = 1 : 3
    N_el = nof_element(i);
    d_e = d(0.5*(N_el+1):0.5*(N_el+1)+1,i);
    B_mid = B(h(i),0,R_mid);  % xi = 0 at the middle of local element
    sigma(:,i) = C*B_mid*d_e;
end
sigma_rr_ex = P*R_i^3*(R_o^3-R_mid^3)/(R_mid^3*(R_i^3-R_o^3));
sigma_tt_ex = -P*R_i^3*(2*R_mid^3+R_o^3)/(2*R_mid^3*(R_i^3-R_o^3));
error_sigma_rr = abs(sigma_rr_ex - sigma(1,:));
error_sigma_tt = abs(sigma_tt_ex - sigma(2,:));

figure; loglog(h,error_sigma_rr,'*-');  
tit = title('Error of $$\sigma_{rr}$$ At Middle Point by 2-node elements'); set(tit,'Interpreter','latex');
xlabel('mesh size'); ylab = ylabel('$$|\sigma_{rr}-\sigma_{rr}^{h}|$$'); set(ylab, "Interpreter",'latex')

figure; loglog(h,error_sigma_tt,'o-'); 
tit = title('Error of $$\sigma_{\theta\theta}$$ At Middle Point by 2-node elements'); set(tit,'Interpreter','latex');
xlabel('mesh size'); ylab = ylabel('$$|u_{r}-u_{r}^{h}|$$'); set(ylab, "Interpreter",'latex');

slop = (log(error_sigma_rr(3))-log(error_sigma_rr(1)))/(log(h(3))-log(h(1)))

