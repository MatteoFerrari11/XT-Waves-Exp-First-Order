% Accompanying code for the Paper:
% "Inf-sup stable space--time discretization of the wave equation based on a first-order-in-time variational formulation"
%
% Authors: M. Ferrari and I. Perugia and E. Zampa
%
% This code computes an approximate solution of the wave problem:
%
%     d^2_t U(x,t) - dx (c(x) dx U(x,t)) = F(x,t),  x,t in Q_T = (0,L)x(0,T)
%     U(x,0) = U0(x),  dt U(x,0) = V0(x),           x   in (0,L)
%     U(0,t) = 0,  U(L,t) = 0,                      t   in (0,T)
%
% in a first order formulation
%
%     d_t V(x,t) - dx (c(x) dx U(x,t)) = F(x,t),    x,t in Q_T = (0,L)x(0,T)
%     d_t U(x,t) - V(x,t) = 0,                      x,t in Q_T = (0,L)x(0,T)
%     U(x,0) = U0(x),  V(x,0) = V0(x),              x in (0,L)
%     U(t,0) = V(t,0) = 0,  V(t,L) = U(t,L) = 0,    t in (0,T)
%
% using the variational formulation:
%
% find U in H^1(0,T;H_0^1(0,L)),  V in H^1(0,T;L^2(0,L)) with U(x,0)=V(x,0)=0
% such that
%
%     (d_t V, lambda)_{L_e^2(Q_T)} + (c dx U, dx lambda)_{L_e^2(Q_T)} = (F,lambda)_{L^2_e(Q_T)} - (c dx U0, dx lambda)_{L^2_e(Q_T)}
%     (d_t U, chi)_{L_e^2(Q_T)} - (V, chi)_{L_e^2(Q_T)} = (V0,chi)_{L^2_e(Q_T)
%
% for all lambda in L^2(0,T;H_0^1(0,L)) and chi in L^2(Q_T)
%
% Here, (U,V)_{L^2_e(Q_T)} = (U,V e^{-./T})_{L^2(Q_T)} 
%
% The code employs arbitrary regularity splines in time and space
% for the numerical solution

% and calculates errors in the norm L^2 of U, V
% Note: this code is far from being optimized

clc
clear
close all
format long

p_x = 3;
p_t = 3;
T = 1;
r_x = 2;
r_t = 2;

N_plot = 100;
Nq = 16;
N1 = 1;
N2 = 3;

h_x = zeros(1,N2-N1+1);
err_U = zeros(1,N2-N1+1);
err_V = zeros(1,N2-N1+1);
err_norm_call_V = zeros(1,N2-N1+1);


% SMOOTH EXAMPLE
a = 0;
b = 1;
U = @(t,x) sin(t*pi*5/4).^2*sin(x*pi) + ones(size(t))*sin(x*pi);
U0 = @(x)  U(0,x);
dU0 = @(x) pi*cos(x*pi);
A = @(t) sin(5/4*pi*t).^2+1;
F = @(t,x) ones(size(t))*(pi*cos(x.*pi).*-2.0)...
    -ones(size(t))*(pi.^2.*sin(x.*pi).*(2.5e+1./8.0))...
    +pi.*cos(t.*pi.*(5.0./4.0)).^2*cos(x.*pi)...
    +ones(size(t))*(x.*pi.^2.*sin(x.*pi).*2.0)...
    +pi.^2.*cos(t.*pi.*(5.0./4.0)).^2*(sin(x.*pi).*(2.5e+1./4.0))...
    -pi.^2.*cos(t.*pi.*(5.0./4.0)).^2*(sin(x.*pi).*x)...
    +A(t)*(pi^2*sin(pi*x));
V0 = @(x) 0*x;
c = @(x) x+1;
dtU = @(t,x) pi*(cos(t*pi*5/4).*sin(t*pi*5/4))*sin(x*pi)*2.5;
dxU = @(t,x) pi*ones(size(t))*cos(x*pi)+pi*sin(t*pi*5/4).^2*cos(x*pi);

% SINGULAR EXAMPLE
% a = -1.5;
% b = 1.5;
% w = @(x) exp(-20*(x-0.1).^2) - exp(-20*(x+0.1).^2);
% dw = @(x) exp(-20*(x+0.1).^2).*(40*x+4) - exp(-20*(x-0.1).^2).*(40*x-4);
% U = @(t,x) w(ones(size(t))*x-t*ones(size(x))+1).*(ones(size(t))*x-t*ones(size(x))+1>0);
% U0 = @(x) w(x+1).*(x+1>0);
% F = @(t,x) 0*t*x;
% V0 = @(x) -dw(x+1).*(x+1>0);
% dU0 = @(x) dw(x+1).*(x+1>0);
% c = @(x)  1+0*x;
% dtU =  @(t,x) -dw(ones(size(t))*x-t*ones(size(x))+1).*(ones(size(t))*x-t*ones(size(x))+1>0);

cont = 0;

for ii = N1:N2

    tic

    cont = cont + 1;
    N_x = 2^(ii+1);
    N_t = 2^(ii+1);

    knots_t = [zeros(r_t+1,1); repelem(linspace(0,T,N_t+1),p_t-r_t)'; ones(r_t+1,1)*T]';
    knots_x = [a*ones(r_x+1,1); repelem(linspace(a,b,N_x+1),p_x-r_x)'; ones(r_x+1,1)*b]';

    t = ((0:N_t)/N_t)*T;
    x = a+((0:N_x)/N_x)*(b-a);
    h_t = t(2);
    h_x(cont) = x(2)-x(1);


    B_t = mat_splines_exp(N_t,p_t,r_t,1,1,knots_t,T);
    C_t = mat_splines_exp(N_t,p_t,r_t,0,1,knots_t,T);
    M_t = mat_splines_exp(N_t,p_t,r_t,0,0,knots_t,T);
    M_x = mat_splines(N_x,p_x,r_x,0,0,knots_x,a,b);
    B_x_c = mat_splines_with_c(N_x,p_x,r_x,1,1,knots_x,a,b,c);

    time_matrices = toc

    tic
    siz_t = N_t*(p_t-r_t)+r_t+1;
    siz_x = N_x*(p_x-r_x)+r_x+1;
    F_vec_1 = zeros(siz_x*siz_t,1);
    F_vec_2 = zeros(siz_x*siz_t,1);
    for j_t = 1 : siz_t
        for j_x = 1 : siz_x
            for k_t = max(1,floor(j_t/(p_t-r_t))-(p_t+j_t)):min(j_t,N_t)
                phi_t = zeros(Nq,1);
                for k_x = max(1,floor(j_x/(p_x-r_x))-(p_x+r_x)):min(j_x,N_x)
                    phi_x = zeros(Nq,1);
                    dphi_x = zeros(Nq,1);
                    [xs_t,ps_t] = lgwt(Nq,t(k_t),t(k_t+1));
                    [xs_x,ps_x] = lgwt(Nq,x(k_x),x(k_x+1));
                    for iii = 1:Nq
                        phi_t(iii) = sp_and_der(p_t,knots_t,j_t-1,xs_t(iii),1)*exp(-xs_t(iii)/T);
                        phi_x(iii) = sp_and_der(p_x,knots_x,j_x-1,xs_x(iii),0);
                        dphi_x(iii) = sp_and_der(p_x,knots_x,j_x-1,xs_x(iii),1);
                    end
                    F_vec_1((j_t-1)*siz_x+j_x) = F_vec_1((j_t-1)*siz_x+j_x) + (phi_t.*ps_t)'*F(xs_t,xs_x')*(phi_x.*ps_x);
                    F_vec_1((j_t-1)*siz_x+j_x) = F_vec_1((j_t-1)*siz_x+j_x) - (phi_t.*ps_t)'*kron(ones(size(xs_t)),(c(xs_x).*dU0(xs_x))')*(dphi_x.*ps_x);
                    F_vec_2((j_t-1)*siz_x+j_x) = F_vec_2((j_t-1)*siz_x+j_x) + (phi_t.*ps_t)'*kron(ones(size(xs_t)),(V0(xs_x))')*(phi_x.*ps_x);
                end
            end
        end
    end

    F_vec_for_N = zeros(siz_x,1);
    for j_x = 2 : siz_x-1
        for k_x = max(1,floor(j_x/(p_x-r_x))-(p_x+r_x)):min(j_x,N_x)
            phi_x = zeros(Nq,1);
            [xs_x,ps_x] = lgwt(Nq,x(k_x),x(k_x+1));
            for iii = 1:Nq
                phi_x(iii) = sp_and_der(p_x,knots_x,j_x-1,xs_x(iii),0);
            end
            F_vec_for_N(j_x) = F_vec_for_N(j_x) + sum(sin(xs_x*pi).*phi_x.*ps_x);
        end
    end

    time_F = toc

    tic

    %removing initial condition in time and boundary in space
    cont_v = [];
    for j_t = 2 : siz_t
        for j_x = 2 : siz_x-1
            cont_v = [cont_v (j_t-1)*siz_x+j_x];
        end
    end
    F_vec_1 = F_vec_1(cont_v);
    F_vec_2 = F_vec_2(cont_v);

    M1 = kron(sparse(C_t),sparse(B_x_c));
    M1 = M1(cont_v,cont_v);
    M2 = kron(sparse(B_t),sparse(M_x));
    M2 = M2(cont_v,cont_v);
    M3 = kron(sparse(B_t),sparse(M_x));
    M3 = M3(cont_v,cont_v);
    M4 = -kron(sparse(C_t),sparse(M_x));
    M4 = M4(cont_v,cont_v);
    S = [M1 M2 ; M3 M4];

    SOL = S \ [F_vec_1 ; F_vec_2];

    V_app_coeff = SOL(length(SOL)/2+1:end);
    U_app_coeff = SOL(1:length(SOL)/2);

    V_N_app_coeff = -M1\(M4*V_app_coeff);
    dtV_N_coeff = B_x_c(2:end-1,2:end-1)\F_vec_for_N(2:end-1);

    tic
    U_app_coeff_w_b = zeros(siz_x*siz_t,1);
    V_app_coeff_w_b = zeros(siz_x*siz_t,1);
    V_N_app_coeff_w_b = zeros(siz_x*siz_t,1);
    cont_v = 0;
    for j_t = 1 : siz_t
        for j_x = 1 : siz_x
            if j_t == 1
                U_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
                V_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
                V_N_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
            elseif j_x == 1
                U_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
                V_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
                V_N_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
            elseif j_x == siz_x
                U_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
                V_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
                V_N_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
            else
                cont_v = cont_v+1;
                U_app_coeff_w_b((j_t-1)*siz_x+j_x) = U_app_coeff(cont_v);
                V_app_coeff_w_b((j_t-1)*siz_x+j_x) = V_app_coeff(cont_v);
                V_N_app_coeff_w_b((j_t-1)*siz_x+j_x) = V_N_app_coeff(cont_v);
            end
        end
    end

    U_app_coeff = U_app_coeff_w_b;
    V_app_coeff = V_app_coeff_w_b;
    V_N_app_coeff = V_N_app_coeff_w_b;
    dtV_N_coeff = [0; dtV_N_coeff; 0];

    [err_U(cont),err_V(cont)] = err_L2_1D(U_app_coeff,V_app_coeff,N_x,N_t,siz_x,siz_t,p_x,p_t,Nq,x,t,knots_x,knots_t,U,dtU,U0,V0)

    %if first example
    %err_norm_call_V(cont) = err_Vh(U_app_coeff,V_app_coeff,dtV_N_coeff,V_N_app_coeff,N_x,N_t,siz_x,siz_t,p_x,p_t,Nq,x,t,knots_x,knots_t,dtU,dxU,V0,dU0,c)

    errors = toc

end

%error estimates
if N1 < N2 && N_x == N_t
    plot_error({}, err_U,T./2.^(N1+1:N2+1),1);
    plot_error({}, err_V, T./2.^(N1+1:N2+1),1);
    plot_error({'L2 error U','L2 error V','norm Vh'}, err_norm_call_V, T./2.^(N1+1:N2+1),1);
    title(['Errors for splines with px=' num2str(p_x) ' rx=' num2str(r_x) 'and pt=' num2str(p_t) ' rt=' num2str(r_t)],'FontSize',16,'Interpreter','latex')
    xlabel('h')
end
