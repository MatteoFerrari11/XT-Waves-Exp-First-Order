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

% and calculates errors in the norm L^2 of U, V, dt U, dt V, dx U, dt U

% Note: this code is far from being optimized

clc
clear
close all
format long

px = 2; % Degree of spline functions in space
pt = 2; % Degree of spline functions in time
T = 1; % Final time
L = 1; % Space interval [0,L]
rx = 1; % Regularity splines in space
rt = 1; % Regularity splines in time

N_plot = 100; % Number of points for reconstructing the solution
Nq = 32; % Number of quadrature points
N1 = 2; % Starting resolution level (time and space intervals 2^(N1+1))
N2 = 5; % Ending resolution level (time and space intervals 2^(N2+1))
%also different mesh sizes for space and time can be employed

err_U = zeros(1,N2-N1+1);
err_V = zeros(1,N2-N1+1);
err_dtU = zeros(1,N2-N1+1);
err_dtV = zeros(1,N2-N1+1);
err_dxV = zeros(1,N2-N1+1);
err_dxU = zeros(1,N2-N1+1);

syms t x
%Analytical solution
%c = x;
%U = sin(pi*x)*sin(5/4*pi*t)^2 + sin(x*pi); % U(0,t)=U(1,t)=0, U(x,0) = U0(x), dtU(x,0)=V0(x)
%F = matlabFunction(simplify(diff(U,2,t) - diff(c*diff(U,x),x)));
%dtxU = matlabFunction(diff(diff(U,t),x))
%dxU = matlabFunction(diff(U,x));
%U = matlabFunction(U);

%in good vectorial notation
U = @(t,x) sin(t*pi*5/4).^2*sin(x*pi) + ones(size(t))*sin(x*pi);
U0 = @(x)  U(0,x);
dU0 = @(x) pi*cos(x*pi);
F = @(t,x) ones(size(t))*(pi*cos(x.*pi).*-2.0)...
           -ones(size(t))*(pi^2*sin(x*pi).*(2.5e+1./8.0))...
           +pi*cos(t.*pi.*5/4).^2*cos(x*pi)...
           +ones(size(t))*(x*pi.^2.*sin(x*pi).*2.0)...
           +pi^2.*cos(t*pi*5/4).^2*(sin(x*pi).*(2.5e+1./4.0))...
           -pi^2.*cos(t*pi*5/4).^2*(sin(x*pi).*x);
V0 = @(x) 0*x;
c = @(x)  x;
dtU =  @(t,x) pi*(cos(t*pi*5/4).*sin(t*pi*5/4))*sin(x*pi)*2.5;
dttU =  @(t,x)  pi.^2.*cos(t.*pi.*(5.0./4.0)).^2*sin(x*pi).*(2.5e+1./8.0)-pi.^2.*sin(t.*pi.*(5.0./4.0)).^2*sin(x.*pi).*(2.5e+1./8.0);
dxU = @(t,x)  pi*ones(size(t))*cos(x*pi)+pi*sin(t*pi*5/4).^2*cos(x*pi);
dtxU = @(t,x) pi^2*(sin(t.*pi*(5.0./4.0)).*cos(t*pi*5/4))*cos(x*pi)*5/2;

cont = 0;

for ii = N1:N2

    cont = cont + 1;
    Nx = 2^(ii+1);
    Nt = 2^(ii+1);

    knots_t = [zeros(rt+1,1); repelem(linspace(0,T,Nt+1),pt-rt)'; ones(rt+1,1)*T]';
    knots_x = [zeros(rx+1,1); repelem(linspace(0,L,Nx+1),px-rx)'; ones(rx+1,1)*L]';

    t = ((0:Nt)/Nt)*T;
    x = ((0:Nx)/Nx)*L;
    h_t = t(2);
    h_x = x(2);

    B_t = mat_splines_exp(Nt,pt,rt,1,1,knots_t,T); 
    C_t = mat_splines_exp(Nt,pt,rt,0,1,knots_t,T); 
    M_x = mat_splines(Nx,px,rx,0,0,knots_x,0,L); 
    B_x_c = mat_splines_with_c(Nx,px,rx,1,1,knots_x,0,L,c); 

    siz_t = Nt*(pt-rt)+rt+1;
    siz_x = Nx*(px-rx)+rx+1;
    F_vec_1 = zeros(siz_x*siz_t,1);
    F_vec_2 = zeros(siz_x*siz_t,1);
    for jt = 1 : siz_t
        for jx = 1 : siz_x
            for k_t = max(1,floor(jt/(pt-rt))-(pt+jt)):min(jt,Nt)
                phi_t = zeros(Nq,1);
                for k_x = max(1,floor(jx/(px-rx))-(px+rx)):min(jx,Nx)
                    phi_x = zeros(Nq,1);
                    dphi_x = zeros(Nq,1);
                    [xs_t,ps_t] = lgwt(Nq,t(k_t),t(k_t+1));
                    [xs_x,ps_x] = lgwt(Nq,x(k_x),x(k_x+1));
                    for iii = 1:Nq
                        phi_t(iii) = sp_and_der(pt,knots_t,jt-1,xs_t(iii),1)*exp(-xs_t(iii)/T);
                        phi_x(iii) = sp_and_der(px,knots_x,jx-1,xs_x(iii),0);
                        dphi_x(iii) = sp_and_der(px,knots_x,jx-1,xs_x(iii),1);
                    end
                    F_vec_1((jt-1)*siz_x+jx) = F_vec_1((jt-1)*siz_x+jx) + (phi_t.*ps_t)'*F(xs_t,xs_x')*(phi_x.*ps_x);
                    F_vec_1((jt-1)*siz_x+jx) = F_vec_1((jt-1)*siz_x+jx) - (phi_t.*ps_t)'*kron(ones(size(xs_t)),(c(xs_x).*dU0(xs_x))')*(dphi_x.*ps_x);
                    F_vec_2((jt-1)*siz_x+jx) = F_vec_2((jt-1)*siz_x+jx) + (phi_t.*ps_t)'*kron(ones(size(xs_t)),(V0(xs_x))')*(phi_x.*ps_x);
                end
            end
        end
    end

    %removing initial condition in time and boundary in space
    cont_v = [];
    for jt = 2 : siz_t
        for jx = 2 : siz_x-1
            cont_v = [cont_v (jt-1)*siz_x+jx];
        end
    end
    F_vec_1 = F_vec_1(cont_v);
    F_vec_2 = F_vec_2(cont_v);

    B_t = B_t(2:end,2:end); %zero initial condition
    C_t = C_t(2:end,2:end); %zero initial condition
    B_x_c = B_x_c(2:end-1,2:end-1); %zero boundary conditions
    M_x = M_x(2:end-1,2:end-1); %zero boundary conditions

    S = [kron(sparse(C_t),sparse(B_x_c)) kron(sparse(B_t),sparse(M_x)); ...
         kron(sparse(B_t),sparse(M_x)) -kron(sparse(C_t),sparse(M_x))];

    SOL = S \ [F_vec_1 ; F_vec_2];

    V_app_coeff = SOL(length(SOL)/2+1:end);
    U_app_coeff = SOL(1:length(SOL)/2);

    t_plot = linspace(0.0001,T-0.0001,N_plot);
    x_plot = linspace(0.0001,L-0.0001,N_plot);

    U_ex = U(t_plot',x_plot);
    dtU_ex = dtU(t_plot',x_plot);
    dttU_ex = dttU(t_plot',x_plot);
    dtxU_ex = dtxU(t_plot',x_plot);
    dxU_ex = dxU(t_plot',x_plot);
    U_app = zeros(size(U_ex));
    V_app = zeros(size(U_ex));
    dtU_app = zeros(size(U_ex));
    dxU_app = zeros(size(U_ex));
    dxV_app = zeros(size(U_ex));
    dtV_app = zeros(size(U_ex));

    %we add zero values at zero in time and at the boundary in space
    U_app_coeff_w_b = zeros(siz_x*siz_t,1);
    V_app_coeff_w_b = zeros(siz_x*siz_t,1);
    cont_v = 0;
    for jt = 1 : siz_t
        for jx = 1 : siz_x
            if jt == 1
                U_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                V_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
            elseif jx == 1
                U_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                V_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
            elseif jx == siz_x
                U_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                V_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
            else
                cont_v = cont_v+1;
                U_app_coeff_w_b((jt-1)*siz_x+jx) = U_app_coeff(cont_v);
                V_app_coeff_w_b((jt-1)*siz_x+jx) = V_app_coeff(cont_v);
            end
        end
    end

    U_app_coeff = U_app_coeff_w_b;
    V_app_coeff = V_app_coeff_w_b;

    for i_plot_t = 1 : N_plot
        for i_plot_x = 1 : N_plot
            for ind_t = 1 : siz_t
                if t_plot(i_plot_t) >= knots_t(ind_t) && t_plot(i_plot_t) < knots_t(ind_t+pt+1)
                    for ind_x = 1 : siz_x
                        if x_plot(i_plot_x) >= knots_x(ind_x) && x_plot(i_plot_x) < knots_x(ind_x+px+1)

                            U_app(i_plot_t,i_plot_x) = U_app(i_plot_t,i_plot_x)...
                                + U_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                sp_and_der(px,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                sp_and_der(px,knots_x,ind_x-1,x_plot(i_plot_x),0);

                            V_app(i_plot_t,i_plot_x) = V_app(i_plot_t,i_plot_x)...
                                + V_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                sp_and_der(px,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                sp_and_der(px,knots_x,ind_x-1,x_plot(i_plot_x),0);

                            dxU_app(i_plot_t,i_plot_x) = dxU_app(i_plot_t,i_plot_x)...
                                + U_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                sp_and_der(px,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                sp_and_der(px,knots_x,ind_x-1,x_plot(i_plot_x),1);

                            dxV_app(i_plot_t,i_plot_x) = dxV_app(i_plot_t,i_plot_x)...
                                + V_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                sp_and_der(px,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                sp_and_der(px,knots_x,ind_x-1,x_plot(i_plot_x),1);

                            dtU_app(i_plot_t,i_plot_x) = dtU_app(i_plot_t,i_plot_x)...
                                + U_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                sp_and_der(px,knots_t,ind_t-1,t_plot(i_plot_t),1)* ...
                                sp_and_der(px,knots_x,ind_x-1,x_plot(i_plot_x),0);

                            dtV_app(i_plot_t,i_plot_x) = dtV_app(i_plot_t,i_plot_x)...
                                + V_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                sp_and_der(px,knots_t,ind_t-1,t_plot(i_plot_t),1)* ...
                                sp_and_der(px,knots_x,ind_x-1,x_plot(i_plot_x),0);

                        end
                    end
                end
            end
        end
    end
    U_app = U_app + kron(ones(size(t_plot))',U0(x_plot));
    V_app = V_app + kron(ones(size(t_plot))',V0(x_plot));
    dxU_app = dxU_app + kron(ones(size(t_plot))',dU0(x_plot));

    % Errori in time
    err_U(cont) = sqrt(sum(sum(abs(U_app-U_ex).^2))/(N_plot^2))/sqrt(sum(sum(abs(U_ex).^2))/(N_plot^2))

    err_V(cont) = sqrt(sum(sum(abs(V_app-dtU_ex).^2))/(N_plot^2))/sqrt(sum(sum(abs(dtU_ex).^2))/(N_plot^2))

    err_dtV(cont) = sqrt(sum(sum(abs(dtV_app-dttU_ex).^2))/(N_plot^2))/sqrt(sum(sum(abs(dttU_ex).^2))/(N_plot^2))

    err_dtU(cont) = sqrt(sum(sum(abs(dtU_app-dtU_ex).^2))/(N_plot^2))/sqrt(sum(sum(abs(dtU_ex).^2))/(N_plot^2))

    err_dxV(cont) = sqrt(sum(sum(c(x_plot).*abs(dxV_app-dtxU_ex).^2))/(N_plot^2))/sqrt(sum(sum(c(x_plot).*abs(dtxU_ex).^2))/(N_plot^2))

    err_dxU(cont) = sqrt(sum(sum(c(x_plot).*abs(dxU_app-dxU_ex).^2))/(N_plot^2))/sqrt(sum(sum(c(x_plot).*abs(dxU_ex).^2))/(N_plot^2))
    
end

%error estimates
if N1 < N2 && Nx == Nt
    plot_error({}, err_U,T./2.^(N1+1:N2+1),1);
    plot_error({}, err_dxU,T./2.^(N1+1:N2+1),1);
    plot_error({'L2 error U','L2 error dx U', 'L2 error dt U'}, err_dtU, T./2.^(N1+1:N2+1),1);
    title(['Various errors in U for splines with px=' num2str(px) 'and pt=' num2str(pt)],'FontSize',16,'Interpreter','latex')
    xlabel('h')
    plot_error({}, err_V,T./2.^(N1+1:N2+1),2);
    plot_error({}, err_dxV,T./2.^(N1+1:N2+1),2);
    plot_error({'L2 error V','L2 error dx V', 'L2 error dt V'}, err_dtV, T./2.^(N1+1:N2+1),2);
    title(['Various errors in V for splines with px=' num2str(pt) 'and pt=' num2str(pt)],'FontSize',16,'Interpreter','latex')
    xlabel('h')
end
