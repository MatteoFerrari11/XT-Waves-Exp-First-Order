% Accompanying code for the Paper:
% "Inf-sup stable space--time discretization of the wave equation based on a first-order-in-time variational formulation"
% Authors: M. Ferrari and I. Perugia and E. Zampa
%
% This code generates Figure 1 with max regularity splines

clc
clear
close all
format long

p1 = 1;
p2 = 5;

T = 1; % Final time
L = 1; % Space interval [0,L]

N_plot = 100; % Number of points for reconstructing the solution
Nq = 16; % Number of quadrature points
N1 = 2; % Starting resolution level (time and space intervals 2^(N1+1))
N2 = 4; % Ending resolution level (time and space intervals 2^(N2+1))

err_U = zeros(p2-p1+1,N2-N1+1);
err_V = zeros(p2-p1+1,N2-N1+1);
err_norm_call_V = zeros(p2-p1+1,N2-N1+1);

U = @(t,x) sin(t*pi*5/4).^2*sin(x*pi) + ones(size(t))*sin(x*pi);
U0 = @(x) U(0,x);
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
dV0 = @(x) 0*x;
c = @(x) (x+1);
dtU = @(t,x) pi*(cos(t*pi*5/4).*sin(t*pi*5/4))*sin(x*pi)*2.5;
dttU = @(t,x) pi.^2.*cos(t.*pi.*(5.0./4.0)).^2*sin(x*pi).*(2.5e+1./8.0)-pi.^2.*sin(t.*pi.*(5.0./4.0)).^2*sin(x.*pi).*(2.5e+1./8.0);
dxU = @(t,x) pi*ones(size(t))*cos(x*pi)+pi*sin(t*pi*5/4).^2*cos(x*pi);
dtxU = @(t,x) pi^2*(sin(t*pi*5/4).*cos(t*pi*5/4))*cos(x*pi)*5/2;


cont_p = 0;

for p = p1:p2

    r = p-1;

    cont_p = cont_p + 1;

    cont_x = 0;

    for ii = N1:N2

        cont_x = cont_x + 1;
        Nx = 2^(ii+1);
        Nt = 2^(ii+1);

        knots_t = [zeros(r+1,1); repelem(linspace(0,T,Nt+1),p-r)'; ones(r+1,1)*T]';
        knots_x = [zeros(r+1,1); repelem(linspace(0,L,Nx+1),p-r)'; ones(r+1,1)*L]';

        t = ((0:Nt)/Nt)*T;
        x = ((0:Nx)/Nx)*L;
        h_t = t(2);
        h_x = x(2);

        B_t = mat_splines_exp(Nt,p,r,1,1,knots_t,T);
        C_t = mat_splines_exp(Nt,p,r,0,1,knots_t,T);
        M_t = mat_splines_exp(Nt,p,r,0,0,knots_t,T);
        M_x = mat_splines(Nx,p,r,0,0,knots_x,0,L);
        B_x_c = mat_splines_with_c(Nx,p,r,1,1,knots_x,0,L,c);

        siz_t = Nt*(p-r)+r+1;
        siz_x = Nx*(p-r)+r+1;
        F_vec_1 = zeros(siz_x*siz_t,1);
        F_vec_2 = zeros(siz_x*siz_t,1);
        for jt = 1 : siz_t
            for jx = 1 : siz_x
                for k_t = max(1,floor(jt/(p-r))-(p+jt)):min(jt,Nt)
                    phi_t = zeros(Nq,1);
                    for k_x = max(1,floor(jx/(p-r))-(p+r)):min(jx,Nx)
                        phi_x = zeros(Nq,1);
                        dphi_x = zeros(Nq,1);
                        [xs_t,ps_t] = lgwt(Nq,t(k_t),t(k_t+1));
                        [xs_x,ps_x] = lgwt(Nq,x(k_x),x(k_x+1));
                        for iii = 1:Nq
                            phi_t(iii) = sp_and_der(p,knots_t,jt-1,xs_t(iii),1)*exp(-xs_t(iii)/T);
                            phi_x(iii) = sp_and_der(p,knots_x,jx-1,xs_x(iii),0);
                            dphi_x(iii) = sp_and_der(p,knots_x,jx-1,xs_x(iii),1);
                        end
                        F_vec_1((jt-1)*siz_x+jx) = F_vec_1((jt-1)*siz_x+jx) + (phi_t.*ps_t)'*F(xs_t,xs_x')*(phi_x.*ps_x);
                        F_vec_1((jt-1)*siz_x+jx) = F_vec_1((jt-1)*siz_x+jx) - (phi_t.*ps_t)'*kron(ones(size(xs_t)),(c(xs_x).*dU0(xs_x))')*(dphi_x.*ps_x);
                        F_vec_2((jt-1)*siz_x+jx) = F_vec_2((jt-1)*siz_x+jx) + (phi_t.*ps_t)'*kron(ones(size(xs_t)),(V0(xs_x))')*(phi_x.*ps_x);
                    end
                end
            end
        end

        %special for this example
        F_vec_for_N = zeros(siz_x,1);
        for jx = 2 : siz_x-1
            for k_x = max(1,floor(jx/(p-r))-(p+r)):min(jx,Nx)
                phi_x = zeros(Nq,1);
                [xs_x,ps_x] = lgwt(Nq,x(k_x),x(k_x+1));
                for iii = 1:Nq
                    phi_x(iii) = sp_and_der(p,knots_x,jx-1,xs_x(iii),0);
                end
                F_vec_for_N(jx) = F_vec_for_N(jx) + sum(sin(xs_x*pi).*phi_x.*ps_x);
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

        dtV_N_coeff = [0; B_x_c\F_vec_for_N(2:end-1); 0];

        x_plot = linspace(0.0001,L-0.0001,N_plot);
        dtV_N_ex_x = zeros(length(x_plot),1);

        for i_plot_x = 1 : N_plot
            for ind_x = 1 : siz_x
                if x_plot(i_plot_x) >= knots_x(ind_x) && x_plot(i_plot_x) < knots_x(ind_x+p+1)

                    dtV_N_ex_x(i_plot_x) = dtV_N_ex_x(i_plot_x)+ dtV_N_coeff(ind_x)*sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),1);

                end
            end
        end

        V_app_coeff = SOL(length(SOL)/2+1:end);
        U_app_coeff = SOL(1:length(SOL)/2);
        V_N_app_coeff = kron(sparse(eye(size(C_t))),sparse(B_x_c))\(kron(sparse(eye(size(C_t))),sparse(M_x))*V_app_coeff);

        t_plot = linspace(0.0001,T-0.0001,N_plot);

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
        dtV_N_app = zeros(size(U_ex));

        %we add zero values at zero in time and at the boundary in space
        U_app_coeff_w_b = zeros(siz_x*siz_t,1);
        V_app_coeff_w_b = zeros(siz_x*siz_t,1);
        V_N_app_coeff_w_b = zeros(siz_x*siz_t,1);
        cont_v = 0;
        for jt = 1 : siz_t
            for jx = 1 : siz_x
                if jt == 1
                    U_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                    V_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                    V_N_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                elseif jx == 1
                    U_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                    V_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                    V_N_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                elseif jx == siz_x
                    U_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                    V_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                    V_N_app_coeff_w_b((jt-1)*siz_x+jx) = 0;
                else
                    cont_v = cont_v+1;
                    U_app_coeff_w_b((jt-1)*siz_x+jx) = U_app_coeff(cont_v);
                    V_app_coeff_w_b((jt-1)*siz_x+jx) = V_app_coeff(cont_v);
                    V_N_app_coeff_w_b((jt-1)*siz_x+jx) = V_N_app_coeff(cont_v);
                end
            end
        end

        U_app_coeff = U_app_coeff_w_b;
        V_app_coeff = V_app_coeff_w_b;
        V_N_app_coeff = V_N_app_coeff_w_b;

        for i_plot_t = 1 : N_plot
            for i_plot_x = 1 : N_plot
                for ind_t = 1 : siz_t
                    if t_plot(i_plot_t) >= knots_t(ind_t) && t_plot(i_plot_t) < knots_t(ind_t+p+1)
                        for ind_x = 1 : siz_x
                            if x_plot(i_plot_x) >= knots_x(ind_x) && x_plot(i_plot_x) < knots_x(ind_x+p+1)

                                U_app(i_plot_t,i_plot_x) = U_app(i_plot_t,i_plot_x)...
                                    + U_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);

                                V_app(i_plot_t,i_plot_x) = V_app(i_plot_t,i_plot_x)...
                                    + V_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);

                                dxU_app(i_plot_t,i_plot_x) = dxU_app(i_plot_t,i_plot_x)...
                                    + U_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),1);

                                dtU_app(i_plot_t,i_plot_x) = dtU_app(i_plot_t,i_plot_x)...
                                    + U_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),1)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),0);

                                dtV_N_app(i_plot_t,i_plot_x) = dtV_N_app(i_plot_t,i_plot_x)...
                                    + V_N_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                                    sp_and_der(p,knots_t,ind_t-1,t_plot(i_plot_t),1)* ...
                                    sp_and_der(p,knots_x,ind_x-1,x_plot(i_plot_x),1);

                            end
                        end
                    end
                end
            end
        end

        U_app = U_app + kron(ones(size(t_plot))',U0(x_plot));
        V_app = V_app + kron(ones(size(t_plot))',V0(x_plot));
        dxU_app = dxU_app + kron(ones(size(t_plot))',dU0(x_plot));
        dxV_app = dxV_app + kron(ones(size(t_plot))',dV0(x_plot));
        dtV_N_ex = kron(((25*pi^2*cos((5*pi*t_plot)/2))/8)',dtV_N_ex_x');

        % Errori in time
        err_U(cont_p,cont_x) = sqrt(sum(sum(abs(U_app-U_ex).^2))/(N_plot^2)/(sum(sum(abs(U_ex).^2))/(N_plot^2)))
        err_V(cont_p,cont_x) = sqrt(sum(sum(abs(V_app-dtU_ex).^2))/(N_plot^2)/(sum(sum(abs(dtU_ex).^2))/(N_plot^2)))

        err_V_loc = sum(sum(abs(V_app-dtU_ex).^2))/N_plot^2;
        rel_V = sum(sum(abs(dtU_ex).^2))/N_plot^2;
        err_dtVN = sum(sum(abs(dtV_N_app-dtV_N_ex).^2))/N_plot^2;
        rel_dtVN = sum(sum(abs(dtV_N_ex).^2))/N_plot^2;
        err_dtU = sum(sum(abs(dtU_app-dtU_ex).^2))/N_plot^2;
        rel_dtU = sum(sum(abs(dtU_ex).^2))/N_plot^2;
        err_dxU = sum(sum(c(x_plot).*abs(dxU_app-dxU_ex).^2))/N_plot^2;
        rel_dxU = sum(sum(c(x_plot).*abs(dxU_ex).^2))/N_plot^2;

        err_norm_call_V(cont_p,cont_x) = sqrt((err_V_loc+err_dtVN+err_dtU+err_dxU)./(rel_V+rel_dtVN+rel_dtU+rel_dxU))

    end

end

%%
close all


plot_error({}, err_norm_call_V(1,:),T./2.^(N1+1:N2+1),1);
plot_error({}, err_norm_call_V(2,:),T./2.^(N1+1:N2+1),1);
plot_error({}, err_norm_call_V(3,:),T./2.^(N1+1:N2+1),1);
plot_error({}, err_norm_call_V(4,:),T./2.^(N1+1:N2+1),1);
plot_error({'p=1','p=2','p=3','p=4','p=5'}, err_norm_call_V(5,:),T./2.^(N1+1:N2+1),1);
legend({'$p=1$','$p=2$','$p=3$','$p=4$','$p=5$'}, 'Interpreter', 'latex','FontSize', 16.7);
set(gca, 'XTick', [0.05, 0.07, 0.09,0.11]); 
set(gca,'FontSize', 14);
xlabel('$h_t=h_x$', 'Interpreter', 'latex','FontSize', 18);
set(gca, 'YTick', [0.000001,0.0001,0.01,1]); 

plot_error({}, err_U(1,:),T./2.^(N1+1:N2+1),2);
plot_error({}, err_U(2,:),T./2.^(N1+1:N2+1),2);
plot_error({}, err_U(3,:),T./2.^(N1+1:N2+1),2);
plot_error({}, err_U(4,:),T./2.^(N1+1:N2+1),2);
plot_error({'p=1','p=2','p=3','p=4','p=5'}, err_U(5,:),T./2.^(N1+1:N2+1),2);
legend({'$p=1$','$p=2$','$p=3$','$p=4$','$p=5$'}, 'Interpreter', 'latex','FontSize', 16.7);
set(gca, 'XTick', [0.05, 0.07, 0.09,0.11]); 
set(gca,'FontSize', 14);
xlabel('$h_t=h_x$', 'Interpreter', 'latex','FontSize', 18);

plot_error({}, err_V(1,:),T./2.^(N1+1:N2+1),3);
plot_error({}, err_V(2,:),T./2.^(N1+1:N2+1),3);
plot_error({}, err_V(3,:),T./2.^(N1+1:N2+1),3);
plot_error({}, err_V(4,:),T./2.^(N1+1:N2+1),3);
plot_error({'p=1','p=2','p=3','p=4','p=5'}, err_V(5,:),T./2.^(N1+1:N2+1),3);
legend({'$p=1$','$p=2$','$p=3$','$p=4$','$p=5$'}, 'Interpreter', 'latex','FontSize', 16.7);
set(gca, 'XTick', [0.05, 0.07, 0.09,0.11]); 
set(gca,'FontSize', 14);
xlabel('$h_t=h_x$', 'Interpreter', 'latex','FontSize', 18);