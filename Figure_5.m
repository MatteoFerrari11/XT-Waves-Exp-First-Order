% Accompanying code for the Paper:
% "Inf-sup stable space--time discretization of the wave equation based on a first-order-in-time variational formulation"
% Authors: M. Ferrari and I. Perugia and E. Zampa
%
% This code generates Figure 5

clc
clear
close all
format long

cont_p = 0;

err_U = zeros(2,6);
err_V = zeros(2,6);

for px = 2:3 % Degree of spline functions in space
    cont_p = cont_p+1;

    pt = px; % Degree of spline functions in time
    T = 1; % Final time
    L = 1.5; % Space interval [-L,L]
    rx = px-1; % Regularity splines in space
    rt = px-1; % Regularity splines in time

    N_plot = 100; % Number of points for reconstructing the solution
    Nq = 32; % Number of quadrature points
    N1 = 2; % Starting resolution level (time and space intervals 2^(N1+1))
    N2 = 4; % Ending resolution level (time and space intervals 2^(N2+1))

    w = @(x) exp(-20*(x-0.1).^2) - exp(-20*(x+0.1).^2);
    dw = @(x) exp(-20*(x+0.1).^2).*(40*x+4) - exp(-20*(x-0.1).^2).*(40*x-4);
    U = @(t,x) w(ones(size(t))*x-t*ones(size(x))+1).*(ones(size(t))*x-t*ones(size(x))+1>0);
    U0 = @(x) w(x+1).*(x+1>0);
    F = @(t,x) 0*t*x;
    V0 = @(x) -dw(x+1).*(x+1>0);
    dU0 = @(x) dw(x+1).*(x+1>0);
    c = @(x)  1+0*x;
    dtU =  @(t,x) -dw(ones(size(t))*x-t*ones(size(x))+1).*(ones(size(t))*x-t*ones(size(x))+1>0);

    cont = 0;

    for ii = N1:N2

        cont = cont + 1;
        Nx = 3*2^(ii);
        Nt = 2^(ii);

        knots_t = [zeros(rt+1,1); repelem(linspace(0,T,Nt+1),pt-rt)'; ones(rt+1,1)*T]';
        knots_x = [-ones(rx+1,1)*L; repelem(linspace(-L,L,Nx+1),px-rx)'; ones(rx+1,1)*L]';

        t = ((0:Nt)/Nt)*T;
        x = -L+((0:Nx)/Nx)*2*L;
        h_t = t(2);
        h_x = x(2)-x(1);

        B_t = mat_splines_exp(Nt,pt,rt,1,1,knots_t,T);
        C_t = mat_splines_exp(Nt,pt,rt,0,1,knots_t,T);
        M_x = mat_splines(Nx,px,rx,0,0,knots_x,-L,L);
        B_x_c = mat_splines_with_c(Nx,px,rx,1,1,knots_x,-L,L,c);

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
        x_plot = linspace(-L+0.0001,L-0.0001,N_plot);

        U_ex = U(t_plot',x_plot);
        dtU_ex = dtU(t_plot',x_plot);
        U_app = zeros(size(U_ex));
        V_app = zeros(size(U_ex));

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

                            end
                        end
                    end
                end
            end
        end
        U_app = U_app + kron(ones(size(t_plot))',U0(x_plot));
        V_app = V_app + kron(ones(size(t_plot))',V0(x_plot));

        % Errori in time
        err_U(cont_p,cont) = sqrt(sum(sum(abs(U_app-U_ex).^2))/N_plot^2)/sqrt(sum(sum(abs(U_ex).^2))/N_plot^2)
        err_V(cont_p,cont) = sqrt(sum(sum(abs(V_app-dtU_ex).^2))/N_plot^2)/sqrt(sum(sum(abs(dtU_ex).^2))/N_plot^2)

    end

end

%%%
imagesc(x_plot, t_plot, U_app);
axis xy;
colorbar;
colormap(jet);
xlabel('x');
ylabel('t');
clim([-1 1]);
%%%

%%%
plot_error({}, err_U(1,:),T./2.^(N1+1:N2+1),1);
plot_error({'p=2','p=3'}, err_U(2,:),T./2.^(N1+1:N2+1),1);
legend({'$p=2$','$p=3$'}, 'Interpreter', 'latex','FontSize', 16.7);
set(gca,'FontSize', 14);
xlabel('$h_t/h_x$', 'Interpreter', 'latex','FontSize', 18);
%%%

%%%
plot_error({}, err_V(1,:),T./2.^(N1+1:N2+1),2);
plot_error({'p=2','p=3'}, err_V(2,:),T./2.^(N1+1:N2+1),2);
legend({'$p=2$','$p=3$'}, 'Interpreter', 'latex','FontSize', 16.7);
set(gca,'FontSize', 14);
xlabel('$h_t/h_x$', 'Interpreter', 'latex','FontSize', 18);
%%%

