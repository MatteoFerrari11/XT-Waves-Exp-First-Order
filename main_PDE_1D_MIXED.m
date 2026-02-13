clc
clear
close all
format long

p_x = 2;
p_t = 2;
T = 1;
a = 0;
b = 1/2;
r_x = 1;
r_t = 1;
version = 'v2';

N_plot = 100;
Nq = 16;
N1 = 2;
N2 = 3;

h_x = zeros(1,N2-N1+1);
err_U = zeros(1,N2-N1+1);
err_V = zeros(1,N2-N1+1);

U = @(t,x) exp(-t*pi/2)*sin(x*pi/2);
dx_U = @(t,x) pi/2*exp(-t*pi/2)*cos(x*pi/2);
dt_U = @(t,x) -pi/2*exp(-t*pi/2)*sin(x*pi/2);
U0 = @(x)  sin(x*pi/2);
dU0 = @(x) pi/2*cos(x*pi/2);
F = @(t,x) (pi^2*exp(-(pi*t)/2)*sin((pi*x)/2))/2;
V0 = @(x) -pi/2*sin(x*pi/2);
c = @(x) 1;

cont = 0;

for ii = N1:N2

    tic

    cont = cont + 1;
    N_x = 2^(ii+1);
    N_t = 2^(ii+1);

    knots_t = [zeros(r_t+1,1); repelem(linspace(0,T,N_t+1),p_t-r_t)'; ones(r_t+1,1)*T]';
    knots_x = [a*ones(r_x+1,1); repelem(linspace(a,b,N_x+1),p_x-r_x)'; ones(r_x+1,1)*b]';

    t = ((0:N_t)/N_t)*T;
    x =  a+((0:N_x)/N_x)*(b-a);
    h_t = t(2);
    h_x(cont) = x(2)-x(1);

    B_t = mat_splines_exp(N_t,p_t,r_t,1,1,knots_t,T);
    C_t = mat_splines_exp(N_t,p_t,r_t,0,1,knots_t,T);
    M_x = mat_splines(N_x,p_x,r_x,0,0,knots_x,a,b);
    R_x = mat_splines_ROBIN(N_x,p_x,r_x,0,0,knots_x,a,b);
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
                tt = 1;
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
                    if tt == 1 & strcmp(version,'v2')
                        phi_x_end = sp_and_der(p_x,knots_x,j_x-1,b-eps,0);
                        F_vec_1((j_t-1)*siz_x+j_x) = F_vec_1((j_t-1)*siz_x+j_x) - (phi_t.*ps_t)'*kron(ones(size(xs_t)),(V0(b-eps))')*phi_x_end;
                        tt = 0;
                    end
                end
            end
        end
    end

    time_F = toc

    tic
    cont_v = [];
    for j_t = 2 : siz_t
        for j_x = 2 : siz_x
            cont_v = [cont_v (j_t-1)*siz_x+j_x];
        end
    end
    F_vec_1 = F_vec_1(cont_v);
    F_vec_2 = F_vec_2(cont_v);

    if strcmp(version,'v1')
        M1 = kron(sparse(C_t),sparse(B_x_c)) + kron(sparse(B_t),sparse(R_x)); %V1
        M2 = kron(sparse(B_t),sparse(M_x)); %V1
    else
        M1 = kron(sparse(C_t),sparse(B_x_c)); %V2
        M2 = kron(sparse(B_t),sparse(M_x)) + kron(sparse(C_t),sparse(R_x)); %V2
    end
    M1 = M1(cont_v,cont_v);
    M2 = M2(cont_v,cont_v); 
    M3 = kron(sparse(B_t),sparse(M_x));
    M3 = M3(cont_v,cont_v);
    M4 = -kron(sparse(C_t),sparse(M_x));
    M4 = M4(cont_v,cont_v);
    S = [M1 M2; M3 M4];

    SOL = S \ [F_vec_1 ; F_vec_2];

    V_app_coeff = SOL(length(SOL)/2+1:end);
    U_app_coeff = SOL(1:length(SOL)/2);

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
            elseif j_x == 1
                U_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
                V_app_coeff_w_b((j_t-1)*siz_x+j_x) = 0;
            else
                cont_v = cont_v+1;
                U_app_coeff_w_b((j_t-1)*siz_x+j_x) = U_app_coeff(cont_v);
                V_app_coeff_w_b((j_t-1)*siz_x+j_x) = V_app_coeff(cont_v);
            end
        end
    end

    U_app_coeff = U_app_coeff_w_b;
    V_app_coeff = V_app_coeff_w_b;

    [err_U(cont),err_V(cont)] = err_L2_1D(U_app_coeff,V_app_coeff,N_x,N_t,siz_x,siz_t,p_x,p_t,Nq,x,t,knots_x,knots_t,U,dt_U,U0,V0)

    errors = toc

end

%error estimates
if N1 < N2 && N_x == N_t
    plot_error({}, err_U,T./2.^(N1+1:N2+1),1);
    plot_error({'L2 error U','L2 error V'}, err_V, T./2.^(N1+1:N2+1),1);
    title(['Errors for splines with px=' num2str(p_x) ' rx=' num2str(r_x) 'and pt=' num2str(p_t) ' rt=' num2str(r_t)],'FontSize',16,'Interpreter','latex')
    xlabel('h')
end
