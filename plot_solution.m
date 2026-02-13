N_plot = 100;
t_plot = linspace(0+0.001,T-0.0001,N_plot);
x_plot = linspace(a+0.001,b-0.0001,N_plot);

U_app = zeros(N_plot,N_plot);
V_app = zeros(N_plot,N_plot);

for i_plot_t = 1 : N_plot
    for i_plot_x = 1 : N_plot
        for ind_t = 1 : siz_t
            if t_plot(i_plot_t) >= knots_t(ind_t) && t_plot(i_plot_t) < knots_t(ind_t+p_t+1)
                for ind_x = 1 : siz_x
                    if x_plot(i_plot_x) >= knots_x(ind_x) && x_plot(i_plot_x) < knots_x(ind_x+p_x+1)

                        U_app(i_plot_t,i_plot_x) = U_app(i_plot_t,i_plot_x)...
                            + U_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                            sp_and_der(p_x,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                            sp_and_der(p_x,knots_x,ind_x-1,x_plot(i_plot_x),0);

                        V_app(i_plot_t,i_plot_x) = V_app(i_plot_t,i_plot_x)...
                            + V_app_coeff((ind_t-1)*siz_x+ind_x)* ...
                            sp_and_der(p_x,knots_t,ind_t-1,t_plot(i_plot_t),0)* ...
                            sp_and_der(p_x,knots_x,ind_x-1,x_plot(i_plot_x),0);

                    end
                end
            end
        end
    end
end
U_app = U_app + kron(ones(size(t_plot))',U0(x_plot));
V_app = V_app + kron(ones(size(t_plot))',V0(x_plot));

%close all
surf(x_plot,t_plot,U_app)

%hold on
%surf(x_plot,t_plot,U(t_plot',x_plot))
az = 170;
el = 30;
view(az, el);
shg