function phi = phi_FEM(x,p)

switch p
    case(1)

        phi = [1-x(1)-x(2), x(1), x(2)]';
        
    case(2)
        
        phi = [(1-x(1)-x(2))*(1-2*x(1)-2*x(2)), x(1)*(2*x(1)-1), x(2)*(2*x(2)-1), 4*x(1)*(1-x(1)-x(2)), 4*x(1)*x(2), 4*x(2)*(1-x(1)-x(2))]';
        bphi = phi;
        phi(4) = bphi(5);
        phi(6) = bphi(4);
        phi(5) = bphi(6);

    case(3)
        L1 = 1 - x(1) - x(2);
        L2 = x(1);
        L3 = x(2);

        phi = [
            % --- Vertici ---
            0.5 * L1 * (3*L1 - 1) * (3*L1 - 2);      % Nodo 1
            0.5 * L2 * (3*L2 - 1) * (3*L2 - 2);      % Nodo 2
            0.5 * L3 * (3*L3 - 1) * (3*L3 - 2);      % Nodo 3

            % --- Lato orizzontale (y=0) ---
            4.5 * L1 * L2 * (3*L1 - 1);              % Nodo 4 (vicino a 1)
            4.5 * L1 * L2 * (3*L2 - 1);              % Nodo 5 (vicino a 2)

            % --- Ipotenusa ---
            4.5 * L2 * L3 * (3*L2 - 1);              % Nodo 6 (vicino a 2)
            4.5 * L2 * L3 * (3*L3 - 1);              % Nodo 7 (vicino a 3)

            % --- Lato verticale (x=0) ---
            4.5 * L3 * L1 * (3*L3 - 1);              % Nodo 8 (vicino a 3)
            4.5 * L3 * L1 * (3*L1 - 1);              % Nodo 9 (vicino a 1)

            % --- Bubble interno ---
            27 * L1 * L2 * L3                        % Nodo 10 (baricentro)
            ]';
        bphi = phi;
        phi(4:5) = bphi(6:7);
        phi(8:9) = bphi(4:5);
        phi(6:7) = bphi(8:9);
end