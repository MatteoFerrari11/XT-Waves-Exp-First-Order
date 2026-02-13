function [U,V] = solver(B_t,C_t,M_x,K_x,F_vec_1)


[Q_t,R_t] = schur(full(B_t\C_t),'complex');
Q_t = sparse(Q_t);
R_t = sparse(R_t);

eye_x = sparse(eye(size(M_x)));
Y = -kron(B_t,eye_x)\sparse(F_vec_1);
Y = kron(Q_t',eye_x)*Y;

Z = (kron(R_t,K_x)+kron(sparse(inv(R_t)),M_x))\Y;
U = -(kron(Q_t,eye_x))*Z;
V = kron(C_t,eye_x)\(kron(B_t,eye_x)*U);

U = full(U);
V = full(V);