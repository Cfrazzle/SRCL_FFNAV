syms q qm q1 q2 qm1 qm2
syms q_dot qm_dot q1_dot q2_dot qm1_dot qm2_dot
syms q_ddot qm_ddot q1_ddot q2_ddot qm1_ddot qm2_ddot

q = [q1 q2].';
qm = [qm1 qm2].';
q_dot = [q1_dot q2_dot].';
qm_dot = [qm1_dot qm2_dot].';
q_ddot = [q1_ddot q2_ddot].';
qm_ddot = [qm1_ddot qm2_ddot].';

x = [ q; q_dot; qm; qm_dot];

q_ddot = @(q, qm, q_dot, qm_dot) inv(M(q))*(k(q, qm)(qm-q) - f(q_dot) - C(q,q_dot)*q_dot - S*qm_ddot);
qm_ddot = @(q, qm, q_dot, qm_dot) inv(J_m)* (tau-k(q, qm)*(qm-q) - S'*q_ddot) ;

H = [   zeros(2,2)  zeros(2,2)  eye(2,2)    zeros(2,2) 
        zeros(2,2)  zeros(2,2)  zeros(2,2)  eye(2,2) ];
    
f = [ q_dot;  q_ddot;   qm_dot;   qm_ddot];

%==========================================================================
%Constructing the Observability Matrix

%Zeroth-Order Lie Derivative Gradient
G0 = H;

%First-Order Lie Derivative Gradient
L1 = H*f;
G1 = jacobian(L1, x);

%Second-Order Lie Derivative Gradient
L2 = G1*f;
G2 = jacobian(L2, x);


%Third-Order Lie Derivative Gradient
L3 = G2*f;
G3 = jacobian(L3, x);

OM = [G0; G1; G2; G3];
OM_rank = rank(OM)

       
