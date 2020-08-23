function [vr,wr,tr,Cr,fr,k,err] = r7p_lin_all_focal(X,u,r0,maxiter,v_init,C_init)
% INPUT
%
% X - 3x7 matrix of 3D points
% u - 2x7 matrix of 2D correspondences
% r0 - the image coordinate around which the RS model is linearized
% default: 0 (center of the image) 
% maxiter - maximum number of iterations of the linear solver
% v_init - initial value of camera orientation (I + skew(v))
% C_init - initial value of camera center
%
% OUTPUT
%
% vr - camera orientation (I + skew(v))
% wr - camera rotational velocity (I + u(1,i)*skew(w))
% tr - camera translational velocity u(1,i)*t
% fr - focal length
% k - number of performed iterations
% err - algebraic error at each iteration

% 2020, Zuzana Kukelova, kukelova@gmail.com

if nargin < 5
    v(:,1) = zeros(3,1);
    w(:,1) = zeros(3,1);
    t(:,1) = zeros(3,1);
    C(:,1) = zeros(3,1);
    f(:,1) = 1;
else
    v(:,1) = v_init;
    w(:,1) = zeros(3,1);
    t(:,1) = zeros(3,1);
    C(:,1) = C_init;
    f(:,1) = 1;
end

%maxiter = 5;
sample = 7;
eps = 10e-4;
not_found = 1;
err = 0;
k = 1;

while not_found && k<=maxiter
    
    %[w(:,k+1),t(:,k+1) v(:,k+1),C(:,k+1),f(:,k+1)] = lin_w_t_v_C_focal(X',u',v(:,k),r0,sample);
    wa = [];
    va = [];
    ta = [];
    Ca = [];
    fa = [];
    
    try
        [wa,ta va,Ca,fa] = lin_w_t_v_C_focal_6(X',u',v(:,k),r0,sample);
    catch
        vr = v(:,k);
        wr = w(:,k);
        tr = t(:,k);
        Cr = C(:,k);
        fr = 1/f(:,k);
        return;
    end
    if isempty(wa)
        %if polyeig didn't return any feasible solution return solution
        %from previous iteration
        vr = v(:,k);
        wr = w(:,k);
        tr = t(:,k);
        Cr = C(:,k);
        fr = 1/f(:,k);
        return;
    end
    % from all solutions return the solution which gives the smallest error
    % on original equations
    err(:,k+1) = 100000;
    for kk = 1:size(wa,2)
        erra(kk) = test_R6P_focal_2lin_equations_focal(va(:,kk),wa(:,kk),ta(:,kk),Ca(:,kk),fa(kk),r0, X, u);
        if erra(kk)< err(:,k+1)
            err(:,k+1) = erra(kk);
            w(:,k+1) = wa(:,kk);
            v(:,k+1) = va(:,kk);
            t(:,k+1) = ta(:,kk);
            C(:,k+1) = Ca(:,kk);
            f(:,k+1) = fa(kk);
        end
    end
    % if solutions from polyeig failed, return solution from previous
    % iteration
    if err(:,k+1) == 100000
        vr = v(:,k);
        wr = w(:,k);
        tr = t(:,k);
        Cr = C(:,k);
        fr = 1/f(:,k);
        
        return;
    end
    
    % if error on original equations is small return the solution
    if err(:,k+1) < eps;
        not_found = 0;    
        vr = v(:,k+1);
        wr = w(:,k+1);
        tr = t(:,k+1);
        Cr = C(:,k+1);
        fr = 1/f(:,k+1);
        return;
    else
        k = k+1;
    end
end

vr = v(:,k);
wr = w(:,k);
tr = t(:,k);
Cr = C(:,k);
fr = 1/f(:,k);
end


%------------------------------------------------------------------------
% solve equations using polyeig returned ff is 1/f
function [w,t,v,C,ff] = lin_w_t_v_C_focal(X,u,vk,r0,sample)

% how the equations and matrices A0 and A1 were created
symbolic = 0;
if symbolic
    syms X1 X2 X3 u1 u2 v1 v2 v3 w1 w2 w3 vk1 vk2 vk3 t1 t2 t3 C1 C2 C3 r0 f
    K = [1 0 0; 0 1 0; 0 0 f];
    eq = [0 -1 u2; 1 0 -u1; -u2 u1 0]*K*(((eye(3)+[0 -v3 v2; v3 0 -v1; -v2 v1 0])+(u1-r0)*[0 -w3 w2; w3 0 -w1; -w2 w1 0]+(u1-r0)*[0 -w3 w2; w3 0 -w1; -w2 w1 0]*[0 -vk3 vk2; vk3 0 -vk1; -vk2 vk1 0])*[X1;X2;X3]+[C1;C2;C3]+(u1-r0)*[t1;t2;t3]);
    [c , m ] = coeffs(eq(1), [v1, v2, v3, w1, w2, w3, C1, C2, C3 t1 t2 t3])
end

%[ v1, v2, v3, w1, w2, w3, C2, C3, t2, t3, 1]
%[ X3 + X2*f*u2, -X1*f*u2, -X1, X1*vk2*(r0 - u1) - f*u2*(X2*(r0 - u1) + X1*vk3*(r0 - u1) - X3*vk1*(r0 - u1)) - X3*(r0 - u1) - X2*vk1*(r0 - u1), f*u2*(X1*(r0 - u1) - X2*vk3*(r0 - u1) + X3*vk2*(r0 - u1)), X1*(r0 - u1) - X2*vk3*(r0 - u1) + X3*vk2*(r0 - u1), -1, f*u2, r0 - u1, -f*u2*(r0 - u1), X3*f*u2 - X2]
%[ -X2*f*u1, X3 + X1*f*u1, -X2, f*u1*(X2*(r0 - u1) + X1*vk3*(r0 - u1) - X3*vk1*(r0 - u1)), X1*vk2*(r0 - u1) - f*u1*(X1*(r0 - u1) - X2*vk3*(r0 - u1) + X3*vk2*(r0 - u1)) - X3*(r0 - u1) - X2*vk1*(r0 - u1), X2*(r0 - u1) + X1*vk3*(r0 - u1) - X3*vk1*(r0 - u1), 1, -f*u1, u1 - r0, f*u1*(r0 - u1), X1 - X3*f*u1]
%[ -X3*u1, -X3*u2, X1*u1 + X2*u2, u1*(X3*(r0 - u1) - X1*vk2*(r0 - u1) + X2*vk1*(r0 - u1)), u2*(X3*(r0 - u1) - X1*vk2*(r0 - u1) + X2*vk1*(r0 - u1)), - u1*(X1*(r0 - u1) - X2*vk3*(r0 - u1) + X3*vk2*(r0 - u1)) - u2*(X2*(r0 - u1) + X1*vk3*(r0 - u1) - X3*vk1*(r0 - u1)), -u2, u1, u2*(r0 - u1), -u1*(r0 - u1), X2*u1 - X1*u2]

A0 = [ X(:,3), zeros(sample,1), -X(:,1), X(:,1).*vk(2).*(r0 - u(:,1)) - X(:,3).*(r0 - u(:,1)) - X(:,2).*vk(1).*(r0 - u(:,1)), zeros(sample,1), X(:,1).*(r0 - u(:,1)) - X(:,2).*vk(3).*(r0 - u(:,1)) + X(:,3).*vk(2).*(r0 - u(:,1)), zeros(sample,1), -ones(sample,1), zeros(sample,1), zeros(sample,1), r0 - u(:,1), zeros(sample,1), - X(:,2); ...
    zeros(sample,1), X(:,3), -X(:,2),  zeros(sample,1), X(:,1).*vk(2).*(r0 - u(:,1)) - X(:,3).*(r0 - u(:,1)) - X(:,2).*vk(1).*(r0 - u(:,1)), X(:,2).*(r0 - u(:,1)) + X(:,1).*vk(3).*(r0 - u(:,1)) - X(:,3).*vk(1).*(r0 - u(:,1)), ones(sample,1), zeros(sample,1) ,zeros(sample,1), u(:,1) - r0, zeros(sample,1), zeros(sample,1), X(:,1)];

A1 = [X(:,2).*u(:,2), -X(:,1).*u(:,2), zeros(sample,1), - u(:,2).*(X(:,2).*(r0 - u(:,1))+ X(:,1).*vk(3).*(r0 - u(:,1)) - X(:,3).*vk(1).*(r0 - u(:,1))), u(:,2).*(X(:,1).*(r0 - u(:,1))- X(:,2).*vk(3).*(r0 - u(:,1)) + X(:,3).*vk(2).*(r0 - u(:,1))), zeros(sample,1), zeros(sample,1), zeros(sample,1), u(:,2), zeros(sample,1), zeros(sample,1), -u(:,2).*(r0 - u(:,1)), X(:,3).*u(:,2); ...
    -X(:,2).*u(:,1), X(:,1).*u(:,1), zeros(sample,1), u(:,1).*(X(:,2).*(r0 - u(:,1))+ X(:,1).*vk(3).*(r0 - u(:,1)) - X(:,3).*vk(1).*(r0 - u(:,1))) , - u(:,1).*(X(:,1).*(r0 - u(:,1))- X(:,2).*vk(3).*(r0 - u(:,1)) + X(:,3).*vk(2).*(r0 - u(:,1))) , zeros(sample,1), zeros(sample,1), zeros(sample,1) ,-u(:,1), zeros(sample,1), zeros(sample,1), u(:,1).*(r0 - u(:,1)),- X(:,3).*u(:,1)];

% solve using polyeig
[n, f] = polyeig(A0(1:13,:),A1(1:13,:));

% select feasible solutions
indx = find(~isinf(f) & f>0 & abs(f)>10e-6 & abs(imag(f))<10e-6);

ff = f(indx);

v= [];
w= [];
C= [];
t= [];
j = 1;
% extract solutions
for i=indx'
    v(:,j) = n(1:3,i)./n(13,i);
    w(:,j) = n(4:6,i)./n(13,i);
    C(:,j) = n(7:9,i)./n(13,i);
    t(:,j) = n(10:12,i)./n(13,i);
    j = j+1;
end
end

%--------------------------------------------------
%compute error on original equations (not with fixed v) 
function err = test_R6P_focal_2lin_equations_focal(vr,wr,tr,Cr,f,r0, X, u)
err= 0;
for i = 1:size(vr,2)
    for j = 1:size(X,2)
       eq(:,j,i) = skewx([u(:,j);1])*[1,0,0;0,1,0;0,0,f(i)]*((eye(3)+(u(1,j)-r0)*skewx(wr(:,i)))*(eye(3)+skewx(vr(:,i)))*X(:,j)+Cr(:,i)+(u(1,j)-r0)*tr(:,i));        
    end
    err(i) = sum(sum(abs(eq(:,:,i))));
end

end



%------------------------------------------------------------------------
% solve equations using polyeig returned ff is 1/f
% version with 6x6 matrix - first eliminate some unknowns using equation
% without the focal length
function [w,t,v,C,ff,k] = lin_w_t_v_C_focal_6(X,u,vk,r0,sample)

% how the equations and matrices A0 and A1 were created
symbolic = 0;
if symbolic
    syms X1 X2 X3 u1 u2 v1 v2 v3 w1 w2 w3 vk1 vk2 vk3 t1 t2 t3 C1 C2 C3 r0 k r2 f
    K = [1 0 0; 0 1 0; 0 0 f];
    %r2 = u1^2+u2^2
    %w = 1+k*(u1^2+u2^2);
    w = 1+k*r2;
   
    % radial + focal
    eq = [0 -w u2; w 0 -u1; -u2 u1 0]*K*(((eye(3)+[0 -v3 v2; v3 0 -v1; -v2 v1 0])+(u1-r0)*[0 -w3 w2; w3 0 -w1; -w2 w1 0]+(u1-r0)*[0 -w3 w2; w3 0 -w1; -w2 w1 0]*[0 -vk3 vk2; vk3 0 -vk1; -vk2 vk1 0])*[X1;X2;X3]+[C1;C2;C3]+(u1-r0)*[t1;t2;t3]);
    [c , m ] = coeffs(eq(3), [v1, v2, v3, w1, w2, w3, C1, C2, C3 t1 t2 t3])
    
    syms t3 C3 r0 k r2 f
    K = [1 0 0; 0 1 0; 0 0 f];
    syms a b c
    syms n11 n21 n31 n41 n51 n61 n71 n81 n91 n101
    syms n12 n22 n32 n42 n52 n62 n72 n82 n92 n102
    syms n13 n23 n33 n43 n53 n63 n73 n83 n93 n103
    syms n14 n24 n34 n44 n54 n64 n74 n84 n94 n104
    v1 = a*n11+b*n12+c*n13+n14
    v2 = a*n21+b*n22+c*n23+n24
    v3 = a*n31+b*n32+c*n33+n34
    w1 = a*n41+b*n42+c*n43+n44
    w2 = a*n51+b*n52+c*n53+n54
    w3 = a*n61+b*n62+c*n63+n64
    C1 = a*n71+b*n72+c*n73+n74
    C2 = a*n81+b*n82+c*n83+n84
    t1 = a*n91+b*n92+c*n93+n94
    t2 = a*n101+b*n102+c*n103+n104
    
    eq = [0 -1 u2; 1 0 -u1; -u2 u1 0]*K*(((eye(3)+[0 -v3 v2; v3 0 -v1; -v2 v1 0])+(u1-r0)*[0 -w3 w2; w3 0 -w1; -w2 w1 0]+(u1-r0)*[0 -w3 w2; w3 0 -w1; -w2 w1 0]*[0 -vk3 vk2; vk3 0 -vk1; -vk2 vk1 0])*[X1;X2;X3]+[C1;C2;C3]+(u1-r0)*[t1;t2;t3]);
    [c , m ] = coeffs(eq(1), [a b c C3 t3])
end

%[ v1, v2, v3, w1, w2, w3, C1, C2, t1, t2, 1]
%[ -X3*u1, -X3*u2, X1*u1 + X2*u2, u1*(X3*(r0 - u1) - X1*vk2*(r0 - u1) + X2*vk1*(r0 - u1)), u2*(X3*(r0 - u1) - X1*vk2*(r0 - u1) + X2*vk1*(r0 - u1)), - u1*(X1*(r0 - u1) - X2*vk3*(r0 - u1) + X3*vk2*(r0 - u1)) - u2*(X2*(r0 - u1) + X1*vk3*(r0 - u1) - X3*vk1*(r0 - u1)), -u2, u1, u2*(r0 - u1), -u1*(r0 - u1), X2*u1 - X1*u2]
 
% equation corresponding to the 3rd row (without focal length)

% 7x11 matrix
A = [ -X(:,3).*u(:,1), -X(:,3).*u(:,2), X(:,1).*u(:,1) + X(:,2).*u(:,2), u(:,1).*(X(:,3).*(r0 - u(:,1)) - X(:,1).*vk(2).*(r0 - u(:,1)) + X(:,2).*vk(1).*(r0 - u(:,1))), u(:,2).*(X(:,3).*(r0 - u(:,1)) - X(:,1).*vk(2).*(r0 - u(:,1)) + X(:,2).*vk(1).*(r0 - u(:,1))), - u(:,1).*(X(:,1).*(r0 - u(:,1)) - X(:,2).*vk(3).*(r0 - u(:,1)) + X(:,3).*vk(2).*(r0 - u(:,1))) - u(:,2).*(X(:,2).*(r0 - u(:,1)) + X(:,1).*vk(3).*(r0 - u(:,1)) - X(:,3).*vk(1).*(r0 - u(:,1))), -u(:,2), u(:,1), u(:,2).*(r0 - u(:,1)), -u(:,1).*(r0 - u(:,1)), X(:,2).*u(:,1) - X(:,1).*u(:,2)];
 
nn = null(A);

for i = 1:size(nn,1)-1
    n(i,1) = nn(i,1)-(nn(i,4)*nn(11,1)/nn(11,4));
    n(i,2) = nn(i,2)-(nn(i,4)*nn(11,2)/nn(11,4));
    n(i,3) = nn(i,3)-(nn(i,4)*nn(11,3)/nn(11,4));
    n(i,4) = nn(i,4)/nn(11,4);
end

 
A0 = [X(:,1).*(n(6,1).*(r0 - u(:,1)) - n(3,1) + n(4,1).*vk(2).*(r0 - u(:,1))) - n(8,1) + n(10,1).*(r0 - u(:,1)) + X(:,3).*(n(1,1) - n(4,1).*(r0 - u(:,1)) + n(6,1).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,1).*vk(1).*(r0 - u(:,1)) + n(6,1).*vk(3).*(r0 - u(:,1))),...
      X(:,1).*(n(6,2).*(r0 - u(:,1)) - n(3,2) + n(4,2).*vk(2).*(r0 - u(:,1))) - n(8,2) + n(10,2).*(r0 - u(:,1)) + X(:,3).*(n(1,2) - n(4,2).*(r0 - u(:,1)) + n(6,2).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,2).*vk(1).*(r0 - u(:,1)) + n(6,2).*vk(3).*(r0 - u(:,1))),...
      X(:,1).*(n(6,3).*(r0 - u(:,1)) - n(3,3) + n(4,3).*vk(2).*(r0 - u(:,1))) - n(8,3) + n(10,3).*(r0 - u(:,1)) + X(:,3).*(n(1,3) - n(4,3).*(r0 - u(:,1)) + n(6,3).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,3).*vk(1).*(r0 - u(:,1)) + n(6,3).*vk(3).*(r0 - u(:,1))),...
      zeros(sample,1),zeros(sample,1),...
      X(:,1).*(n(6,4).*(r0 - u(:,1)) - n(3,4) + n(4,4).*vk(2).*(r0 - u(:,1))) - n(8,4) + n(10,4).*(r0 - u(:,1)) + X(:,3).*(n(1,4) - n(4,4).*(r0 - u(:,1)) + n(6,4).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,4).*vk(1).*(r0 - u(:,1)) + n(6,4).*vk(3).*(r0 - u(:,1))+1)
      ];


A1 = [-u(:,2).*(X(:,2).*(n(4,1).*(r0 - u(:,1)) - n(1,1) + n(5,1).*vk(3).*(r0 - u(:,1))) + X(:,1).*(n(2,1) - n(5,1).*(r0 - u(:,1)) + n(4,1).*vk(3)*(r0 - u(:,1))) - X(:,3).*(n(4,1).*vk(1)*(r0 - u(:,1)) + n(5,1).*vk(2).*(r0 - u(:,1)))),...
      -u(:,2).*(X(:,2).*(n(4,2).*(r0 - u(:,1)) - n(1,2) + n(5,2).*vk(3).*(r0 - u(:,1))) + X(:,1).*(n(2,2) - n(5,2).*(r0 - u(:,1)) + n(4,2).*vk(3)*(r0 - u(:,1))) - X(:,3).*(n(4,2).*vk(1)*(r0 - u(:,1)) + n(5,2).*vk(2).*(r0 - u(:,1)))),...
      -u(:,2).*(X(:,2).*(n(4,3).*(r0 - u(:,1)) - n(1,3) + n(5,3).*vk(3).*(r0 - u(:,1))) + X(:,1).*(n(2,3) - n(5,3).*(r0 - u(:,1)) + n(4,3).*vk(3)*(r0 - u(:,1))) - X(:,3).*(n(4,3).*vk(1)*(r0 - u(:,1)) + n(5,3).*vk(2).*(r0 - u(:,1)))),...
      u(:,2), -u(:,2).*(r0 - u(:,1)), ...
      -u(:,2).*(X(:,2).*(n(4,4).*(r0 - u(:,1)) - n(1,4) + n(5,4).*vk(3).*(r0 - u(:,1))) + X(:,1).*(n(2,4) - n(5,4).*(r0 - u(:,1)) + n(4,4).*vk(3)*(r0 - u(:,1))) - X(:,3).*(n(4,4).*vk(1)*(r0 - u(:,1)) + n(5,4).*vk(2).*(r0 - u(:,1))+1))];
      
 
% solve using polyeig
[xx, f] = polyeig(A0(1:6,:),A1(1:6,:));

% select feasible solutions
indx = find(~isinf(f) & f>0 & abs(f)>10e-6 & abs(imag(f))<10e-6);

ff = f(indx);

v= [];
w= [];
C= [];
t= [];
j = 1;
% extract solutions
for i=indx'
    a(j) = xx(1,i)./xx(6,i);
    b(j) = xx(2,i)./xx(6,i);
    c(j) = xx(3,i)./xx(6,i);
    C(3,j) = xx(4,i)./xx(6,i);
    t(3,j) = xx(5,i)./xx(6,i);
    j = j+1;
end

for i = 1:length(a)
    v(1,i) = a(i)*n(1,1)+b(i)*n(1,2)+c(i)*n(1,3)+n(1,4);
    v(2,i) = a(i)*n(2,1)+b(i)*n(2,2)+c(i)*n(2,3)+n(2,4);
    v(3,i) = a(i)*n(3,1)+b(i)*n(3,2)+c(i)*n(3,3)+n(3,4);
    w(1,i) = a(i)*n(4,1)+b(i)*n(4,2)+c(i)*n(4,3)+n(4,4);
    w(2,i) = a(i)*n(5,1)+b(i)*n(5,2)+c(i)*n(5,3)+n(5,4);
    w(3,i) = a(i)*n(6,1)+b(i)*n(6,2)+c(i)*n(6,3)+n(6,4);
    C(1,i) = a(i)*n(7,1)+b(i)*n(7,2)+c(i)*n(7,3)+n(7,4);
    C(2,i) = a(i)*n(8,1)+b(i)*n(8,2)+c(i)*n(8,3)+n(8,4);
    t(1,i) = a(i)*n(9,1)+b(i)*n(9,2)+c(i)*n(9,3)+n(9,4);
    t(2,i) = a(i)*n(10,1)+b(i)*n(10,2)+c(i)*n(10,3)+n(10,4);

end
end

function mx = skewx(x)
    mx = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
end