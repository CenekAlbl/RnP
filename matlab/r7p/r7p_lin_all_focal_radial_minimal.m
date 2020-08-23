function [vr,wr,tr,Cr,fr,rdr,err] = r7p_lin_all_focal_radial_minimal(X,u,r0,maxiter,v_init,C_init)
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
% rdr - radial distortion coefficient
% err - algebraic error at each iteration

% 2020, Zuzana Kukelova, kukelova@gmail.com


if nargin < 5
    v(:,1) = zeros(3,1);
    w(:,1) = zeros(3,1);
    t(:,1) = zeros(3,1);
    C(:,1) = zeros(3,1);
    f(:,1) = 1;
    rd(:,1) = 0;
else
    v(:,1) = v_init;
    w(:,1) = zeros(3,1);
    t(:,1) = zeros(3,1);
    C(:,1) = C_init;
    f(:,1) = 1;
    rd(:,1) = 0;
end

sample = 7;
eps = 10e-4;
not_found = 1;
err = 0;
k = 1;

while not_found && k<=maxiter
   
    wa = [];
    va = [];
    ta = [];
    Ca = [];
    fa = [];
    rda = [];
    
    try
        [wa,ta va,Ca,fa, rda] = lin_w_t_v_C_focal_radial_minimal(X',u',v(:,k),r0,sample);
    catch
        vr = v(:,k);
        wr = w(:,k);
        tr = t(:,k);
        Cr = C(:,k);
        fr = 1/f(:,k);
        rdr = rd(:,k);
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
        rdr = rd(:,k);
        return;
    end
    % from all solutions return the solution which gives the smallest error
    % on original equations
    err(:,k+1) = 100000;
    for kk = 1:size(wa,2)
        erra(kk) = test_R6P_focal_2lin_equations_focal_radial(va(:,kk),wa(:,kk),ta(:,kk),Ca(:,kk),fa(kk),rda(kk),r0, X, u);
        if erra(kk)< err(:,k+1)
            err(:,k+1) = erra(kk);
            w(:,k+1) = wa(:,kk);
            v(:,k+1) = va(:,kk);
            t(:,k+1) = ta(:,kk);
            C(:,k+1) = Ca(:,kk);
            f(:,k+1) = fa(kk);
            rd(:,k+1) = rda(kk);
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
        rdr = rd(:,k);
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
        rdr = rd(:,k+1);
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
rdr = rd(:,k);
end


%--------------------------------------------------
%compute error on original equations (not with fixed v and k) 
function err = test_R6P_focal_2lin_equations_focal_radial(vr,wr,tr,Cr,f,rd,r0, X, u)
err= 0;
for i = 1:size(vr,2)
    for j = 1:size(X,2)
       eq(:,j,i) = skewx([u(1:2,j);1+rd(i)*(u(1,j)^2+u(2,j)^2)])*[1,0,0;0,1,0;0,0,f(i)]*((eye(3)+(u(1,j)-r0)*skewx(wr(:,i)))*(eye(3)+skewx(vr(:,i)))*X(:,j)+Cr(:,i)+(u(1,j)-r0)*tr(:,i));        
    end
    err(i) = sum(sum(abs(eq(:,:,i))));
end

end



%------------------------------------------------------------------------
% solve equations using polyeig returned ff is 1/f
% version with 6x6 matrix - first eliminate some unknowns using equation
% without the focal length and radial
function [w,t,v,C,ff,rd] = lin_w_t_v_C_focal_radial_minimal(X,u,vk,r0,sample)

% how the equations and matrices were created
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
    
    syms X1 X2 X3 u1 u2 vk1 vk2 vk3
    syms t3 C3 r0 k r2 w f
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
    
    w = 1+k*r2;
    eq = [0 -w u2; w 0 -u1; -u2 u1 0]*K*(((eye(3)+[0 -v3 v2; v3 0 -v1; -v2 v1 0])+(u1-r0)*[0 -w3 w2; w3 0 -w1; -w2 w1 0]+(u1-r0)*[0 -w3 w2; w3 0 -w1; -w2 w1 0]*[0 -vk3 vk2; vk3 0 -vk1; -vk2 vk1 0])*[X1;X2;X3]+[C1;C2;C3]+(u1-r0)*[t1;t2;t3]);
    [cc , mm ] = coeffs(eq(1), [a b c C3 t3 k f])
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


%equations corresponding to the first row
% monomials

%[ a*f, a*k, a, b*f, b*k, b, c*f, c*k, c, C3*f, f*t3, f, k, 1]
r2 = u(:,1).^2+u(:,2).^2;

A = [ -u(:,2).*(X(:,2).*(n(4,1).*(r0 - u(:,1)) - n(1,1) + n(5,1).*vk(3).*(r0 - u(:,1))) + X(:,1).*(n(2,1) - n(5,1).*(r0 - u(:,1)) + n(4,1).*vk(3).*(r0 - u(:,1))) - X(:,3).*(n(4,1).*vk(1).*(r0 - u(:,1)) + n(5,1).*vk(2).*(r0 - u(:,1)))), ...
      r2(:).*(X(:,1).*(n(6,1).*(r0 - u(:,1)) - n(3,1) + n(4,1).*vk(2).*(r0 - u(:,1))) - n(8,1) + n(10,1).*(r0 - u(:,1)) + X(:,3).*(n(1,1) - n(4,1).*(r0 - u(:,1)) + n(6,1).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,1).*vk(1).*(r0 - u(:,1)) + n(6,1).*vk(3).*(r0 - u(:,1)))),...
      X(:,1).*(n(6,1).*(r0 - u(:,1)) - n(3,1) + n(4,1).*vk(2).*(r0 - u(:,1))) - n(8,1) + n(10,1).*(r0 - u(:,1)) + X(:,3).*(n(1,1) - n(4,1).*(r0 - u(:,1)) + n(6,1).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,1).*vk(1).*(r0 - u(:,1)) + n(6,1).*vk(3).*(r0 - u(:,1))),...
      -u(:,2).*(X(:,2).*(n(4,2).*(r0 - u(:,1)) - n(1,2) + n(5,2).*vk(3).*(r0 - u(:,1))) + X(:,1).*(n(2,2) - n(5,2).*(r0 - u(:,1)) + n(4,2).*vk(3).*(r0 - u(:,1))) - X(:,3).*(n(4,2).*vk(1).*(r0 - u(:,1)) + n(5,2).*vk(2).*(r0 - u(:,1)))),... 
      r2(:).*(X(:,1).*(n(6,2).*(r0 - u(:,1)) - n(3,2) + n(4,2).*vk(2).*(r0 - u(:,1))) - n(8,2) + n(10,2).*(r0 - u(:,1)) + X(:,3).*(n(1,2) - n(4,2).*(r0 - u(:,1)) + n(6,2).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,2).*vk(1).*(r0 - u(:,1)) + n(6,2).*vk(3).*(r0 - u(:,1)))),...
      X(:,1).*(n(6,2).*(r0 - u(:,1)) - n(3,2) + n(4,2).*vk(2).*(r0 - u(:,1))) - n(8,2) + n(10,2).*(r0 - u(:,1)) + X(:,3).*(n(1,2) - n(4,2).*(r0 - u(:,1)) + n(6,2).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,2).*vk(1).*(r0 - u(:,1)) + n(6,2).*vk(3).*(r0 - u(:,1))),...
      -u(:,2).*(X(:,2).*(n(4,3).*(r0 - u(:,1)) - n(1,3) + n(5,3).*vk(3).*(r0 - u(:,1))) + X(:,1).*(n(2,3) - n(5,3).*(r0 - u(:,1)) + n(4,3).*vk(3).*(r0 - u(:,1))) - X(:,3).*(n(4,3).*vk(1).*(r0 - u(:,1)) + n(5,3).*vk(2).*(r0 - u(:,1)))),...
      r2(:).*(X(:,1).*(n(6,3).*(r0 - u(:,1)) - n(3,3) + n(4,3).*vk(2).*(r0 - u(:,1))) - n(8,3) + n(10,3).*(r0 - u(:,1)) + X(:,3).*(n(1,3) - n(4,3).*(r0 - u(:,1)) + n(6,3).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,3).*vk(1).*(r0 - u(:,1)) + n(6,3).*vk(3).*(r0 - u(:,1)))),...
      X(:,1).*(n(6,3).*(r0 - u(:,1)) - n(3,3) + n(4,3).*vk(2).*(r0 - u(:,1))) - n(8,3) + n(10,3).*(r0 - u(:,1)) + X(:,3).*(n(1,3) - n(4,3).*(r0 - u(:,1)) + n(6,3).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,3).*vk(1).*(r0 - u(:,1)) + n(6,3).*vk(3).*(r0 - u(:,1))),...
      u(:,2),...
      -u(:,2).*(r0 - u(:,1)),...
      -u(:,2).*(X(:,2).*(n(4,4).*(r0 - u(:,1)) - n(1,4) + n(5,4).*vk(3).*(r0 - u(:,1))) + X(:,1).*(n(2,4) - n(5,4).*(r0 - u(:,1)) + n(4,4).*vk(3).*(r0 - u(:,1))) - X(:,3).*(n(4,4).*vk(1).*(r0 - u(:,1)) + n(5,4).*vk(2).*(r0 - u(:,1)) + 1)),...
      r2(:).*(X(:,1).*(n(6,4).*(r0 - u(:,1)) - n(3,4) + n(4,4).*vk(2).*(r0 - u(:,1))) - n(8,4) + n(10,4).*(r0 - u(:,1)) + X(:,3).*(n(1,4) - n(4,4).*(r0 - u(:,1)) + n(6,4).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,4).*vk(1).*(r0 - u(:,1)) + n(6,4).*vk(3).*(r0 - u(:,1)) + 1)),...
      X(:,1).*(n(6,4).*(r0 - u(:,1)) - n(3,4) + n(4,4).*vk(2).*(r0 - u(:,1))) - n(8,4) + n(10,4).*(r0 - u(:,1)) + X(:,3).*(n(1,4) - n(4,4).*(r0 - u(:,1)) + n(6,4).*vk(2).*(r0 - u(:,1))) - X(:,2).*(n(4,4).*vk(1).*(r0 - u(:,1)) + n(6,4).*vk(3).*(r0 - u(:,1)) + 1)];
size(A)
  % reorder
%[ a*f, a*k, a, b*f, b*k, b, c*f, c*k, c, C3*f, f*t3, f, k, 1]  
% to
  % [ a*f, a*k, a, b*f, b*k, C3*f, f*t3, b, c*f, c*k, c, f, k, 1]
A = A(:,[1,2,3,4,5, 10,11,6,7,8,9,12,13,14]); 
%Ar = gj(A); 
Ar = A(:,1:7) \ A(:,1:14);
 
data = [Ar(1,8:14),Ar(2,8:14),Ar(3,8:14),Ar(4,8:14),Ar(5,8:14)];

% solve using GB
sols = solver_R7Pfr_gb(data);

f = sols(3,:);

% select feasible solutions
indx = find(~isinf(f) & f>0 & abs(f)>10e-6 & abs(imag(f))<10e-6);

ff = f(indx);
b = sols(1,indx);
c = sols(2,indx);
rd= sols(4,indx);


v= [];
w= [];
C= [];
t= [];

% extract solutions
for i=1:length(ff)
    mon = [b(i), c(i)*ff(i), c(i)*rd(i), c(i), ff(i), rd(i), 1]';
    a(i) = -Ar(3,8:end)*mon;
    C(3,i) = (-Ar(6,8:end)*mon)/ff(i);
    t(3,i) = (-Ar(7,8:end)*mon)/ff(i);
    
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

function sols = solver_R7Pfr_gb(data)
display("size data")
size(data)
[C0,C1] = setup_elimination_template(data);
display("size C0 C1")
size(C0)
size(C1)
C1 = C0 \ C1;
RR = [-C1(end-6:end,:);eye(10)];
AM_ind = [1,2,9,3,4,11,5,6,14,7];
AM = RR(AM_ind,:);
display("size AM")
size(AM)
[V,D] = eig(AM);
scale = (diag(D).' ./ ((V(10,:))));
V = V .* (ones(size(V,1),1)*scale);
sols(1,:) = V(3,:);
sols(2,:) = V(9,:);
sols(3,:) = diag(D).';
sols(4,:) = V(5,:) ./ (sols(2,:).*sols(2,:));
end
% Action =  z
% Quotient ring basis (V) = x*y,x*z,x,y^2*z,y^2*w,y^2,y*z,y*w,y,z,
% Available monomials (RR*V) = x*y*z,x*z^2,y^2*z^2,y^2*z*w,y*z^2,y*z*w,z^2,x*y,x*z,x,y^2*z,y^2*w,y^2,y*z,y*w,y,z,

function [coeffs] = compute_coeffs(data)
coeffs(1) = -data(16);
coeffs(2) = -data(17);
coeffs(3) = -data(15);
coeffs(4) = data(2) - data(18);
coeffs(5) = -data(19);
coeffs(6) = data(3);
coeffs(7) = -data(20);
coeffs(8) = data(1);
coeffs(9) = data(4);
coeffs(10) = data(5) - data(21);
coeffs(11) = data(6);
coeffs(12) = data(7);
coeffs(13) = data(9);
coeffs(14) = data(10) - data(18);
coeffs(15) = data(8);
coeffs(16) = data(11);
coeffs(17) = data(12);
coeffs(18) = data(13) - data(21);
coeffs(19) = data(14);
coeffs(20) = 1;
coeffs(21) = data(23);
coeffs(22) = data(24);
coeffs(23) = data(22);
coeffs(24) = data(25);
coeffs(25) = data(26);
coeffs(26) = data(27);
coeffs(27) = data(28);
coeffs(28) = data(30);
coeffs(29) = data(31);
coeffs(30) = data(29);
coeffs(31) = data(32);
coeffs(32) = data(33);
coeffs(33) = data(34);
coeffs(34) = data(35);
end

function [C0,C1] = setup_elimination_template(data)
[coeffs] = compute_coeffs(data);
coeffs0_ind = [1,1,1,1,1,13,20,4,3,2,1,14,20,3,20,5,1,5,1,6,2,20,6,2,22,29,3,15,20,8,3,20,9,8,6,14,3,23,22,20,30,29,10,7,5,3,18,20,3,20,...
11,7,20,11,7,26,22,29,6,2,33,8,15,23,30,12,11,18,23,30,8,26,3,20,33,18,25,32,10,5,26,7,33,26,33,11,7,27,34,12,18,33,26,11,19,34,27,12,4,13,...
16,23,15,21,3,20,30,28,5,17,20,21,13,1,28,4,1,21,22,14,2,28,29,13,25,17,21,1,5,28,32,10,14,5,25,26,18,21,28,4,1,22,2,7,29,32,33,17,25,5,...
32];
coeffs1_ind = [9,16,15,24,8,23,30,31,10,17,15,19,25,23,20,3,30,32,12,19,27,15,30,23,8,34,13,24,16,4,21,28,31,9,14,24,6,22,29,31,16,9,24,31,16,17,27,19,13,24,...
28,21,4,10,31,25,32,34,12,18,27,24,31,9,14,29,22,6,11,26,33,34,19,16,31,24,9,12,27,34,19,17,27,32,25,10,34];
C0_ind = [1,6,28,36,55,58,61,79,80,81,82,84,86,88,103,106,109,114,117,131,134,154,158,163,164,180,185,188,196,209,212,234,235,236,237,238,241,242,248,257,258,260,261,263,264,265,266,271,273,281,...
287,290,298,314,319,320,323,324,325,327,336,341,342,352,364,365,367,368,375,376,377,378,379,381,390,395,401,402,403,405,406,409,411,427,428,429,431,453,454,455,457,459,460,461,483,485,486,487,497,498,...
500,503,504,508,514,516,519,520,523,526,536,555,556,566,571,574,579,580,581,582,592,596,597,603,607,608,614,617,618,619,623,626,629,631,632,633,634,635,636,637,639,640,643,644,645,648,649,655,666,669,...
671];
C1_ind = [3,4,7,14,20,22,23,26,29,30,31,32,40,42,44,45,47,52,55,56,66,67,69,70,71,78,85,87,88,98,100,101,103,106,111,112,124,126,127,128,137,150,152,153,161,163,165,166,171,172,...
173,174,175,176,177,178,179,181,184,189,190,193,194,195,197,199,200,201,202,204,205,206,215,223,225,226,227,228,230,231,239,249,250,251,252,253,255];
C0 = zeros(26,26);
C1 = zeros(26,10);
C0(C0_ind) = coeffs(coeffs0_ind);
C1(C1_ind) = coeffs(coeffs1_ind);
end

function mx = skewx(x)
    mx = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
end