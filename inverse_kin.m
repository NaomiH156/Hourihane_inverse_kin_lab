function output_thetas = inverse_kin(g_des)
%Solves the inverse kinematics of the intellidex robot

%THIS ISN'T WORKING BECAUSE g1*Pwhatever is trying to multiply a 4x4 matrix
    %by a 3x1 vector.


%Segment lengths, in inches
L0 = -13;
L1 = 14.7;
L2 = 12;
L3 = 12;
L4 = 9;

%this vector keeps coming up so i'm just gonna name it
down1 = [0; 0; -1];

%defining the positions along the arm AT ZERO CONFIGURATION:
Psh = [0; -L0; L1];     %position of "shoulder"
Psh_p = Psh + down1;    %point that lies on w2 but not on w0 or w1. "Psh prime"
Pe = Psh + [L2; 0; 0];  %position of "elbow"
Pe_p = Pe + down1;      %"Pe prime," a point that lies on w3 but no other axes.
Pw = Pe + [L3; 0; 0];   %position of "wrist"
Pw_p = Pe + down1;      %"Pw prime", a point that lies only on w4.

%defining directions of axes AT ZERO CONFIGURATION:
w0 = [0 0 1].';
w1 = [0 -1 0].';
w2 = [0 0 -1].';
w3 = [0 0 -1].';
w4 = [0 0 -1].';
w5 = [1 0 0].';

%Zero-angle configuration, gst_0
gst_0 = eye(4);
gst_0(1:3, 4) = [L2+L3+L4, -1*L0, L1];
g1 = g_des*(gst_0^-1);
theta = zeros(1,6); %initializing vector of angles

%start passing stuff to subproblems to find joint angles
d1 = norm(g1*[Pw;1] - [Psh;1]);
theta(4) = SP3_solver(w3, Pw, Psh, d1, Pe_p);
    %note that this is actually theta_3, since MATLAB array indexing starts
        %at 1 and not at 0.
    %thus, theta(3) is actually theta_2 in my notes, etc.
d2 = norm(g1*[Pw;1] - [Psh_p;1]);
theta(3) = SP3_solver(w2, Pw, Psh_p, d2, Pw_p);

%define some exponential twists before we can continue
e2 = eye(4);
e2(1:3, 1:3) = Rz(-theta(3));
e2(1:3, 4) = Psh;
e3 = eye(4);
e3(1:3, 1:3) = Rz(-theta(4));
e3(1:3,4) = Pe;

P01_arg = e2*e3*[Pw;1];
P01 = P01_arg(1:3);
q01_arg = g1*[Pw;1];
q01 = q01_arg(1:3);


[theta(1), theta(2)] = SP2_solver(w0, w1, P01, q01, Psh);

%define a couple more exponential twists
e0 = eye(4);
e0(1:3, 1:3) = Rz(theta(1));
e0(1:3, 4) = Psh;
e1 = eye(4);
e1(1:3, 1:3) = Ry(-1*theta(2));
e1(1:3, 4) = Psh;

g2 = (e0*e1*e2*e3)^-1 * g1;

[theta(5), theta(6)] = SP2_solver(w4, w5, Pe, g2*Pe, Pw);

%account for offset in joint angles
offset = [-0.003 -0.002 0 0 0 -1.571];
output_thetas = theta - offset;


end



%solving inverse kinematics subproblems
function theta = SP1_solver(w, p, q, r)
%solves subproblem 1.
    %w = "omega" = axis of rotation.
    %p = "starting point" - the point before the rotation.
    %q = "ending point" - the point after the rotation.
    %r = another point you need to define vectors U and V. r is an
        %arbitrary point which is along the axis of rotation.
    %all of the above arguments are 3 x 1 vectors.
    u = p - r;
    v = q - r;
    uP = u - w*w.'*u; %"u-prime"
    vP = v - w*w.'*v; %"v-prime"
    %all of the above are 3x1 vectors
    
    Y = w.'*(cross(uP,vP));
    X = uP.' * vP;
    
    theta = atan2(Y,X);
end
function [theta1, theta2] = SP2_solver(w1, w2, p, q, r)
%Solves subproblem 2, for intersecting axes.
    %w1 = first axis of rotation.
    %w2 = second axis of rotation.
    %p = "starting position."
    %q = "ending position."
    %r = position where the axes intersect.

    %C = position after rotation 1 but before rotation 2.
    %Z = vector from r to C.
    I = eye(3);
    u = p-r;
    v = q-r;
    %Z = C - r, Z and C unknown
    %solve for alpha, beta, and gamma to calculate Z:
    %Z = alpha*w1 + beta*w2 + gamma*(w1 x w2)
    beta = w2.'*u * (I - (w1.'*w2)*(w1.'*w2))^(-1);
    alpha = w1.'*v - (w1.'*w2)*beta;
    top = (norm(u))^2 - alpha^2 - beta^2 -2*alpha*beta*(w1.'*w2);
    bottom = (norm(cross(w1,w2)))^2;
    gamma = sqrt(top/bottom);
    Z = alpha*w1 + beta*w2 + gamma*(cross(w1,w2));
    C = Z + r;
    
    %call to subproblem 1 solver:
    theta2 = SP1_solver(w2, p, C, r);
    theta1 = SP1_solver(-1*w1, q, C, r);
end
function theta = SP3_solver(w, p, q, d, r)
%Solves subproblem 3.
    %w = axis of rotation
    %p = starting position
    %q = "ending" position - the position you are trying to be distance
        %d away from.
    %d = "delta" = distance
    %w, p, and q are all 3x1 vectors. d is a scalar.

    u = p - r
    v = q - r
    uP = u - w*w.'*u %"u-prime"
    vP = v - w*w.'*v %"v-prime"
    dP2 = norm(d)^2 - norm(w.'*w*(q-p))^2
    theta_0 = atan2(w.'*cross(uP,vP), uP.'*vP)
    top = norm(uP)^2 + norm(vP)^2 - dP2
    bottom = 2*norm(uP)*norm(vP)
    phi = acos(top/bottom)    %THIS IS COMING OUT COMPLEX BECAUSE cos-1(top/bottom)
                                %IS COMPLEX IF TOP > BOTTOM
    
    theta = theta_0 + phi
    %theta = theta_0 - phi is also mathematically valid
    

end

%defining rotation matrices
function R = Rx(theta)
    R = eye(3);
    R(2:3, 2:3) = [cos(theta) -sin(theta); sin(theta) cos(theta)];
end
function R = Ry(theta)
    R = eye(3);
    R(1,:) = [cos(theta) 0 sin(theta)];
    R(3,:) = [-sin(theta) 0 cos(theta)];
end
function R = Rz(theta)
    R = eye(3);
    R(1:2, 1:2) = [cos(theta) -sin(theta); sin(theta) cos(theta)];
end


%NB: g_des = [-0.0256 -0.4496 -0.8929 -4.197;...
% 0.9758 0.1829 -0.12 15.369;...
% 0.2173 -0.8743 0.430 13.931;...
% 0 0 0 1]
