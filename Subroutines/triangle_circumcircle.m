function [R, I, r] = triangle_circumcircle(A, B, C, nb_samples, option_display)
%% triangle_circumcircle : function to compute and display the circumcircle of one given triangle
%
% Author : nicolas.douillet9 (at) gmail.com, 2022-2024.
%
% Syntax
%
% triangle_circumcircle(A, B, C);
% triangle_circumcircle(A, B, C, nb_samples);
% triangle_circumcircle(A, B, C, nb_samples, option_display);
% [R, I, r] = triangle_circumcircle(A, B, C, nb_samples, option_display);
%
% Description
%
% triangle_circumcircle(A, B, C) computes and displays the circumcircle of ABC triangle.
% triangle_circumcircle(A, B, C, nb_samples) uses nb_samples to draw the circle.
% triangle_circumcircle(A, B, C, nb_samples, option_display) displays the circle when option_display is set either to
% logical true or real numeric 1, and doesn't when it is set to logical false or real numeric 0.
% [R, I, r] = triangle_circumcircle(A, B, C, nb_samples, option_display) stores the results in [R, I, r] vector.
%
% See also CIRCUMCENTER INCENTER
%
% Input arguments
%
%       [Ax]
% - A = [Ay] : real column vector double. 2 <= numel(A) <= 3. One of the three ABC vertices. 
%       [Az]
%
%       [Bx]
% - B = [By] real column vector double. 2 <= numel(A) <= 3. One of the three ABC vertices.
%       [Bz]
%
%       [Cx]
% - C = [Cy] real column vector double. 2 <= numel(A) <= 3. One of the three ABC vertices.
%       [Cz]
%
% - nb_samples : integer scalar double. The number of samples to draw the
%                circumcircle. nb_samples >= 3.
%
% - option_display : logical *true(1) / false(0), to enable/disable the display mode.
%
%
% Output arguments
%
%       [- Rx -]
% - R = [- Ry -]: real matrix double. The circumcircle coordinates. size(R) = [size(A,1), nb_samples].
%       [- Rz -]
%
%       [Ix]
% - I = [Iy] : real column vector double. 2 <= numel(I) <= 3. The circumcircle centre.
%       [Iz]
%
% - r : real scalar double. the circumcircle radius.
%
%
% Example #1
% From a triangle of the 3D space
% A = 2*(rand(3,1)-0.5);
% B = 2*(rand(3,1)-0.5);
% C = 2*(rand(3,1)-0.5);
% nb_samples = 30;
% option_display = true;
% triangle_circumcircle(A,B,C,nb_samples,option_display);
%
% Example #2
% From a triangle of the 2D space
% A = 2*(rand(2,1)-0.5);
% B = 2*(rand(2,1)-0.5);
% C = 2*(rand(2,1)-0.5);
% [R,I,r] = triangle_circumcircle(A,B,C);
%% Input parsing
assert(nargin > 2, 'Not enought input arguments. Three points required to define one triangle.');
assert(nargin < 6, 'Too many input arguments.');
if nargin < 5
    
    option_display = false;
    
    if nargin < 4
       
        nb_samples = 60;
        
    end
    
end
assert(isequal(size(A),size(B),size(C)),'All inputs points must have the same size.');
assert(isequal(ndims(A),ndims(B),ndims(C),2),'All inputs points must have the same number of dimensions (2).');
assert(isreal(A) && isreal(B) && isreal(C),'All inputs points must contain real numbers only.');
dimension = numel(A);
assert(dimension > 1 && dimension < 4,'Input points must have 2 or 3 elements.');
%% Body
if dimension < 3 % one padding in 2D case    
    
    A = cat(1,A,1);
    B = cat(1,B,1);
    C = cat(1,C,1);
    
end
% ABC triangle director and normal vectors computation 
AB = (B-A)/norm(B-A);
AC = (C-A)/norm(C-A);
n = cross(AB,AC);
% AB and AC segment mediatrices computation
I_AB = 0.5*(A+B);
I_AC = 0.5*(A+C);
% Circle centre I computation
[D,u] = planes_intersection(AB,I_AB,AC,I_AC); 
I = line_plane_intersection(u,D,n,A);
% Circle radius r computation
r = norm(I-A);
% Compute the circle point coordinates
angle_step = 2*pi/nb_samples;
theta = linspace(0,2*pi-angle_step,nb_samples);
Cx = r*cos(theta);
Cy = r*sin(theta);
Cz = zeros(1,nb_samples);
if dimension > 2
    
    % Vector u to rotate around
    k = [0 0 1]';
    u = cross(k,n)/norm(cross(k,n));
    
    % Angle between k and u
    alpha = atan2(norm(cross(k,n)),dot(k,n));
    
    % 3D rotation matrix around u vector
    Rm = @(delta)[u(1,1)^2+cos(delta)*(1-u(1,1)^2) (1-cos(delta))*u(1,1)*u(2,1)-u(3,1)*sin(delta) (1-cos(delta))*u(1,1)*u(3,1)+u(2,1)*sin(delta);
                  (1-cos(delta))*u(1,1)*u(2,1)+u(3,1)*sin(delta) u(2,1)^2+cos(delta)*(1-u(2,1)^2) (1-cos(delta))*u(2,1)*u(3,1)-u(1,1)*sin(delta);
                  (1-cos(delta))*u(1,1)*u(3,1)-u(2,1)*sin(delta) (1-cos(delta))*u(2,1)*u(3,1)+u(1,1)*sin(delta) u(3,1)^2+cos(delta)*(1-u(3,1)^2)];
    
    R = (Rm(alpha) * cat(1,Cx,Cy,Cz))' + I';
    
else % if dimension == 2
    
    R = cat(1,Cx,Cy,Cz)' + I';
    
    % one simplifications in 2D case
    R = R(:,1:2); 
    I = I(1:2);
    
end
%% Display
if option_display
    
    figure    
    
    if dimension > 2
        
        line([A(1,1) B(1,1) C(1,1) A(1,1)],[A(2,1) B(2,1) C(2,1) A(2,1)],[A(3,1) B(3,1) C(3,1) A(3,1)],'Color',[1 0 0],'Linewidth',2), hold on;
        line([R(:,1); R(1,1)],[R(:,2); R(1,2)],[R(:,3); R(1,3)],'Color',[0 0 1],'Linewidth',2), hold on;    
        view(3);
    
    else % if dimension == 2
        
        line([A(1,1) B(1,1) C(1,1) A(1,1)],[A(2,1) B(2,1) C(2,1) A(2,1)],'Color',[1 0 0],'Linewidth',2), hold on;
        line([R(:,1); R(1,1)],[R(:,2); R(1,2)],'Color',[0 0 1],'Linewidth',2), hold on;    
        view(2);
        
    end
    
    axis equal, axis tight;    
    
end
end % triangle_circumcircle

%% planes_intersection subfunction
function [I, u, rc] = planes_intersection(n1, M1, n2, M2)
%
% Author : nicolas.douillet9 (at) gmail.com, 2019-2024.
d1 = -dot(n1,M1); % -a1*x1 - b1*y1 - c1*z1
d2 = -dot(n2,M2); % -a2*x2 - b2*y2 - c2*z2
u = cross(n1,n2);
if norm(u) == 0 % (M1,n1) = (M2,n2) or (M1,n1) // (M2,n2)
   
    if (dot(n1,M2) + d1) == 0 && (dot(n2,M1) + d2) == 0 % (a1*M2(1) + b1*M2(2) + c1*M2(3) + d1) == 0              
        
        I = M1;
        u = M2 - M1;
        rc = 2;
        
    else                
        
        I = [];
        u = [];
        rc = 0;
        
    end
    
else 
          
     dir = find((abs(u) == max(abs(u))));     
     dir = dir(1);
     
     % => the line does exist in this direction, and then it can be set to t = 0.
     
     switch dir
         
         case 1 % setting : x = 0
             
             dx0y = (n1(3)*d2 - n2(3)*d1); % c1*d2 - c2*d1
             dx0z = (n2(2)*d1 - n1(2)*d2); % b2*d1 - b1*d2
             
             xI = 0;           
             yI = dx0y/u(1); 
             zI = dx0z/u(1);
             
         case 2 % setting : y = 0
             
             dxy0 = (n1(3)*d2 - n2(3)*d1); % c1*d2 - c2*d1
             dy0z = (n2(1)*d1 - n1(1)*d2); % a2*d1 - a1*d2
             
             xI = -dxy0/u(2);
             yI = 0;
             zI = -dy0z/u(2);
             
         case 3 % setting : z = 0
             
             dxz0 = (n1(2)*d2 - n2(2)*d1); % b1*d2 - b2*d1
             dyz0 = (n2(1)*d1 - n1(1)*d2); % a2*d1 - a1*d2
             
             xI = dxz0/u(3);
             yI = dyz0/u(3);
             zI = 0;                         
             
     end
     
     I = zeros(size(M1));
     I(1) = xI;
     I(2) = yI;
     I(3) = zI;
     
     rc = 1;
     
end
end % planes_intersection
%% line_plane_intersection subfunction
function [I,rc] = line_plane_intersection(u, N, n, M)
%
% Author : nicolas.douillet9 (at) gmail.com, 2019-2024.
% Plane offset parameter
d = -dot(n,M);
% Specific cases treatment
if ~dot(n,u) % n & u perpendicular vectors
    
    if dot(n,N) + d == 0 % N in P => line belongs to the plane  
        
        I = M;
        rc = 2;
        
    else % line // to the plane  
        
        I = [];
        rc = 0;
        
    end
    
else
    
    % Parametric line parameter t
    t = - (d + dot(n,N)) / dot(n,u);
    
    % Intersection coordinates
    I = N + u*t;
    
    rc = 1;
    
end
end % line_plane_intersection
