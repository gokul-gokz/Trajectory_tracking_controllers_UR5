function T = dhparam2matrices(theta, d, a, alpha, flag)
%dhparam2matrix returns the DH parameters defined homogeneous
%transformation
%   inputs nx1 DH parameter vectors, theta, d, a, alpha, and computes the
%   transformation

% get the number of links
n = size(theta,1);
n1 = size(theta,1);
n2 = size(d,1);
n3 = size(a,1);
n4 = size(alpha,1);

if (n1 == n2) && (n2 == n3) && (n3 == n4)
    % continue into the function
else
    disp("please provide same size nx1 vectors as params");
    return;
end

% initialize T
Ti = eye(4,4);

% check if the thetas are symbolic or not
if prod(size(class(theta), 2) == size('sym', 2)) == 1
    T = sym('a',[4, 4, n]);
else
    T = zeros(4,4,n);
end


% a loop that goes through all n links and generates the corresponding
% transformation matrix
    for i = 1:n
        %Ai = rot_theta * tran_z * tran_x * rot_alpha;
        Ai = [cos(theta(i)), -sin(theta(i))*cos(alpha(i)), sin(theta(i))*sin(alpha(i)), a(i)*cos(theta(i)); 
            sin(theta(i)), cos(theta(i))*cos(alpha(i)), -cos(theta(i))*sin(alpha(i)), a(i)*sin(theta(i)); 
            0, sin(alpha(i)), cos(alpha(i)), d(i); 
            0, 0, 0, 1];
        if flag == true
            disp('A = ');
            disp(Ai);
        end
        
        % compute output
        Ti  = Ti  * Ai;
        if flag == true
            disp('Ti = ');
            disp(Ti);
        end
        
        % store intermediate results
        T(:,:,i) = Ti;
    end
end

