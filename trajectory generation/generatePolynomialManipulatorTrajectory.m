function [trajectories, parameters] = generatePolynomialManipulatorTrajectory(polyOrder, initCondMat, finalCondMat, timeSpan)
%generatePolynomialManipulatorTrajectory A function to generate a polynomial 
%trajectory for rigid manipulators
%   This function take as input:
%
%           polyOrder: the order of the polygon of the trajectory. It must
%           be larger than the number of the initial and final constraints
%           - 1.
%
%           initCondMat: inital condition values of q and its derivatives.
%           Each column corresponds to a joint, and each row corresponds to
%           the resepective derivatives or the q(s) themselves.
%
%           finalCondMat: final condition values of q and its derivatives.
%           Each column corresponds to a joint, and each row corresponds to
%           the resepective derivatives or the q(s) themselves.
%
%           timeSpan: the final time value as a constraint.
%
%   This function returns as output: 
%
%           trajectories: a matrix that contains trajectories for each
%           joint in terms of time. The columns correspond to a desired 
%           trajectory and its desired input. The rows correspond to y, dy,
%           ddy, and u of each trajectory.
%           
%           parameters: the parameters of the desired trajectory, instead
%           of functions in terms of time. Each column is a set of parameters
%           for the trajectory of each joint.
%


    %% check if the constraints are input correctly
    if size(initCondMat,2) == 0 || size(finalCondMat,2) == 0
        disp('You should enter at least one initial and final condition')
        return;
    end

    % if size(initCondMat,1) > 1 && size(finalCondMat,1) > 1
    %     disp('You should enter row vectors for the initial and final conditions')
    %     return;
    % end

    if size(initCondMat,2) ~= size(finalCondMat,2) || size(initCondMat,1) ~= size(finalCondMat,1)
        disp('You should have the same number of initial and final conditions')
        return;
    end

    if polyOrder + 1 < 2*size(initCondMat,1)
        disp('You should have the at least as many bases as constraints, so that you do not over constrain the system; (i.e. polyOrder + 1 >= size(initCondMat,1) + size(finalCondMat,1))')
        return;
    end

    %% variables

    % x = [q1, q2, q3, ..., qn, qdot1, qdot2, qdot3, ..., qdotn, ...]
    % conclude the number of joints from the size of the constraints
    % the number of columns of the constraints is the number of joints
    n = size(initCondMat, 2);
    % the number of rows of the constraints is the number of derivatives of qi,
    % qidot, qiddot, ..., etc, and qi itself.
    cond = size(initCondMat, 1);

    %% build bases vector based on the polynomial order

    % declare the symbolic time variable t
    syms t
    % allocate a row vector of bases with n + 1 variables
    bases = sym('basis',[1,polyOrder + 1]);
    % build the bases vector
    basis = 1;
    for i = 1:polyOrder + 1
        bases(i) = basis;
        basis = basis * t;
    end

    %% make the system

    % build the A matrix
    A = zeros(cond*n, cond*n);
    for i = 1:cond
        % build the ith column of x
        iColumn = [zeros(n,i*n)];
        if i < cond
            iColumn = [iColumn,eye(n,n)];
            iColumn = [iColumn,zeros(n,n*(cond-i-1))];
        end
        % add the ith column to the matrix
        A((n * (i - 1)) + 1 : (n * (i - 1)) + n,:) = iColumn;
    end

    % build the B matrix
    B = [zeros((cond - 1)*n,n);eye(n,n)];

    %% compute the bases coefficients
    % consider the case of underconstraints
    % build the bases matrix needed for the computation
    basesMatrix = sym('p',[cond,polyOrder + 1]);
    for i = 1:cond
        basesMatrix((i - 1) + 1,:) = diff(bases,i-1);
    end

    % rename the timeSpan
    T = timeSpan;

    % prepare the output
    trajectories = [];
    parameters = [];

    % compute the trajectories
    for joint = 1:n
        % allocate the coefficient vector
        b = sym('b',[polyOrder + 1,1]);

        % rename the initial and final conditions
        z0 = initCondMat(:,joint);
        zf = finalCondMat(:,joint);

        % build the equation sets for the initial and final conditions
        eq1 = subs(basesMatrix, t, 0)*b;
        eq2 = subs(basesMatrix, t, T)*b;

        % solve the equations
        solution = solve([eq1 - z0, eq2 - zf],b);

        % build the coefficient vector
        bVec = zeros(polyOrder + 1, 1);
        for i = 1:polyOrder + 1
            % programatically access the content of the solution structure
            bVec(i) = getfield(solution,strcat('b',num2str(i)));
        end

        % compute the desired trajectory polynomials
        y = basesMatrix(1,:)*bVec;
        dy = diff(y,t);
        ddy = diff(dy,t);

        % compute the desired input
        % manipulators are already in the control canonical form with the last row all zeros
        u = ddy;

        % add the resulting coefficients to the output
        desiredTraj = [y;dy;ddy;u];
        trajectories = [trajectories, desiredTraj];
        parameters = [parameters, bVec]; 
    end

end

