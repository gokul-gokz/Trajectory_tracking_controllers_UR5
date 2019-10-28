if exist('q1')
    X = zeros(3, size(q1, 1));
    for i = 1:size(q1,1)
        q = [q1(i,2) q2(i,2) q3(i,2) q4(i,2) q5(i,2) q6(i,2)];
        X(:,i) = computeFK(q.');
    end
end

plot3(X(1,:), X(2,:), X(3,:));