addpath('other_p3p_solvers')
solvers = {@p3p_gao, @p3p_kneip, @(m,X)p3p_ke(m,X,true), @(m,X)p3p_lambdatwist(m,X,true), @(m,X)p3p_nakano_bmvc2019(m,X,1)};
methods = {'Gao', 'Kneip', 'Ke', 'LambdaTwist', 'Nakano'};


npts  = 3;


% projective depth
min_depth = 5;
max_depth = 10;
d = min_depth + (max_depth - min_depth)*rand(1, npts);

% image points
w = 640;
h = 480;
f = 800;
K = [f, 0, w/2
     0, f, h/2
     0, 0,  1];
pts2d = [w*rand(1, npts)
         h*rand(1, npts)
         ones(1, npts)];
normalized_pts2d = K\pts2d;
     
% 3D points d*m = K*(R*X+t) -> X = R'*(d*K\m - t)
[U,~,V] = svd(rand(3));
R       = U * diag([1,1,det(U*V')]) * V';
t       = min_depth/2 * rand(3,1);
pts3d   = R' * (d.*normalized_pts2d - t);



% run solvers
for i=1:length(solvers)
    normalized_pts2d = normalized_pts2d./vecnorm(normalized_pts2d);
    [R_est, t_est] = solvers{i}(normalized_pts2d, pts3d);
    R_err = calc_R_err(R, R_est);
    t_err = calc_t_err(t, t_est);
    disp([ methods{i} '   --> R_err: ' num2str(min(R_err)) ', t_err: ' num2str(min(t_err))])
end



function R_err = calc_R_err(R_gt, R_est)
    for i=1:size(R_est,3)
        R_err(i) = norm(R_gt'*R_est(:,:,i) - eye(3), 'fro');
    end
end

function t_err = calc_t_err(t_gt, t_est)
    for i=1:size(t_est,2)
        t_err(i) = norm(t_gt - t_est(:,i)) / norm(t_gt) * 100;
    end
end