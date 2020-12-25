% A MATLAB code of Gao's p3p solver [1] implemented by G. Nakano.
% To use this code, follow the licence term of the original code.
% 
% Compatibility issue: 
%   1. A quartic equation solver, solve_deg4() in solve_for_length(), is
%   replaced with SolveQuartic.m for a fiar comparison with other solvers.
%   2. An eigenvalue solver, jacobi_4x4() in align(), is replaced with
%   MATLAB's eig() for simplicity.
%
% USAGE: 
%   [R, t] = p3p_gao(m, X)
%   solves P3P problem given by m_i ~ R*X_i + t (i={1,2,3}). 
%
% INPUTS:
%   m - 3x3 matrix of 2D points represented by homogeneous coordinates.
%       Each column m(:,i) corresponds to the 3D point X(:,i),
%       [u1, u2, u3
%        v1, v2, v3
%        w1, w2, w3].
%
%   X - 3x3 matrix of 3D points.
%       Each column X(:,i) corresponds to the 2D point m(:,i),
%       [x1, x2, x3
%        y1, y2, y3
%        z1, z2, z3].
%
% OUTPUS:
%   R - 3x3x4 rotation matrix. 
%       R(:,:,i) corresponds to t(:,i). 
%   t - 3x4 translation vector.
%
% REFERENCE:
%   [1] X. Gao et al., "Completesolution classification for the
%   perspective-three-point problem", IEEE PAMI, 25(8):930-943, 2003.
%   https://github.com/opencv/opencv/blob/4.1.0/modules/calib3d/src/p3p.cpp#L132

% By downloading, copying, installing or using the software you agree to this license.
% If you do not agree to this license, do not download, install,
% copy or use the software.
% 
% 
%                           License Agreement
%                For Open Source Computer Vision Library
%                        (3-clause BSD License)
% 
% Copyright (C) 2000-2019, Intel Corporation, all rights reserved.
% Copyright (C) 2009-2011, Willow Garage Inc., all rights reserved.
% Copyright (C) 2009-2016, NVIDIA Corporation, all rights reserved.
% Copyright (C) 2010-2013, Advanced Micro Devices, Inc., all rights reserved.
% Copyright (C) 2015-2016, OpenCV Foundation, all rights reserved.
% Copyright (C) 2015-2016, Itseez Inc., all rights reserved.
% Third party copyrights are property of their respective owners.
% 
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
% 
%   * Redistributions of source code must retain the above copyright notice,
%     this list of conditions and the following disclaimer.
% 
%   * Redistributions in binary form must reproduce the above copyright notice,
%     this list of conditions and the following disclaimer in the documentation
%     and/or other materials provided with the distribution.
% 
%   * Neither the names of the copyright holders nor the names of the contributors
%     may be used to endorse or promote products derived from this software
%     without specific prior written permission.
% 
% This software is provided by the copyright holders and contributors "as is" and
% any express or implied warranties, including, but not limited to, the implied
% warranties of merchantability and fitness for a particular purpose are disclaimed.
% In no event shall copyright holders or contributors be liable for any direct,
% indirect, incidental, special, exemplary, or consequential damages
% (including, but not limited to, procurement of substitute goods or services;
% loss of use, data, or profits; or business interruption) however caused
% and on any theory of liability, whether in contract, strict liability,
% or tort (including negligence or otherwise) arising in any way out of
% the use of this software, even if advised of the possibility of such damage.
function [R, t] = p3p_gao(m, X)

	m = m./sqrt(sum(m.^2));
	
	distances(1) = norm(X(:,2) - X(:,3));
	distances(2) = norm(X(:,1) - X(:,3));
	distances(3) = norm(X(:,1) - X(:,2));
	
	cosines(1) = m(:,2)'*m(:,3);
	cosines(2) = m(:,1)'*m(:,3);
	cosines(3) = m(:,1)'*m(:,2);
	
	lengths = solve_for_length(distances, cosines);

	R = zeros(3,3,4);
	t = zeros(3,4);
	for i=1:4
	
		M_orig = zeros(3,3);
		M_orig(:,1) = lengths(i,1) * m(:,1);
		M_orig(:,2) = lengths(i,2) * m(:,2);
		M_orig(:,3) = lengths(i,3) * m(:,3);
		
		[R(:,:,i), t(:,i)] = align(M_orig, X);
	
	end
end


function lengths = solve_for_length(distances, cosines)

	lengths = zeros(4,3);

	p = cosines(1) * 2;
    q = cosines(2) * 2;
    r = cosines(3) * 2;

    inv_d22 = 1 / (distances(3) * distances(3));
    a = inv_d22 * (distances(1) * distances(1));
    b = inv_d22 * (distances(2) * distances(2));

    a2 = a * a; b2 = b * b; p2 = p * p; q2 = q * q; r2 = r * r;
    pr = p * r; pqr = q * pr;


    ab = a * b; a_2 = 2*a;
    A = -2 * b + b2 + a2 + 1 + ab*(2 - r2) - a_2;


    a_4 = 4*a;

    B = q*(-2*(ab + a2 + 1 - b) + r2*ab + a_4) + pr*(b - b2 + ab);
    C = q2 + b2*(r2 + p2 - 2) - b*(p2 + pqr) - ab*(r2 + pqr) + (a2 - a_2)*(2 + q2) + 2;
    D = pr*(ab-b2+b) + q*((p2-2)*b + 2 * (ab - a2) + a_4 - 2);
    E = 1 + 2*(b - a - ab) + b2 - b*p2 + a2;

    temp = (p2*(a-1+b) + r2*(a-1-b) + pqr - a*pqr);
    b0 = b * temp * temp;
    

    sols = solveQuartic([A, B, C, D, E]);
	n = length(sols);

    
    r3 = r2*r; pr2 = p*r2; r3q = r3 * q;
    inv_b0 = 1 / b0;

    % For each solution of x
    for i = 1:n
        x = real(sols(i));

        x2 = x*x;
        b1 =...
            ((1-a-b)*x2 + (q*a-q)*x + 1 - a + b) * ...
            (((r3*(a2 + ab*(2 - r2) - a_2 + b2 - 2*b + 1)) * x + ...
            (r3q*(2*(b-a2) + a_4 + ab*(r2 - 2) - 2) + pr2*(1 + a2 + 2*(ab-a-b) + r2*(b - b2) + b2))) * x2 + ...
            (r3*(q2*(1-2*a+a2) + r2*(b2-ab) - a_4 + 2*(a2 - b2) + 2) + r*p2*(b2 + 2*(ab - b - a) + 1 + a2) + pr2*q*(a_4 + 2*(b - ab - a2) - 2 - r2*b)) * x + ...
            2*r3q*(a_2 - b - a2 + ab - 1) + pr2*(q2 - a_4 + 2*(a2 - b2) + r2*b + q2*(a2 - a_2) + 2) + ...
            p2*(p*(2*(ab - a - b) + a2 + b2 + 1) + 2*q*r*(b + a_2 - a2 - ab - 1)));

        y = inv_b0 * b1;
        v = x2 + y*y - x*y*r;

        Z = distances(3) / sqrt(v);
        X = x * Z;
        Y = y * Z;

        lengths(i,1) = X;
        lengths(i,2) = Y;
        lengths(i,3) = Z;
    end

end


function [R, t] = align(M_end, X)
    % Centroids:
	C_end   = mean(M_end, 2);
	C_start = mean(X, 2);

    % Covariance matrix s:
	X0 = X(1,1); X1 = X(1,2); X2 = X(1,3);
	Y0 = X(2,1); Y1 = X(2,2); Y2 = X(2,3);
	Z0 = X(3,1); Z1 = X(3,2); Z2 = X(3,3);
    s = zeros(9,1);
    for j=1:3
        s(0 * 3 + j) = (X0 * M_end(j,1) + X1 * M_end(j,2) + X2 * M_end(j,3)) / 3 - C_end(j) * C_start(1);
        s(1 * 3 + j) = (Y0 * M_end(j,1) + Y1 * M_end(j,2) + Y2 * M_end(j,3)) / 3 - C_end(j) * C_start(2);
        s(2 * 3 + j) = (Z0 * M_end(j,1) + Z1 * M_end(j,2) + Z2 * M_end(j,3)) / 3 - C_end(j) * C_start(3);
    end
	
    Qs = zeros(16,1);
    Qs(0 * 4 + 1) = s(0 * 3 + 1) + s(1 * 3 + 2) + s(2 * 3 + 3);
    Qs(1 * 4 + 2) = s(0 * 3 + 1) - s(1 * 3 + 2) - s(2 * 3 + 3);
    Qs(2 * 4 + 3) = s(1 * 3 + 2) - s(2 * 3 + 3) - s(0 * 3 + 1);
    Qs(3 * 4 + 4) = s(2 * 3 + 3) - s(0 * 3 + 1) - s(1 * 3 + 2);
      
    Qs(1 * 4 + 1) = s(1 * 3 + 3) - s(2 * 3 + 2);
	Qs(0 * 4 + 2) = Qs(1 * 4 + 1);
    Qs(2 * 4 + 1) = s(2 * 3 + 1) - s(0 * 3 + 3);
	Qs(0 * 4 + 3) = Qs(2 * 4 + 1);
    Qs(3 * 4 + 1) = s(0 * 3 + 2) - s(1 * 3 + 1);
	Qs(0 * 4 + 4) = Qs(3 * 4 + 1);
    Qs(2 * 4 + 2) = s(1 * 3 + 1) + s(0 * 3 + 2);
	Qs(1 * 4 + 3) = Qs(2 * 4 + 2);
    Qs(3 * 4 + 2) = s(2 * 3 + 1) + s(0 * 3 + 3);
	Qs(1 * 4 + 4) = Qs(3 * 4 + 2);
    Qs(3 * 4 + 3) = s(2 * 3 + 2) + s(1 * 3 + 3);
	Qs(2 * 4 + 4) = Qs(3 * 4 + 3);
	
    %jacobi_4x4(Qs, evs, U);
	Qs = reshape(Qs,4,4);
	[U, evs] = eig(Qs);
	
    % Looking for the largest eigen value:
	[~, i_ev] = max(diag(evs));

    % Quaternion:
	q = U(:,i_ev);

    q02  = q(1) * q(1); q12  = q(2) * q(2); q22  = q(3) * q(3); q32 = q(4) * q(4);
    q0_1 = q(1) * q(2); q0_2 = q(1) * q(3); q0_3 = q(1) * q(4);
    q1_2 = q(2) * q(3); q1_3 = q(2) * q(4);
    q2_3 = q(3) * q(4);

	R = zeros(3,3);
    R(1,1) = q02 + q12 - q22 - q32;
    R(1,2) = 2 * (q1_2 - q0_3);
    R(1,3) = 2 * (q1_3 + q0_2);

    R(2,1) = 2 * (q1_2 + q0_3);
    R(2,2) = q02 + q22 - q12 - q32;
    R(2,3) = 2 * (q2_3 - q0_1);

    R(3,1) = 2 * (q1_3 - q0_2);
    R(3,2) = 2 * (q2_3 + q0_1);
    R(3,3) = q02 + q32 - q12 - q22;

	t = C_end - R*C_start;

end