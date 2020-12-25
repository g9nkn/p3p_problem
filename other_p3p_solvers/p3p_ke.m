% A MATLAB code of Ke's p3p solver [1] implemented by G. Nakano.
% To use this code, follow the licence term of the original code.
% 
% Compatibility issue: 
%   1. A quartic equation solver is replaced with SolveQuartic.m for a fiar
%   comparison with other solvers.
%
% USAGE: 
%   [R, t] = p3p_ke(m, X)
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
%   polishing - (optional) a boolean whether to use the root polishing.
%               If false, the root polishing is not performed. (default: true)
% OUTPUS:
%   R - 3x3x4 rotation matrix. 
%       R(:,:,i) corresponds to t(:,i). 
%   t - 3x4 translation vector.
%
% REFERENCE:
%   [1] T. Ke and S. I. Roumeliotis, "An Efficient Algebraic Solution to the
%   Perspective-Three-Point Problem," CVPR2017.
%   https://github.com/opencv/opencv/blob/4.1.0/modules/calib3d/src/ap3p.cpp#L157

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
function [R, t] = p3p_ke(m, X, polishing)

    narginchk(2,3);
    if nargin < 3
        polishing = true;
    end

    m = m./sqrt(sum(m.^2));
    
    w1 = X(:,1);
    w2 = X(:,2);
    w3 = X(:,3);
    
    u0 = w1 - w2;
    nu0 = norm(u0);
    k1 = u0 / nu0;
    
    b1 = m(:,1);
    b2 = m(:,2);
    b3 = m(:,3);

    k3 = cross(b1,b2);
    nk3 = norm(k3);
    k3 = k3/nk3;
    
    tz = cross(b1,k3);
    v1 = cross(b1,b3);
    v2 = cross(b2,b3);
    
    u1 = w1 - w3;
    u1k1 = u1'*k1;
    k3b3 = k3'*b3;
    
    f11 = k3'*b3;
    f13 = k3'*v1;
    f15 = -u1k1*f11;
    nl = cross(u1,k1);
    delta = norm(nl);
    nl = nl/delta;
    f11 = delta * f11;
    f13 = delta * f13;
    
    u2k1 = u1k1 - nu0;
    f21 = tz'*v2;
    f22 = nk3 * k3b3;
    f23 = k3'*v2;
    f24 = u2k1 * f22;
    f25 = -u2k1 * f21;
    f21 = delta * f21;
    f22 = delta * f22;
    f23 = delta * f23;
    
    g1 = f13 * f22;
    g2 = f13 * f25 - f15 * f23;
    g3 = f11 * f23 - f13 * f21;
    g4 = -f13 * f24;
    g5 = f11 * f22;
    g6 = f11 * f25 - f15 * f21;
    g7 = -f15 * f24;
    alpha = [g5 * g5 + g1 * g1 + g3 * g3
             2 * (g5 * g6 + g1 * g2 + g3 * g4)
             g6 * g6 + 2 * g5 * g7 + g2 * g2 + g4 * g4 - g1 * g1 - g3 * g3
             2 * (g6 * g7 - g1 * g2 - g3 * g4)
             g7 * g7 - g2 * g2 - g4 * g4];

         
    sols = solveQuartic(alpha);
    if polishing
        sols = rootpolishing(alpha, sols);
    end
    
    Ck1nl = [k1, nl, cross(k1,nl)];
    Cb1k3tzT = [b1, k3, tz]';
    b3p = (delta/k3b3) * b3;
    
	R = zeros(3,3,4);
    t = zeros(3,4);
    for i = 1:4
        
        ctheta1p = real( sols(i) );
        if abs(ctheta1p) > 1
            continue
        end
        stheta1p = sqrt(1 - ctheta1p * ctheta1p);
        if k3b3 < 0
            stheta1p =  -stheta1p;
        end
        ctheta3 = g1 * ctheta1p + g2;
        stheta3 = g3 * ctheta1p + g4;
        ntheta3 = stheta1p / ((g5 * ctheta1p + g6) * ctheta1p + g7);
        ctheta3 = ntheta3 * ctheta3;
        stheta3 = ntheta3 * stheta3;

        C13 =[ctheta3,            0,         -stheta3
              stheta1p * stheta3, ctheta1p,  stheta1p * ctheta3
              ctheta1p * stheta3, -stheta1p, ctheta1p * ctheta3];
        
        R(:,:,i) = ( Ck1nl * C13 * Cb1k3tzT )';
        pxstheta1p = stheta1p * b3p;
        t(:,i)   = pxstheta1p - R(:,:,i)*w3; 
    end
end

function x = rootpolishing(alpha, x, maxitr)
    if nargin < 3
        maxitr = 2;
    end
    
    for i = 1:maxitr
        fx  = alpha(1)*x.^4 + alpha(2)*x.^3 + alpha(3)*x.^2 + alpha(4)*x + alpha(5);
        dfx = 4*alpha(1)*x.^3 + 3*alpha(2)*x.^2 + 2*alpha(3)*x + alpha(4);
        
        x = x - fx./dfx;
    end

end

