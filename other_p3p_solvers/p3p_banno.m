% A MATLAB code of Banno's p3p solver [1] implemented by G. Nakano.
% To use this code, follow the licence term of the original code.
%
% Compatibility issue: 
%   1. A quartic equation solver is replaced with SolveQuartic.m for a fiar
%   comparison with other solvers.
%
% USAGE: 
%   [R, t] = p3p_banno(m, X)
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
%   [1] A. Banno, "A P3P problem solver representing all parameters as a
%   linear combination", Image Vision and Computing, 70 (2018) 55-62.
%   https://github.com/atsuhikobanno/p3p

% BSD 3-Clause License
% 
% Copyright (c) 2019, Atsuhiko Banno (AIST)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
function [R, t] = p3p_banno(m, X)

    basescale = 10;
    X = X/10;
    
    nx = X(:,2) - X(:,1);
    nx = nx/norm(nx);
    nz = cross( nx, X(:,3)-X(:,1) );
    nz = nz/norm(nz);
    ny = cross(nz,nx);
    Rg = [nx, ny, nz];
    Xg = Rg'*( X - X(:,1) );
    
    
    u0 = zeros(12,1); 
    u1 = zeros(12,1); 
    u2 = zeros(12,1); 
    
    u0(10) = 1;
    u1(11) = 1;
    u2(12) = 1;

    XY  = 1 / Xg(1,2) / Xg(2,3);
    XmY = (Xg(1,2)-Xg(1,3)) * XY;
    
    m = m ./ sqrt(sum(m.^2));
    
    u0(1:6) = [-m(1,1) / Xg(1,2)
               -m(1,1) * XmY
               -m(2,1) / Xg(1,2)
               -m(2,1) * XmY
               -m(3,1) / Xg(1,2)
               -m(3,1) * XmY];
    u0(7:9) = m(:,1);
    
    u0_1 = u0([1,3,5]);
    u0_2 = u0([2,4,6]);
    u0_3 = u0(7:9);
    
    u1(1:6) = [m(1,2) / Xg(1,2)
               -m(1,2) * Xg(1,3) * XY
               m(2,2) / Xg(1,2)
               -m(2,2) * Xg(1,3) * XY
               m(3,2) / Xg(1,2)
               -m(3,2) * Xg(1,3) * XY];

    u1_1 = u1([1,3,5]);
    u1_2 = u1([2,4,6]);
    
	u2([2,4,6]) = m(:,3) / Xg(2,3);
    u2_2 = u2([2,4,6]);
    
    c12 = m(:,1)'*m(:,2);
    c23 = m(:,2)'*m(:,3);
    c31 = m(:,3)'*m(:,1);
    
    y_32 = Xg(2,3)^2;
    r_12 = Xg(1,2) - Xg(1,3);
    f1 = y_32 - r_12^2;
    f2 = -2.0*(y_32 + r_12*Xg(1,3)) * c12;
    f3 = 2.0*Xg(1,2) * r_12*c31;
    f4 = y_32 - Xg(1,3)^2;
    f5 = 2.0*Xg(1,2) * Xg(1,3) * c23;
    f6 = - Xg(1,2)^2;

    g1 = r_12;
    g2 = (Xg(1,3) - r_12)*c12;
    g3 = -Xg(1,2) * c31;
    g4 = -Xg(1,3);
    g5 = Xg(1,2) * c23;

    f11 = f1*f1;
    f12 = f1*f2;
    f13 = f1*f3;
    f14 = f1*f4;
    f15 = f1*f5;
    f16 = f1*f6;

    f22 = f2*f2;
    f23 = f2*f3;
    f24 = f2*f4;
    f25 = f2*f5;
    f26 = f2*f6;

    f33 = f3*f3;
    f34 = f3*f4;
    f35 = f3*f5;
    f36 = f3*f6;

    f44 = f4*f4;
    f45 = f4*f5;
    f46 = f4*f6;

    f55 = f5*f5;
    f56 = f5*f6;

    f66 = f6*f6;

    g11 = g1*g1;
    g12 = g1*g2;
    g13 = g1*g3;
    g14 = g1*g4;
    g15 = g1*g5;
    
    g22 = g2*g2;
    g23 = g2*g3;
    g24 = g2*g4;
    g25 = g2*g5;

    g33 = g3*g3;
    g34 = g3*g4;
    g35 = g3*g5;

    g44 = g4*g4;
    g45 = g4*g5;

    g55 = g5*g5;
    
    h = [f44*g11 + f11*g44 - f24*g12 + f22*g14 + f14*g22 - f12*g24 - 2*f14*g14
         2*(f45*g11 + f11*g45 + f23*g14 + f14*g23 - f15*g14 - f14*g15) - f34*g12 - f25*g12 - f24*g13 + f22*g15 + f15*g22 - f13*g24 - f12*g34 - f12*g25
         2*(f46*g11 + f23*g15 + f15*g23 - f16*g14 - f15*g15) + f55*g11 + f11*g55 - f35*g12 - f26*g12 - f34*g13 - f25*g13 + f33*g14 + f16*g22 + f14*g33 - f13*g34 - f13*g25 - f12*g35
         2*(f56*g11 + f16*g23 - f16*g15) - f36*g12 - f35*g13 - f26*g13 + f33*g15 + f15*g33 - f13*g35
         f66*g11 - f36*g13 + f16*g33];

    sols = solveQuartic(h);
    
    n = length(sols);
    if n == 0
        R = [];
        t = [];
        
    else
        R = zeros(3,3,n);
        t = zeros(3,n);
        for i = 1:n
            y = real(sols(i));
            x = (f6*g1 + y*( f5*g1 - f1*g5 + y*(f4*g1 - f1*g4))) / ((f1*g2 - f2*g1)*y + f1*g3 - f3*g1);
            
            nrm = norm(x*u0_1+y*u1_1)^2 + norm(x*u0_2+y*u1_2+u2_2)^2;
            z = 1/sqrt(nrm/2);
            x = x*z;
            y = y*z;
            
            r1 = x*u0_1 + y*u1_1;
            r2 = x*u0_2 + y*u1_2 + z*u2_2;
            
            r1 = r1/norm(r1);
            r2 = r2/norm(r2);
            r3 = cross(r1,r2);
            r3 = r3/norm(r3);
            Rc = [r1, r2, r3]';
              
            tc = x*u0_3;
            
            R(:,:,i) = (Rg*Rc)';
            t(:,i)   = basescale * (tc - R(:,:,i)*X(:,1));            
        end
        
    end
    
end