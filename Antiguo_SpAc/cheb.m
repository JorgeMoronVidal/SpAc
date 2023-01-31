% Copyright (C) 2020 Sandeep Kumar
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.

% -*- texinfo -*- 
% @deftypefn {} {@var{retval} =} prg1 (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: Sandeep Kumar <sk@Sandeeps-MacBook-Pro.local>
% Created: 2020-11-10

% ------------------------------------------------------------------------
% cheb.m - Chebyshev differentiation matrix and grid points
% Aim: To compute the Chebyshev differentiation matrix and grid points for 
%      a given N value  
% ------------------------------------------------------------------------
 function [D,x] = cheb(N, R) 
   if N==0 
     D = 0; x =1 ; 
     return
   end
   x = R*cos(pi*(0:N)/N).' ; 
   c = [2; ones(N-1,1); 2] .* (-1) .^ (0:N)' ;
   X = repmat(x,1,N+1) ; 
   dX = X - X.' ;
   % The off-diagnoal entries 
   D = (c* (1 ./ c)') ./ (dX + eye(N+1)) ; 
   % diagnoal entries 
   D = D - diag(sum(D.')) ; 
 end

