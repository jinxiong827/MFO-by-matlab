function f = initialize_variables(N, M, V, min_range, max_range,name)

%% function f = initialize_variables(N, M, V, min_tange, max_range) 
% This function initializes the chromosomes. Each chromosome has the
% following at this stage
%       * set of decision variables
%       * objective function values
% 
% where,
% N - Population size
% M - Number of objective functions
% V - Number of decision variables
% min_range - A vector of decimal values which indicate the minimum value
% for each decision variable.
% max_range - Vector of maximum possible values for decision variables.
%
%
%     Copyright (C) 2009 Aravind Seshadri
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

min = min_range;
max = max_range;

% K is the total number of array elements. For ease of computation decision
% variables and objective functions are concatenated to form a single
% array. For crossover and mutation only the decision variables are used
% while for selection, only the objective variable are utilized.

K = M + V;

%% Initialize each chromosome
% For each chromosome perform the following (N is the population size)
LHS=lhsdesign(N,V);
for i = 1 : N
    % Initialize the decision variables based on the minimum and maximum
    % possible values. V is the number of decision variable. A random
    % number is picked between the minimum and maximum possible values for
    % the each decision variable.
    for j = 1 : V
        f(i,j) = min(j) + (max(j) - min(j))*LHS(i,j);
    end
    % For ease of computation and handling data the chromosome also has the
    % vlaue of the objective function concatenated at the end. The elements
    % V + 1 to K has the objective function valued. 
    % The function evaluate_objective takes one chromosome at a time,
    % infact only the decision variables are passed to the function along
    % with information about the number of objective functions which are
    % processed and returns the value for the objective functions. These
    % values are now stored at the end of the chromosome itself.
    f(i,V + 1: K) = evaluate_objective(f(i,:),name);
end
