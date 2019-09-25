function A_hat = findPSDS(A)
% FFNAV Extended Kalman Filter ============================================
% Description: This function takes a square matrix and determines the
% nearest positive semi-definite symmetric matrix, based on the Frobenius
% norm.
%
% Inputs:
%   A - a square matrix, which will be converted to the nearest PSDS matrix
%   
% Outputs:
%   A_hat - the nearest PSDS matrix to A
%
% References:
%   Higham - Computing a Nearest Symmetric Positive Semidefinite Matrix
%   D'Errico - "nearestSPD", MathWorks File Exchange
%
% Created by: Cory Fraser  - JUN 06, 2018
% Latest Edit: Cory Fraser - JUN 07, 2018
% Copyright(c) 2018 by Cory Fraser
% =========================================================================

%% Initialize Parameters
    coder.extrinsic('fprintf');
    coder.extrinsic('eig');

% Test if A is square
[r,c] = size(A);
if r ~= c
  error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
  % A was scalar and non-positive, so just return eps
  A_hat = eps;
  return
end

% =========================================================================
%% Nearest PSDS Matrix Algorithm 

% Step 1: Construct the symmetric B matrix, and skew-symmetric C matrix
B = (A + A')/2;
%C = (A - A')/2;    

% Step 2: Compute the symmetric polar factor of B
[U,Sigma,V] = svd(B);
H = V*Sigma*V';

% Step 3: Compute the positive approximant matrix, symmetricize
A_hat = (B+H)/2;
A_hat = (A_hat + A_hat')/2;

% Step 4: Testing for PSD, adjusting A_hat if not
p = 1;
k = 0;

while p ~= 0
  [~,p] = chol(A_hat);
  %p = any(real(eig(A_hat))<0); %Another possible test for PDS
  
  k = k + 1;
  if k > 10000
     fprintf('Forced break after 10000 iterations \n')
     break
  end
  
  if p ~= 0 % A_hat is not PSDS, perturb it a bit
    min_eig = 0;
    min_eig = min(eig(A_hat)); 
    
    %tol = eps;
    %min_eig_test = max(tol,min(eig(A_hat))); %New Check

    A_hat = A_hat + (-min_eig*k.^2 + eps(min_eig))*eye(size(A));
  end
end
     
end