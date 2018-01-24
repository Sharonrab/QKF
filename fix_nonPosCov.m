%https://www.mathworks.com/matlabcentral/answers/320134-make-sample-covariance-correlation-matrix-positive-definite
%Repair non-Positive Definite Correlation Matrix
function fixed_P = fix_nonPosCov(P)
% [V,D] = eig(P);       % Calculate the eigendecomposition of your matrix (A = V*D*V') 
%                         % where "D" is a diagonal matrix holding the eigenvalues of your matrix "A"
% d= diag(D);           % Get the eigenvalues in a vector "d" 
% d(d <= 1e-7) = 1e-7;  % Set any eigenvalues that are lower than threshold "TH" ("TH" here being 
% d(imag(d) <= 1e-7) = 1e-7;                        % equal to 1e-7) to a fixed non-zero "small" value (here assumed equal to 1e-7)
% D_c = diag(d);        % Built the "corrected" diagonal matrix "D_c"
% fixed_P = V*D_c*V';      % Recalculate your matrix "A" in its PD variant "A_PD"
P(abs(P) <= 1e-7) = 0;
fixed_P=P;