function [S_tabasco, be_bst, k_bst, Wbst, betas] = tabasco(X,vec_k,method,type,SCM,mu_known)
% TABASCO estimator (TApered or BAnded Shrinkage COvariance estimator) 
% 
% TABASCO is a shrinkage covariance matrix estimator that shrinks the 
% tapered  or banded sample covariance matrix (SCM), defined by W o S, 
% towards a scaled identity matrix [trace(S)/p]*I, where 
% S is the SCM and W is the tapering or banding matrix. 
% 
% The method chooses the minimum mean squared error (MMSE) optimal shrinkage 
% parameter, beta, and the bandwidth parameter, k, of the tapering
% or banding matrix data adaptively among all possible choises specified in
% the input parameter vec_k.  
%
% INPUT: 
% 
% vec_k     -  a vector of bandwidths.  A vector with possible values ranging 
%              from 1 to p, where p is the dimension of the data. Default 
%              is vec_k = 1:p if type is 'band' and vec_k=2:2:p if type is 
%              'taper'
% method    - 'ell1' or 'ell2'. Default is 'ell1' which means that 
%              Ell1-estimator of sphericity parameter is used.  
% type      - 'band' or 'taper'. Default is 'band' which implies using the
%              banding matrix W(i,j) = 1 if |i - j|<k and 0 otherwise. 
%              'taper' implies using the tapering matrix specified in 
%               equation (4) of our paper (Ollila and Breloy, 2021).
% SCM       -  pxp sample covariance matrix (SCM) of the data set X. If you
%              have already computed the SCM, you can pass it as an
%              optional argument to the function. 
% mu_known  -  true or false (logical). Default is false which implies that
%              the population mean vector of the data set is unknown, and
%              needs to be estimated. Hence the SCM is computed using the 
%              formula
%                   S = 1/(n-1) * SUM_i (x_i - bar_x)(x_i - bar_x)^T, 
%              where barx is the sample mean vector. 
%
% OUTPUT
% 
% S_tabasco  - computes TABASCO estimator (a symmetric pxp matrix)
% be_bst     - the computed MMSE optimal shrinkage parameter (a real scalar
%              in the interval [0,1])
% k_bst      - the computed MMSE optimal bandwidth. Equal to one of the
%              elements in the input vector vec_k. 
% Wbst       - the optimal tapering or banding matrix corresponding to the
%              chosen optimal bandwidth parameter k_bst
% betas      - optimal MMSE shrinkage parameter for each bandwidth in
%              vector vec_k. Thus betas is a vector of same size as vec_k.
% 
% USAGE (examples)
%
% Stab = tabasco(X,1:p,'ell1','band',SCM,true);
% [Stab, be_bst ] = tabasco(X,3,'ell1','taper');
% [Stab, be_bst,k_bst] = tabasco(X);
%
% REFERENCES:
% ------------
% If you use this code in your research, then please cite:
%
% [1] E. Ollila and Arnaud Breloy, "Regularized Tapered or Banded Sample 
%       Covariance Matrix", ArXiv preprint, submitted for publication, 
%       Sept. 2021.
%
% AUTHORS
% -------
% Esa Ollila (Aalto University), Arnaud Breloy (Universite Paris Nantarre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[n,p] = size(X);

assert((isreal(X) || iscomplex(X)) && (numel(size (X))==2),['''X'' must be a ' ...
    'real or complex matrix having at least two columns (variables)']);

if any (any (isnan (X)))
    error ('Input data contains NaN''s.');
end

if ~isa (X, 'double')
    X = double(X);
end

isrealX = isreal(X);

if nargin < 6
    mu_known = false;
end

if nargin < 3 || isempty(method)
    method = 'ell1';
end

assert(any(strcmp(method,{'ell1','ell2'})), ['''method'' ' ...
    'must be string equal to ''ell2'' or ''ell1''']);

if strcmpi(method,'ell1') 
    [SSCM,~,d] = SpatialSCM(X,mu_known);
    d_SSCM = diag(SSCM);
end

if ~mu_known
    xbar = mean(X);
    X = bsxfun(@minus,X,xbar);
    scaling = n-1;
else
    scaling = n;
end
       
if nargin < 5 || isempty(SCM)
    SCM = X'*X/scaling;
end

assert((isrealX==isreal(SCM)) && issymmetric(SCM) && (size(SCM,1)==p),['''SCM'' must be a ' ...
    'real or complex symmetric p x p matrix']);

if nargin < 4 || isempty(type)
    type = 'band';
end

assert(any(strcmp(type,{'band','taper'})), ['''type'' ' ...
    'must be string equal to ''band'' or ''taper''']);

if nargin < 2 || isempty(vec_k)    
    if strcmp(type,'band')
        vec_k = 1:p;
    elseif strcmp(type,'taper')
        vec_k = 2:2:p;
    end
end

assert(isvector(vec_k) && all(vec_k>=1) && all(vec_k<=p)  ... 
    && isequal(round(vec_k),vec_k), ['''vec_k'' must be a vector with ' ... 
    'integer elements between 1 and p']);

%% compute statistics

%-- elliptical kurtosis  / eta / gamma / D_S
kappa = elliptical_kurtosis(X);
trS = trace(SCM);
eta = trS/p;
d_S = diag(SCM);
D_S = diag(d_S);

% Compute sphericity statistic estimator of V o Sigma, where V = 1*1' 
if strcmpi(method,'ell1')
   %if mu_known  
   %    gam = sphericity_ell1(norm(SSCM,'Fro')^2,1,n,p,[]); 
   %else
       gam = sphericity_ell1(norm(SSCM,'Fro')^2,1,n,p,d,isreal(X));
   %end
else   
   gam = sphericity_ell2(norm(SCM,'Fro')^2,trS^2,trS,n,p,kappa,mu_known,isrealX);
end

%% start the for loop over bandwidths 
betas = zeros(1,length(vec_k));
risk_bst = inf;

for j = 1:length(vec_k)
    
    %-- compute the banding / tapering matrix 
    k =  vec_k(j);
    if strcmpi(type,'band')
        W = toeplitz( [ ones(1,k) , zeros(1,p-k) ] ); % banding matrix
    else
        W = taper_matrix(p,k,'taper');
    end    
    
    %-- compute sphericity estimators of W o Sigma and V o Sigma 
    trDsW2 = d_S'*(W.*W)*d_S;
    if strcmpi(method,'ell1') 
        %if mu_known  
        %    gam_W = sphericity_ell1(norm(W.*SSCM,'Fro')^2,d_SSCM'*(W.*W)*d_SSCM,n,p,[]);
        %else
            gam_W = sphericity_ell1(norm(W.*SSCM,'Fro')^2,d_SSCM'*(W.*W)*d_SSCM,n,p,d,isreal(X));
        %end

        if strcmpi(type,'taper')   
         %   if mu_known 
         %       gam_V = sphericity_ell1(norm(sqrt(W).*SSCM,'Fro')^2,d_SSCM'*W*d_SSCM,n,p,[]);   
         %   else
                gam_V = sphericity_ell1(norm(sqrt(W).*SSCM,'Fro')^2,d_SSCM'*W*d_SSCM,n,p,d,isreal(X));  
         %   end
        end
    else
        gam_W = sphericity_ell2(norm(W.*SCM,'Fro')^2,trDsW2,trS,n,p,kappa,mu_known,isrealX);
        if strcmpi(type,'taper')   
            gam_V = compute_sphericity(norm(sqrt(W).*SCM,'Fro')^2,d_S'*W*d_S,trS,n,p,kappa,mu_known);   
        end
    end

    if strcmpi(type,'band')
        gam_V = gam_W;
    end 

    %-- compute the optimal beta
    if mu_known 
        c = 1;
    else
        c = n/(n-1);
    end
    num = n * (gam_V - 1);
    den1 = n * (gam_W - 1);
    if isreal(X) 
        A = trDsW2/(p*eta^2) - 1 + 2*gam_W - 2*gam/p;
        den = den1 +   c*(trDsW2/(p*eta^2) + gam_W - 2*gam/p) + kappa*A;
    else
        A = trDsW2/(p*eta^2) - 1 + gam_W - 2*gam/p;
        den = den1 +   c*(trDsW2/(p*eta^2) - 2*gam/p) + kappa*A; 
    end

    be_opt = num/den;
    be_opt = min(1,max(0,be_opt));

    betas(j) = be_opt;

    % compute the risk for determining the best beta
    risk_k = be_opt*(1-gam_V);

    if risk_k < risk_bst 
        risk_bst = risk_k;
        k_bst = k;
        be_bst = betas(j);
    end
           
end

if strcmpi(type,'band')
    Wbst= toeplitz( [ ones(1,k_bst) , zeros(1,p-k_bst) ] ); % banding matrix
else
    Wbst = taper_matrix(p,k_bst,'taper');
end    
    
S_tabasco = be_bst * (Wbst.*SCM) + (1-be_bst)*eta*eye(p);  
end

%% FUNCTION 
function kappahat = elliptical_kurtosis(X,print_info)
% ELLIPTICAL_KURTOSIS computes the estimate of the elliptical kurtosis 
% parameter of a p-dimensional distribution given the data set X. 
% Function assumes that the data X is already centered (by sample mean when  
% mean vector is not known or by the true mean vector when it is known).  
%
% kappahat = elliptical_kurttos(X,...)
%
% INPUT:
%   X               data matrix of size n x p (rows are observations)
% Optional inputs:
%   print_info      (logical) verbose flag. Default=false
%
% Adapted from toolbox: RegularizedSCM ({esa.ollila,elias.raninen}@aalto.fi)
%--------------------------------------------------------------------------

[n,p] = size(X);

if isreal(X) 
    ka_lb = -2/(p+2); % theoretical lower bound for the kurtosis parameter
else 
    ka_lb = -1/(p+1); % theoretical lower bound for kurtosis parameter
end

if nargin < 2 || isempty(print_info)
    print_info = false;
end

vari =  mean(abs(X).^2);

indx = (vari==0);
if any(indx)
    if print_info
        fprintf('elliptical_kurtosis : found a variable with a zero sample variance\n');
        fprintf('         ...ignoring the variable in the calculation\n');
    end
end

if isreal(X) 
    kurt1n = (n-1)/((n-2)*(n-3));
    g2 = mean(X(:,~indx).^4)./(vari(~indx).^2)-3;
    G2 = kurt1n*((n+1)*g2 + 6);
    kurtest = mean(G2);
    kappahat = (1/3)*kurtest;
else
    g2 = mean(abs(X(:,~indx)).^4)./(vari(~indx).^2)-2;
    kurtest = mean(g2);
    kappahat = (1/2)*kurtest;
end

if kappahat > 1e6
    error('ellkurt: something is worong, too large value for kurtosis\n');
end

if kappahat <= ka_lb + (abs(ka_lb))/40
      kappahat = ka_lb + (abs(ka_lb))/40;
end
end

%% FUNCTION sphericity_ell2
function gam_V = sphericity_ell2(tr2a,tr2b,trS,n,p,kappa,mu_known,isrealX)
% SPHERICITY_ELL2 estimator of V o Sigma based on an estimator  V o S
% where S is the sample covariance matrix and V is the symmetric template 
% matrix whose diagonals are equal to 1. 
%
% INPUTS:
%   tr2_a = || V o S ||^2 
%   tr2_b = d_S'*(V o V)*d_S
%   trS = trace(S)
%   kappa = elliptical kurtosis estimate 
%   n = number of samples 
%   p = dimension 

if mu_known 
    a = (1 + kappa)/(n+kappa);
    b = (n*(kappa+n))/((n-1)*(n+2+3*kappa));
else
    a =  (1/(n+kappa))*( n/(n-1) + kappa); 
    if isrealX 
        b =  ((kappa  + n)*(n-1)^2)/((n-2)*(3*kappa*(n-1) + n*(n+1)));
    else
        b = ( n*(n-1)^2*(kappa + n))/(2*kappa*n*(n^2 - 4*n + 3) - kappa^2*(n-1)^2 + n^2 *( n^2 - 2*n  - 1));
    end
end 
trS_sq = trS^2;
gam_V = p * b * ( tr2a/trS_sq - a * (tr2b/trS_sq) ); 
gam_V = min (p , max(1,gam_V));

end

%%
function gam_V = sphericity_ell1(tr2a,tr2b,n,p,d,isrealX)
% INPUTS:
%   tr2a = || V o SSCM ||^2 
%   tr2b = d_SSCM'*(V o V)*d_SSCM, d_SSCM = diag(SSCM)
%   n = number of samples 
%   p = dimension 
%   d = n x 1 vector with ith component given by || X(i,:) ||
%
% If d is given, then the sphericity estimator uses a 
% correction factor developed in C. Zou et. al.
% "Multivariate sign-based high-dimensional tests for
% sphericity,” Biometrika, vol. 101, no. 1, pp. 229–236,
% 2014, that improves the sphericity estimator when p/n is large.
% 

if nargin < 6 || isempty(isrealX)
    isrealX = true;
end

if ~isempty(d) %&& isrealX
    m3 = mean(d.^(-3));
    m2 = mean(d.^(-2));
    m1 = mean(1./d);
    ratio = m2/(m1^2);
    ratio3 = m3/(m1^3);
    delta = (1/n^2)*(2 - 2*ratio + ratio^2) + ...
       (1/n^3)*(8*ratio - 6*ratio^2 + 2*ratio*ratio3 - 2*ratio3);
else
    delta = 0;
end

gam_V = p*(n/(n-1))*( tr2a - (1/n) * tr2b );
gam_V  = gam_V - p*delta;
gam_V = min (p , max(1,gam_V));
end
%%
function [Csgn,muhat,d] = SpatialSCM(X,mu_known,muhat)
% SPATIALSCM computes the spatial sign covariance matrix.
%
% INPUTS:
%
%   X               data matrix of size n x p (rows are observations)
%                   can be complex-valued or real-valued data.
% OPTIONAL INPUTS:
%
%   is_centered     (logical) is the X already centered. Default=false.
%   muhat           1 x p vector (e.g., spatial median of the data) to be
%                   used for centering the data.
%
% OUTPUTS:
%
%   Csgn            Spatial sign covariance matrix (p x p matrix)
%   muhat           Spatial median (1 x p vector) or the given muhat vector
%   d               Euclidean lengths of (centered) observations.
%
% By E. Ollila and E. Raninen (2020)

print_info = false;
p = size(X,2);

if nargin < 2 || isempty(mu_known)
    mu_known = false;
end

if ~mu_known
    if print_info
        fprintf('centering the data');
    end
    if nargin < 3
        muhat = spatmed(X);
    else
        assert(isequal(size(muhat),[1 p]));
    end
    X = bsxfun(@minus,X,muhat);
else
  muhat = [];
  if print_info
      fprintf('Not centering the data');
  end
end

d = sqrt(sum(X.*conj(X),2));
X = X(d~=0,:); % eliminate observations that have zero length
n = size(X,1);
X = bsxfun(@rdivide, X,d(d~=0));
Csgn = X'*X/n; % Sign covariance matrix
d(d<1.0e-10)= 1.0e-10;
end

%%
function  smed = spatmed(X,print_info)
% SPATMED computes the spatial median of the data set X.
%
% INPUTS:
%   X               data matrix of size n x p (rows are observations)
%                   Can be complex- or real-valued. 
% OPTIONAL INPUTS:
%   print_info      (logical) verbose flag. Default=false
%
% modified from toolbox:
% RegularizedSCM available from http://users.spa.aalto.fi/esollila/regscm/
%
% By E. Ollila and E. Raninen (2020)

if nargin ==1
    print_info = false;
end

if ~islogical(print_info)
    error('Input ''print_info'' needs to be logical');
end

len = sum(X.*conj(X),2); 
X = X(len~=0,:);
n = size(X,1);

if isreal(X)
    smed0 = median(X);
else
    smed0 = mean(X);
end
norm0 = norm(smed0);

iterMAX = 500;
EPS = 1.0e-4;
%TOL = 1.0e-6;
TOL = 1.0e-10;

for iter = 1:iterMAX

   Xc = bsxfun(@minus,X,smed0);
   len = sqrt(sum(Xc.*conj(Xc),2)); 
   len(len<EPS)= EPS;
   Xpsi = bsxfun(@rdivide, Xc, len);
   update = sum(Xpsi)/sum(1./len);
   smed = smed0 + update;

   dis = norm(update)/norm0;
   %fprintf('At iter = %3d, dis=%.6f\n',iter,dis);

   if (dis<=TOL)
       break;
   end
   smed0 = smed;
   norm0 = norm(smed);

end

if print_info
   fprintf('spatmed::convergence at iter = %3d, dis=%.10f\n',iter,dis);
end
end

%%

function [W] = taper_matrix(p,k,type)
%TAPER_MATRIX Summary of this function goes here
%   Detailed explanation goes here

assert(k <= p);

switch type
    case 'band'
        W = toeplitz( [ ones(1,k) , zeros(1,p-k) ] );
    case 'taper'     
        if mod(k,2) 
             k = k + 1;
        end
        W =  toeplitz( [ ones(1,k/2+1) ,  2*(1-((k/2+1):(k-1))/k),  zeros(1,p-k) ] );        
end

end



