function [X,E,Y,s,be_X,be_Y,tau] = quadeig(A,B,C,options)
%quadeig quadratic eigenvalue problem.
%   E = quadeig(A,B,C) is a vector of length 2*n whose elements are the
%   eigenvalues of the n-by-n quadratic matrix polynomial lambda^2*A +
%   lambda*B + C. The input is 3 square matrices, A, B, C, all of the same
%   order, n.
%
%   [X,E] = quadeig(A,B,C) solves the quadratic eigenvalue problem
%        (lambda^2*A + lambda*B + C)*x = 0.
%   The output is an n-by-2*n matrix, X, whose columns are the right
%   eigenvectors, and a vector of length 2*n, E, whose elements are the
%   eigenvalues.
%
%   [X,E,Y] = quadeig(A,B,C) solves the quadratic eigenvalue problems
%   (lambda^2*A + lambda*B + C)*x = 0 and y^T*(lambda^2*A + lambda*B + C) = 0.
%   The output is an n-by-2*n matrix, X, whose columns are the right
%   eigenvectors, a vector of length 2*n, E, whose elements are the
%   eigenvalues and an n-by-2*n matrix, Y, whose columns are the left
%   eigenvectors.
%      for j = 1:2*n
%         lambda = E(j)
%         x = X(:,j)
%         y = Y(:,j)
%         if lambda ~=Inf
%            (lambda^2*A + lambda*B + C)*x is approximately 0.
%            y'*(lambda^2*A + lambda*B + C) is approximately 0.
%         else
%             C*x and y'*C are approximately 0.
%      end
%
%   If both A and C are singular, or nearly singular, the problem is
%   potentially ill-posed.  Theoretically, the solutions might not exist
%   or might not be unique.  Computationally, the computed solutions may
%   be inaccurate.  If one, but not both, of A and C is singular, the
%   problem is well-posed, but some of the eigenvalues may be zero or
%   infinite.
%
%   quadeig(A,B,C,OPTIONS) performs the computation with the default
%   parameters replaced by values in the structure OPTIONS, which has
%   fields: PSCALE, TOL and NAG:
%
%     OPTIONS.PSCALE = 0: no parameter scaling.
%     OPTIONS.PSCALE = 1: Fan, Lin and Van Dooren scaling if
%                         norm(B)/sqrt(norm(A)*norm(C)) < 10 and
%                         no scaling otherwise (default).
%     OPTIONS.PSCALE = 2: Fan, Lin and Van Dooren scaling.
%     OPTIONS.PSCALE = 3: tropical scaling with largest root.
%     OPTIONS.PSCALE = 4: tropical scaling with smallest root.
%
%     OPTIONS.TOL is for rank decisions in the deflation procedure.
%     OPTIONS.TOL defaults to  n*EPS.
%     If OPTIONS.TOL = -1 then no deflation is attempted.
%
%     OPTIONS.NAG = 1: use MATLAB NAG Toolbox if available (default).
%     OPTIONS.NAG = 0: do not use the MATLAB NAG Toolbox.
%
%   [X,E,Y,S,BE_X,BE_Y] = quadeig(A,B,C) also returns
%   S:    a vector of length 2*n of condition numbers for the eigenvalues
%         (large condition numbers imply that the problem is near one with
%         multiple eigenvalues),
%   BE_X: a vector of length 2*n of backward errors for the computed right
%         eigenpairs and
%   BE_Y: a vector of length 2*n of backward errors for the computed left
%         eigenpairs.
%
%   [X,E,Y,S,BE_X,BE_Y,tau] = quadeig(A,B,C) also returns
%   S:    a vector of length 2*n of condition numbers for the eigenvalues
%         (large condition numbers imply that the problem is near one with
%         multiple eigenvalues),
%   BE_X: a vector of length 2*n of backward errors for the computed right
%         eigenpairs and
%   BE_Y: a vector of length 2*n of backward errors for the computed left
%         eigenpairs.
%   tau:  a scalar whose value is
%         tau = norm(C,'fro')/( norm(A,'fro')*norm(B,'fro') ),
%            when A and B are both non-zero and
%         tau = -1  otherwise.
%
%   See also eig, condeig, polyeig.

%   Chris Munro, Sven Hammarling and Francoise Tisseur.
%   $Revision: 4 $  $Date: 2014/09/30$
%   Revised to take advantage of the version of 'eig' introduced in
%   R2014a (version 8.3), which allows for the computation of left eigenvectors.
%   $Revision: 5 $  $Date: 2014/11/18$
%   Modified the tolerances that determine rank
%   Put in test for heavily damped systems
%      (norm(B)/sqrt(norm(A)*norm(C)) >= 10).
%%%%%%%%%%%%%%%%%%%%%
% Parameter settings
%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
  options = struct('pscale', 1, 'tol', [], 'nag', 1);
else
  % Set up default values for empty fields.
  fields = isfield(options,{'pscale', 'tol', 'nag'});
  if ~fields(1), options.pscale = 1; end
  if ~fields(2), options.tol = []; end
  if ~fields(3), options.nag = 1; end
end

n = length(A);
rev = 0;
Qreal = isreal([A B C]); % To perform computation in real arithmetic.

vt = ver('NAG'); % Check for NAG Toolbox
if isempty(vt)
  nagtoolbox = 0;
else
  nagtoolbox = min(strcmp(vt.Name, 'NAG Toolbox'), options.nag);
end
if nagtoolbox, b32 = ~isempty(strfind(vt.Release, '32')); end

% Before version 8.3 'eig' cannot compute left eigenvectors
if verLessThan('matlab','8.3') == 1
  eigleft = 0;
else
  eigleft = 1;
end

toltau = 10;

%%%%%%%%%%%
% Scalings
%%%%%%%%%%%
g2 = norm(A, 'fro'); % Frobenius norms of initial matrices
g1 = norm(B, 'fro');
g0 = norm(C, 'fro');

if options.pscale == 0
  d = 1;
else
  d = max([g0 g1 g2]);
  if d == 0
     d = 1; % All three matrices are the zero matrix
  end
end

g2 = g2/d; g1 = g1/d; g0 = g0/d;
temp = sqrt(g0*g2);
if temp ~= 0
  tau = g1/temp;
  if tau >= toltau
    warning('MATLAB:quadeig:HeavilyDamped', '%s', ...
            'the problem is heavily damped, ', ...
            'small backward errors cannot be guaranteed.')
  end
else
  tau = toltau;
end
nA = g2; nB = g1; nC = g0;

% Eigenvalue parameter scaling
if options.pscale

  if nC == 0 || nA == 0
    warning('MATLAB:quadeig:ZeroCoeffMatrix', '%s', ...
            'one or both of norm(A) and norm(C) is zero, ', ...
            'no parameter scaling performed.')
    gg = 1; delta = 1;

  elseif tau < toltau  || options.pscale == 2
    % Fan Lin and Van Dooren scaling

    gg = sqrt(nC/nA);
    delta = 2/(nC + nB*gg);

  elseif options.pscale == 3
    % tropical scaling largest root

    gg = nB/nA;
    delta = nA/(nB^2);

  elseif options.pscale == 4
    % tropical scaling smallest root

    gg = nC/nB;
    delta = 1/nC;

  else
    % options.pscale == 1  with tau >= 10

    gg = 1; delta = 1;

  end

else

  gg = 1; delta = 1;

end

if gg ~= 1 || delta ~= 1  || d ~= 1
  sb = gg*delta;
  sa = gg*sb;

  A = (sa/d)*A;
  B = (sb/d)*B;
  C = (delta/d)*C;
  % Frobenius norm of scaled A,B,C
  nA = sa*nA;
  nB = sb*nB;
  nC = delta*nC;
end

if isempty(options.tol)
  options.tol = n*eps;
end

%%%%%%%%%%%%%%%%%%%%
% Rank determination
%%%%%%%%%%%%%%%%%%%%
if issparse(A)
  A = full(A);
end
[Q2,R2,P2] = qr(A);
r2 = sum(abs(diag(R2)) > options.tol*max([1,nA]));
R2 = R2(1:r2,:);         % Throw away negligible rows

if issparse(C)
  C = full(C);
end
[Q0,R0,P0] = qr(C);
r0 = sum(abs(diag(R0)) > options.tol*max([1,nB,nC]));
R0 = R0(1:r0,:);         % Throw away negligible rows

rev = r0 > r2;
if rev  % Solve for reversal eigenproblem.
  temp = A; A = C; C = temp;
  temp = Q2; Q2 = Q0; Q0 = temp;
  temp = P2; P2 = P0; P0 = temp;
  temp = R2; R2 = R0; R0 = temp;
  temp = r2; r2 = r0; r0 = temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearization and deflation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if r0 == n &&  r2 == n % no deflation

  AA = [-Q2'*B*P2 Q2'; -C*P2 zeros(n)];
  BB = [R2 zeros(n); zeros(n) eye(n)];
  if nargout > 1
    V = [P2  zeros(n); zeros(n) eye(n)];
  end
  if nargout > 2
    Q = [Q2' zeros(n); zeros(n) eye(n)];
  end

elseif r0 < r2 && r2 ==n

  AA = [-Q2'*B*P2 Q2'*Q0; -R0*P0'*P2 zeros(r0,n); zeros(n-r0,2*n)];
  BB = [R2 zeros(n); zeros(n) eye(n)];
  if nargout > 1
    V = [P2  zeros(n); zeros(n) Q0 ];
  end
  if nargout > 2
    Q = [Q2' zeros(n); zeros(n) Q0'];
  end

else

  warning('MATLAB:quadeig:SingularCoeffMatrices', '%s',  ...
          'A and C are singular, or nearly singular ',...
          'the problem is potentially ill-posed.')
  AA = [-Q2'*B Q2'*Q0; -R0*P0' zeros(r0,n); zeros(n-r0,2*n)];
  BB = [R2*P2' zeros(r2,n); zeros(n-r2,2*n); zeros(n) eye(n)];
  [Q3, R3, P3] = qr(AA(r2+1:n,1:n+r0));

  % Determine r3 = effective rank.
  r3 = sum(abs(diag(R3)) > options.tol*norm(AA, 'fro'));
  r3 = r3(1);        % Fix for case where R is a vector.
  if r3 < n-r2
    error('MATLAB:quadeig:SingularQuadratic', ...
          'The quadratic matrix polynomial is nonregular');
  end

  R3 = R3(1:r3,:);   % Throw away negligible rows (incl. all zero rows, m > n).
  if r3 ~= n+r0
    % Reduce R3 to upper triangular form.
    
    if nagtoolbox

      if b32
        mp0 = int32(n+r0);
        mm0 = int32(n-r0);
        mm2 = int32(n-r2);
        ir3 = int32(r3);
        r0r2 = int32(r0+r2);
      else
        mp0 = int64(n+r0);
        mm0 = int64(n-r0);
        mm2 = int64(n-r2);
        ir3 = int64(r3);
        r0r2 = int64(r0+r2);
      end

      if isreal(R3)
        [R3, tau3, info] = f08bh(R3,'m',ir3,'n',mp0);
        [Z3,info] = f08bk('R','T',r0r2,R3,tau3,P3);
      else
        [R3, tau3, info] = f08bv(R3,'m',ir3,'n',mp0);
        if isreal(P3)
          P3 = complex(P3,zeros(n+r0));
        end
        [Z3,info] = f08bx('R','C',r0r2,R3,tau3,P3);
      end
      R3 = triu(R3(1:r3,1:r3));
    else
      [Z3,R3] = qr(R3');
      Z3 = P3*Z3;
      R3 = R3(1:r3,1:r3)';
    end
  else
    Z3 = eye(n+r0);
  end

  AA(r2+1:n,1:n+r0) = [R3 zeros(n-r2,r0+r2)];
  AA(r2+1:n,n+r0+1:2*n) = Q3'*AA(r2+1:n,n+r0+1:2*n);
  AA(1:r2,1:n+r0) = AA(1:r2,1:n+r0)*Z3;
  AA(n+1:n+r0,1:n+r0) = AA(n+1:n+r0,1:n+r0)*Z3;
  BB(:,1:n+r0) = BB(:,1:n+r0)*Z3;
  p = [(n+1:n+r0) (r2+1:n)];
  AA(r2+1:n+r0,:) = AA(p,:);
  BB(r2+1:n+r0,:) = BB(p,:);
  p = [(n-r2+1:n+r0) (1:n-r2)];
  AA(:,1:n+r0) = AA(:,p);
  BB(:,1:n+r0) = BB(:,p);

  if nargout > 1
    V = [Z3 zeros(n+r0,n-r0); zeros(n-r0,n+r0) eye(n-r0)];
    V(n+1:2*n,:) = Q0*V(n+1:2*n,:);
    ind1 = (1:n+r0);
    ind2 = [(n-r2+1:n+r0) (1:n-r2)];
    V(:,ind1) = V(:,ind2);
  end

  if nargout > 2
    Q = [[Q2(:,1:r2) Q2(:,r2+1:n)*Q3]' zeros(n);
         zeros(n)   Q0'];
    ind1 = (r2+1:n+r0);
    ind2 = [(n+1:n+r0) (r2+1:n)];
    Q(ind1,:) = Q(ind2,:);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigenvalue/eigenvector computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nagtoolbox || nargout < 3
  % Deflated zero and infinite eigenvalues
  E(r0+r2+1:n+r0,1) = Inf*ones(n-r2,1);
  E(n+r0+1:2*n,1) = zeros(n-r0,1);
  iqz = (1:r0+r2);  % Active part of linearization
else
  iqz = 1:2*n;
end

% For the remaining eigenvalues use the QZ algorithm (via eig or qz)
% on active part of the big matrix pair (AA,BB).
if nargout <= 1

  E(iqz) = eig( AA(iqz,iqz), BB(iqz,iqz) );
elseif nargout == 2

  [XX, temp] = eig( AA(iqz,iqz), BB(iqz,iqz) );
  E(iqz) = diag(temp);
elseif eigleft && r0+r2 == 2*n

  [XX, temp, YY] = eig(AA, BB);
  E = diag(temp);

  j = 1;
  jcplx = [];  % contains indices of complex eigenvalues
  while j < length(2*n)
    if imag(E(j)) == 0
      j = j + 1;
    else
      jcplx = [jcplx j];
      j = j + 2;
    end
  end

else
  % Use qz explicitly, since eig does not return left eigenvectors.
  if Qreal

    if nagtoolbox
      % left eigenvectors are computed later by f08yk
      [TA,TB,QQ,ZZ,XX] = qz(AA(iqz,iqz),BB(iqz,iqz),'real');

    else
      % Need to call qz on full linearization so as to get
      % the left eigenvectors.
      [TA,TB,QQ,ZZ,XX,YY] = qz(AA,BB,'real');
    end

    warning('off', 'MATLAB:divideByZero')
    j = 1;
    jcplx = [];  % contains indices of complex eigenvalues
    while j < length(iqz)
      if TA(j+1,j) == 0
        E(j) = TA(j,j)/TB(j,j);
        j = j + 1;
      else
        jcplx = [jcplx j];
        ind  =  [j j+1];
        E(ind) = eig(TA(ind,ind),TB(ind,ind));
        j = j + 2;
      end
    end
    if j == length(iqz)
      E(j) = TA(j,j)/TB(j,j);
    end
    warning('on', 'MATLAB:divideByZero')

  else
    if nagtoolbox
      % left eigenvectors are computed later by f08yx
      [TA,TB,QQ,ZZ,XX] = qz(AA(iqz,iqz),BB(iqz,iqz),'complex');

    else
      % Need to call qz on full linearization so as to get
      % the left eigenvectors.
      [TA,TB,QQ,ZZ,XX,YY] = qz(AA,BB,'complex');
    end
    warning('off', 'MATLAB:divideByZero')
    E(iqz) = diag(TA)./diag(TB);
    warning('on', 'MATLAB:divideByZero')
  end
end

if ~isempty(find(isnan(E),1))
  error('MATLAB:quadeig:SingularQuadratic', ...
        'The quadratic matrix polynomial is nonregular');
end

% Right eigenvectors recovery
if nargout > 1

  X = zeros(n,2*n);

  if ((options.pscale == 1 || options.pscale == 2) && tau < toltau) || r0 < n

    X(:,iqz) = V(1:n,iqz)*XX;

    % Normalize the vectors
    if nargout == 2 || ~Qreal
      for j = iqz
        X(:,j) = X(:,j)/norm(X(:,j));
      end
    else % nargout > 2 && Qreal

      % Note: When called with 'real', qz returns a real matrix X
      % for the eigenvector matrix, so that the jth eigenvector is X(j)
      % when e(j) is real, but when e(j) is complex the corresponding
      % complex conjugate pair is X(j)+/-X(j+1).
      % recall that jcplx contains indices for complex eigenvalues.
      % Form and normalize the vectors and make the element of largest
      % modulus real and positive
    
      for j = jcplx
        X(:,j) = complex(X(:,j),X(:,j+1));
      end
      j = 1;
      % Normalize the vectors
      while j <= length(iqz)
        if imag( E(j) ) == 0
          X(:,j) = X(:,j)/norm(X(:,j));
          j = j + 1;
        else
          X(:,j) = X(:,j)/norm(X(:,j));
          X(:,j+1) = conj(X(:,j));
          j = j + 2;
        end
      end

    end
  else

    XX = V(:,iqz)*XX; % 2n x (r0+r2) matrix.

    if nargout > 2  &&  Qreal
      for j = jcplx
        XX(:,j) = complex(XX(:,j),XX(:,j+1));
        XX(:,j+1) = conj(XX(:,j));
      end
    end

    XX(n+1:2*n,:) = P0*(R0\(Q0'*XX(n+1:2*n,:))); %x2 = A0\z2

    % Normalize the vectors
    j = 1;
    while j <= length(iqz)
      if imag( E(j) ) == 0 || ~Qreal
        XX(1:n,j) = XX(1:n,j)/norm(XX(1:n,j));
        XX(n+1:2*n,j) = XX(n+1:2*n,j)/norm(XX(n+1:2*n,j));
        j = j + 1;
      else
        XX(1:n,j) = XX(1:n,j)/norm(XX(1:n,j));
        XX(1:n,j+1) = conj(XX(1:n,j));
        XX(n+1:2*n,j) = XX(n+1:2*n,j)/norm(XX(n+1:2*n,j));
        XX(n+1:2*n,j+1) = conj(XX(n+1:2*n,j));
        j = j + 2;
      end
    end

    % We now have two approximate right eigenvectors per eigenvalue.
    % Keep the vector with smallest backward error.
    W = zeros(n,2); absW= zeros(1,2);
    for j = iqz
      W(:) = XX(:,j);
      absW = sum(abs(W),1);
      if absW(1) ~= 0 && absW(2) ~= 0
        Z = A;
        if ~isinf(E(j))
          Z = B + E(j)*Z;
          Z = C + E(j)*Z;
        end
        Z = Z*W;
        res = sum(abs(Z))./absW;
        [ignore,ind] = min(res);
        X(1:n,j) = W(:,ind); % Eigenvector with smallest backward error.
      else
        if absW(1) ~= 0
          X(1:n,j) = W(:,1);
        else
          X(1:n,j) = W(:,2);
        end
      end
    end
  end % ((options.pscale == 1 || options.pscale == 2) && tau < toltau) || r0 < n

  % Right eigenvectors corresponding to deflated eigenvalues.
  if r0 < n
    if r0 == 0
      X0 = eye(n);
    elseif nagtoolbox

      if b32
        mm0 = int32(n-r0);
      else
        mm0 = int64(n-r0);
      end

      % Reduce R0 to lower triangular form.
      if isreal(R0)
        [R0, tau0, info] = f08bh(R0);
        [X0,info] = f08bk('R','T',mm0,R0,tau0,P0);
      else
        [R0, tau0, info] = f08bv(R0);
        if isreal(P0)
          P0 = complex(P0,zeros(n));
        end
        [X0,info] = f08bx('R','C',mm0,R0,tau0,P0);
      end

    else
      [X0,R0] = qr(R0');
      X0 = P0*X0;
    end
    X(:,n+r0+1:2*n) = X0(:,r0+1:n);

  end

  if r2 < n
    if r2 == 0
      X2 = eye(n);
    elseif nagtoolbox
      % Reduce R2 to lower triangular form.
      if isreal(R2)
        [R2, tau2, info] = f08bh(R2);
        [X2,info] = f08bk('R','T',mm2,R2,tau2,P2);
      else
        [R2, tau2, info] = f08bv(R2);
        if isreal(P2)
          P2 = complex(P2,zeros(n));
        end
        [X2,info] = f08bx('R','C',mm2,R2,tau2,P2);
      end
    else
      [X2,R2] = qr(R2');
      X2 = P2*X2;
    end
    X(:,r2+r0+1:n+r0) = X2(:,r2+1:n);
  end
  
  % Make the element of largest modulus of each eigenvector real and
  % positive
  if Qreal
    j = 1;
    while j <= 2*n
      if imag( E(j) ) == 0
        [xx,i] = max(abs(X(1:n,j)));
        X(1:n,j) = sign(X(i,j))*X(1:n,j);
        j = j + 1;
      else
        [xx,i] = max(abs(X(1:n,j)));
        X(1:n,j) = X(1:n,j)/(X(i,j)/xx);
        X(i,j) = real(X(i,j));
        X(1:n,j+1) = conj(X(1:n,j));
        j = j + 2;
      end
    end
  else
    for j = 1: 2*n
      [xx,i] = max(abs(X(1:n,j)));
      X(1:n,j) = X(1:n,j)/(X(i,j)/xx);
      X(i,j) = real(X(i,j));
    end
  end

end

% Left eigenvectors recovery.
if nargout > 2
  if nagtoolbox
    % Need to compute left eigenvectors of linearization.
    % Transform 2nx2n linearization to Schur form.
    AA(iqz,iqz) = TA; BB(iqz,iqz) = TB;
    AA(iqz,r0+r2+1:2*n) = QQ*AA(iqz,r0+r2+1:2*n);
    BB(iqz,r0+r2+1:2*n) = QQ*BB(iqz,r0+r2+1:2*n);
    Q(iqz,:) = QQ*Q(iqz,:);

    select = logical(ones(2*n,1));
    for j = r0+r2+1:2*n
      select(j) = false;
    end
    YY = zeros(2*n,r0+r2); temp = zeros(2*n,1);

    if b32
      r0r2 = int32(r0+r2);
    else
      r0r2 = int64(r0+r2);
    end

    if Qreal
      [YY,temp,m,info] = f08yk('L', 'S', select,AA,BB, YY, temp, r0r2);
    else
      YY = complex(YY,YY);
      temp = complex(temp,temp);
      [YY,temp,m,info] = f08yx('L', 'S', select,AA,BB, YY, temp, r0r2);
    end
  end
  YY = Q'*YY;

  % In the real case, form the actual eigenvectors from Y.
  if Qreal
    for j = jcplx
      YY(:,j) = complex(YY(:,j),YY(:,j+1));
      YY(:,j+1) = conj(YY(:,j));
    end
  end

  Y = zeros(n,2*n);

  if ((options.pscale == 1 || options.pscale == 2) && tau < toltau)

    Y(:,iqz) = YY(1:n,:);
    ind = find( abs(E(iqz))<1 );
    Y(:,ind) = YY(n+1:2*n,ind);
    % Normalize the vectors
    for j = iqz
      Y(:,j) = Y(:,j)/norm(Y(:,j));
    end

  else

    % We have two approximate left eigenvectors per eigenvalue.
    % Keep the vector with smallest backward error.
    W = zeros(n,2);
    for j = iqz
      W(:) = YY(:,j);
      absW = sum(abs(W),1);
      if absW(1) ~= 0 && absW(2) ~= 0
        Z = A;
        if ~isinf(E(j))
          Z = B + E(j)*Z;
          Z = C + E(j)*Z;
        end
        Z = W'*Z;
        res = sum(abs(Z),2)./absW'; % normalized residuals.
        [ignore,ind] = min(res);
        Y(1:n,j) = W(:,ind)/norm(W(:,ind)); % Eigenvector with unit 2-norm.
      else
        if absW(1) ~= 0
          Y(1:n,j) = W(:,1)/norm(W(:,1));
        else
          Y(1:n,j) = W(:,2)/norm(W(:,2));
        end
      end
    end
  end % if (options.pscale ~= 1 && tau < toltau)

  % Compute left eigenvectors of deflated zero and infinite eigenvalues.
  if r0 < n
    Y(:,n+r0+1:2*n) = Q0(:,r0+1:n);
  end
  if r2 < n
    Y(:,r2+r0+1:n+r0) = Q2(:,r2+1:n);
  end

  % Make the element of largest modulus of each eigenvector real and positive
  if Qreal
    j = 1;
    while j <= 2*n
      if imag( E(j) ) == 0
        [yy,i] = max(abs(Y(1:n,j)));
        Y(1:n,j) = sign(Y(i,j))*Y(1:n,j);
        j = j + 1;
      else
        [yy,i] = max(abs(Y(1:n,j)));
        Y(1:n,j) = Y(1:n,j)/(Y(i,j)/yy);
        Y(i,j) = real(Y(i,j));
        Y(1:n,j+1) = conj(Y(1:n,j));
        j = j + 2;
      end
    end
  else
    for j = 1: 2*n
      [yy,i] = max(abs(Y(1:n,j)));
      Y(1:n,j) = Y(1:n,j)/(Y(i,j)/yy);
      Y(i,j) = real(Y(i,j));
    end
  end

end

if rev
  warning('off', 'MATLAB:divideByZero')
  E = 1./E;
  warning('on', 'MATLAB:divideByZero');
  temp = A; A = C; C = temp;
end

% If scaling took place, scale eigenvalues back.
if options.pscale
  E = E(:);
  E = gg*E;
end

if nargout < 2
  X = E;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigenvalue condition numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 3
  nA = g2; nB = g1; nC = g0;
  if options.pscale
    A = (1/(gg^2*delta))*A;
    B = (1/(gg*delta))*B;
    C = (1/delta)*C;
  end

  % Reconstruct alpha-beta representation of eigenvalues:
  %    E(i) = alpha(i)/beta(i).
  a = E(:);
  b = ones(size(a));
  k = isinf(a); a(k) = 1; b(k) = 0;
  % Normalize the eigenvalues (a,b) so that they lie on the unit circle
  nab = sqrt(abs(a).^2 + abs(b).^2);
  a = a./nab;
  b = b./nab;
  ab = [abs(a).^2 abs(a.*b) abs(b).^2];

  % Compute condition numbers.
  abnu = [ab(:,1)*nA  ab(:,2)*nB ab(:,3)*nC];
  s = zeros(2*n,1);
  DAX = 2*A*X*diag(a.*conj(b)) + B*X*diag(abs(b).^2-abs(a).^2) - ...
        2*C*X*diag(conj(a).*b);
  for j = 1:2*n
    as =Y(:,j)'*DAX(:,j);
    if as ~= 0
      s(j) = norm(abnu(j,:))/abs(as);
    else
      s(j) = Inf;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Backward error of right/left eigenpairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout > 4
  be_X = zeros(2*n,1);
  p = nA*abs(a.^2)+nB*abs(a.*b)+nC*abs(b.^2);
  tmp = options.tol*max([nA nB nC]);
  if nA <= tmp || nC <= tmp
    ind = find(p<=tmp);
    p(ind) = 1;
  end
  Res = A*X*diag(a.^2)+B*X*diag(a.*b)+C*X*diag(b.^2);
  for j = 1:2*n
    be_X(j) = norm(Res(:,j))/p(j);
  end
  be_X = be_X(:);
end

if nargout > 5
  be_Y = zeros(2*n,1);
  % Compute backward errors of left eigenpairs
  Res = diag(a.^2)*Y'*A+diag(a.*b)*Y'*B+diag(b.^2)*Y'*C;
  for j = 1:2*n
    be_Y(j) = norm(Res(j,:))/p(j);
  end
  be_Y = be_Y(:);
end

if g0*g2 == 0
  tau = -1;
end
