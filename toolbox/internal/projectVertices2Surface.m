function val = projectVertices2Surface(v0,surfaceF,gradSurfaceF,tol,maxIter,display)
%projectVertices2Surface - Project vertices onto implicit surface
%   VAL = projectVertices2Surface(V0,SURFACEF,GRADSURFACEF) returns the
%   normal projections VAL (N-by-3) of the points V0 (N-by-3, one row
%   per point) onto the zero level set of the function SURFACEF with
%   the gradient GRADSURFACEF: V0 = VAL + D*n(VAL) with n the
%   normalized gradient and D the distance, signed as SURFACEF(V0)
%   (a point of V0 on the surface stays where it is). SURFACEF(X)
%   returns an N-by-1 vector and GRADSURFACEF(X) an N-by-3 array for
%   the points X (N-by-3). Each iteration takes a Newton step for
%   SURFACEF and moves the point back onto the line through V0 along
%   the normal at the new point, keeping the distance to V0. The
%   iteration stops when the distance estimate
%   |SURFACEF(VAL)|/|GRADSURFACEF(VAL)| is at most 1e-12 for all points
%   or after 1000 iterations.
%
%   VAL = projectVertices2Surface(...,TOL) uses the tolerance TOL for
%   the distance estimate (default 1e-12).
%
%   VAL = projectVertices2Surface(...,TOL,MAXITER) stops after at most
%   MAXITER iterations (default 1000).
%
%   VAL = projectVertices2Surface(...,TOL,MAXITER,DISPLAY) prints the
%   number of iterations, the largest distance estimate and the number
%   of points when DISPLAY is true (default false).
%
%   An empty TOL, MAXITER or DISPLAY selects the default.
%
%   See also SurfaceTransformation, MeshSurface, ProblemSurfaceTemplate

% Copyright (c) 2026 by the fem4m authors

if nargin < 4 || isempty(tol)
    tol = 1e-12;
end
if nargin < 5 || isempty(maxIter)
    maxIter = 1e3;
end
if nargin < 6 || isempty(display)
    display = false;
end

fVal = surfaceF(v0);
gradFVal = gradSurfaceF(v0);
signV0 = sign(fVal);
val = v0;
normGradSquared = sum(gradFVal.^2,2);
iter = 0;
error = 1;
while max(error) > tol && iter < maxIter
    val = val - fVal(:,[1 1 1]).*gradFVal./normGradSquared(:,[1 1 1]);
    dist = signV0.*sqrt(sum((val-v0).^2,2));
    gradFVal = gradSurfaceF(val);
    normGrad = sqrt(sum(gradFVal.^2,2));
    val = v0 - dist(:,[1 1 1]).*gradFVal./normGrad(:,[1 1 1]);
    %estimate error
    fVal = surfaceF(val);
    gradFVal = gradSurfaceF(val);
    normGradSquared = sum(gradFVal.^2,2);
    error = sqrt(fVal.^2./normGradSquared);
    iter = iter + 1;
end

if display
    disp(['iter = ',num2str(iter),' maxError = ',num2str(max(error)),' nrNodes = ',num2str(size(v0,1))]);
end

end