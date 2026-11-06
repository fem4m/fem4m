function fem4matLogo
%fem4matLogo - Render and save the fem4mat logo
%   fem4matLogo computes the first Dirichlet eigenfunction of the
%   Laplacian on the L-shaped domain with P1 elements (Pk) on the
%   eighth uniform refinement of the initial mesh, plots its absolute
%   value from above in figure 1 with a flat colormap of three colors
%   (blue, sky blue, orange), whose bands follow the level lines, and
%   exports the figure to fem4mat.png in the current folder with 75 dpi
%   and a transparent background.
%
%   See also problemScalarEllipticEigenvalueLshape, Pk, eigs,
%   exportgraphics

% Copyright (c) 2026 by the fem4mat authors

p = problemScalarEllipticEigenvalueLshape;

mesh = Mesh(p.geom);
mesh = mesh.refine('uniform',8);

pk = Pk(1,mesh);
a = BilinearForm('grad(u)*grad(v)dx');
A = a.assemble(pk);
b = BilinearForm('u*v dx');
B = b.assemble(pk);
x = zeros(size(A,1),1);
freeNodes = pk.freeNodes;
opts = [];
opts.disp = 0;
[V,D] = eigs(A(freeNodes,freeNodes),B(freeNodes,freeNodes),1,'sm',opts);
[~,ind] = sort(diag(D));
x(freeNodes) = V(:,ind(1));

figure(1);clf;
plot(pk,abs(x));shading interp;
view(2);
axis equal;
axis off;
colormap([0 76 153; 0 158 219; 255 153 0]/255);
drawnow;

exportgraphics(gcf,'fem4mat.png',...
    'ContentType','image',...
    'Resolution',75,...
    'BackgroundColor','white');
% images get no transparent background from exportgraphics: make the
% white one (no channel below 200, unlike the three colors) transparent
img = imread('fem4mat.png');
imwrite(img,'fem4mat.png','Alpha',double(min(img,[],3) <= 200));


end