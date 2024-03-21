# Getting Started

## Installation
Add the folder `toolbox` of the repository and all its subfolders to the MATLAB path. If you saved fem4m in the folder `~/fem4m`, type

```matlab
addpath(genpath('~/fem4m/toolbox'));
```

Alternatively, change to the folder `toolbox/internal` of the repository and call `fem4mPath.addpath()`, which adds all folders of the repository to the path except those whose name starts with `.`, such as `.git`:

```matlab
cd ~/fem4m/toolbox/internal
fem4mPath.addpath()
```

In the Command Window, `help <foldername>` lists the contents of a folder, e.g. `help elements`, and `help <name>` shows the help text of a class or function.

## Requirements
- MATLAB R2026a; the toolbox uses only core MATLAB.
- The examples are plain-text live scripts, which MATLAB opens as live scripts with formatted text.

## Minimal workflow
The simplest example is `exPoissonP1`, a P1 finite element method for the Poisson problem on the unit square:

```matlab
p = problemScalarEllipticSquareExact();        % problem data and exact solution
mesh = Mesh(p.geom);                           % coarse mesh of the unit square
mesh = mesh.refine('uniform');
for lvl = 1:6
    fe = Pk(1,mesh);                           % P1 finite element space
    a = BilinearForm('grad(u) * grad(v) dx');
    A = a.assemble(fe);                        % stiffness matrix
    l = LinearForm('f(x) v dx');
    b = l.assemble(fe,'f',@p.f);               % load vector
    x = zeros(fe.nrDoF,1);
    free = fe.freeNodes;                       % u = 0 on the boundary
    x(free) = A(free,free)\b(free);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    disp(h1.assemble(x,fe,'exact',@p.exactGradU))  % error in the energy norm
    mesh = mesh.refine('uniform');
end
```

The page of the example `exPoissonP1` derives every step from the Poisson problem and shows the output with the convergence rates. The [Examples](examples/examples_overview.html) page lists all examples by topic and suggests an order for working through them.
