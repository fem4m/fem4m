# forms

## Purpose
This folder provides the classes that parse variational forms from strings and assemble them. `Form` is the abstract base class that parses a form string such as `"grad(u)*grad(v) dx"` into terms with operators, integration domains and scaling powers. `BilinearForm` assembles matrices, `LinearForm` load vectors, `ErrorForm` global error values and `EstimatorForm` local error indicators (one value per element). `FormBasis` provides the operators on basis functions (used by `BilinearForm` and `LinearForm`), `FormValue` those on finite element functions (used by `EstimatorForm` and `ErrorForm`) and `FormEval` the evaluation of coefficients and functions during assembly. `FormType` names the value types these operate on (`scalar`, `vector2d`, `vector3d`, `matrix2d`, `matrix3d`, `symmatrix2d`, `symmatrix3d`, `tensor`).

## Form strings
A form string is a sum of terms `[sign] [h^k*] [p^k*] [coefficient(x)*] integrand domain`, for example

```matlab
A = BilinearForm('grad(u)*grad(v)dx - [u]*n.{grad(v)}ds - {grad(u)}.[v]*nds + h^-1*alpha(x)*[u][v] ds');
A = A.assemble(fe,'alpha',penalty);
b = LinearForm('f(x) v dx + g(x) v dsn');
b = b.assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
eta = EstimatorForm('h^2*(f(x)+div(grad(u)))^2 dx');
err = ErrorForm('(exactGrad(x)-grad(u))^2 dx');
```

- Variables: the trial variable `u` and the test variable `v` by default; other names are set with a second argument, e.g. `BilinearForm('sigma*tau dx',{'sigma','tau'})`; `BilinearForm` and `LinearForm` need a cell `{TRIAL,TEST}` of two different names, `EstimatorForm` and `ErrorForm` also take a single name for `u`, and their second name (`v` by default) refers to a second discrete function passed to `assemble`.
- Operators on a variable: `grad`, `div`, `curl`, `rot`, `eps` (symmetric gradient), `sym`, `tr`, `d2` (Hessian) and `laplace`, in `EstimatorForm` and `ErrorForm` also `skew` and `biharmonic` (`dev` is parsed in `BilinearForm` and `LinearForm`, but no element provides it yet), and their compositions `divdiv(sigma)` (double divergence of a symmetric matrix field) and `div(d2(u))` = grad(laplace(u)) (scalar fields); a trailing `n` or `t` takes the normal or tangential component (e.g. `vn`, `grad(v)n`, `tau.n`), `xn` the product with the normal: on edges the jump `[v]` times the normal (the outer product for vectors), on the faces of `Mesh3d` the cross product (`grad(v)xn` and `[grad(v)]xn` only on faces and only in `BilinearForm` and `LinearForm`).
- Curl and rot follow the notation of Girault and Raviart. In 2d, `curl(u)` = (d_y u, -d_x u) is the vector curl of a scalar field and `rot(v)` = d_x v_2 - d_y v_1 the rot of a vector field, row-wise for matrix fields, so that (curl u, v) = (u, rot v) for u vanishing on the boundary. The second-order operators are `curl(curl(u))` (the Airy operator of a scalar field), `curl(rot(v))` (vector fields) and `rot(rot(sigma))` (symmetric matrix fields), also written `curlcurl`, `curlrot` and `rotrot`. `sym(curl(v))` is the symmetric part of the row-wise curl, with the rows (d_y v_i, -d_x v_i), of a 2d vector field, the only place where curl acts on a 2d vector field. In 3d, `curl(u)` and `curl(curl(u))` apply to vector fields; `rot` is a 2d operator. An element that does not provide one of these curl and rot operators raises an error that names the operators of the dimension.
- Jumps and averages over the sides: `[u]` and `{u}`, e.g. `[u][v]` or `[sigma]*n {mu}` in bilinear forms and `g(x)*[v] ds` in linear forms. The jump is the value on the first neighboring element minus that on the second, with the normal of the first element; on boundary sides jump and average are the one-sided value.
- Domains: `dx` integrates over the elements and `ds` over all sides, i.e. the edges in 2d and the faces in 3d; the suffixes `i`, `d`, `n` and `b` restrict to interior, Dirichlet, Neumann and all boundary sides (`dsb` only in `BilinearForm` and `ErrorForm`), `dss` integrates over the symmetry sides (edges in 2d, faces in 3d).
- Sides: the side operators (suffix `_S`) and the `ds` domains act on the edges of 2d, surface and polygonal meshes and on the faces of `Mesh3d`, decided by the mesh class (`Form.sideOf`), not by the number of coordinates. Inside the forms, `Form.sideOf(mesh,'s4e','e4s',...)` returns the edge or face data under one name (`n4s`, `e4s`, `s4e`, `nrSides`, `s4Db`/`s4Nb`/`s4Sb` and `normal4s` for `n4ed`/`n4f` and so on, `diam4s` for the edge length or the longest edge of the face), and `Form.sideOf(fe,'dof4s','degree4s')` the side dofs and degrees of an element; the meshes and elements keep only their edge and face names.
- Scaling: for an integer k, `h^k*` multiplies a term by the k-th power of the local mesh size, `p^k*` by the k-th power of the local polynomial degree.
- Coefficients and functions: names written as `name(x)` are passed to `assemble` as name-value pairs, either as function handles or as values, e.g. `'f',@p.f`; a function (not a coefficient) given as a handle of `(x,n)` receives the side normals in side terms. Names without a value are allowed: a coefficient without a value is the identity, a function without a value is zero, and a weight `W(x)` of `EstimatorForm` or `ErrorForm` without a value leaves the term unweighted. A coefficient of a vector or matrix operator is a scalar or the matrix acting on its components: d-by-d for vectors, 4-by-4 for 2d matrices and 3-by-3 or 6-by-6 for symmetric matrices with components (11,22,12) or (11,22,33,12,13,23); in `b*grad(u)*v` a vector.
- Subdomains: `assemble` of `BilinearForm` and `LinearForm` takes the indices of the elements (`dx` terms) or sides (`ds` terms) to assemble on, e.g. `BilinearForm('u*v ds').assemble(fe,fe,[],mesh.ed4Db)`.
- Quadrature: the optional `DEGREE` argument of `assemble` is the polynomial degree that the quadrature integrates exactly, by default the sum of the trial and test degrees in `BilinearForm` and twice the element degree in the other forms (plus the mesh degree minus one on curved meshes in `BilinearForm` and `LinearForm`, twice the larger degree with a second function in `EstimatorForm` and `ErrorForm`).
- Powers: for terms `(...)^p`, `EstimatorForm` and `ErrorForm` return the p-th root of the sum of the terms, i.e. the indicator or the error norm itself; a side term of `EstimatorForm` counts one half for each of the two neighboring elements of an interior side.

## Relation to other parts
Assembly needs an element object (from `elements`) on a mesh (from `mesh`); the coefficient and data functions usually come from a problem class (from `problems`).

## Contents

{{MEMBERS}}

## Related examples

{{RELATED_EXAMPLES}}
