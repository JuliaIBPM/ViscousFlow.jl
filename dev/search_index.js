var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#ViscousFlow-1",
    "page": "Home",
    "title": "ViscousFlow",
    "category": "section",
    "text": "a framework for simulating viscous incompressible flowsThe objective of this package is to allow easy setup and fast simulation of incompressible flows, particularly those past bodies in motion. The package provides tools forconstructing grids and body shapes,\nusing the operators on those grids,\nspecifying the relevant parameters and setting their values,\nsolving the problem, and finally,\nvisualizing and analyzing the results.The underlying grids are uniform and Cartesian, allowing the use of the lattice Green\'s function (LGF) for inverting the Poisson equation; the diffusion operators are solved with the integrating factor (Liska and Colonius ref). Many of the core aspects of the fluid-body interaction are based on the immersed boundary projection method, developed by Taira and Colonius (ref). The coupled fluid-body interactions are based on the work of Wang and Eldredge (ref)."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "This package works on Julia 0.6, 0.7 and 1.0. To install, simply runjulia> using Pkg\n\njulia> Pkg.add(\"https://github.com/jdeldre/ViscousFlow.jl\")in the Julia REPL. Since this package is still under heavy development, you should runjulia> Pkg.test(\"ViscousFlow\") # might take some timeto make sure things are working as intended andjulia> Pkg.update()to get the most recent version of the library and its dependencies.The plots in this documentation are generated using Plots.jl. You might want to install that too to follow the examples."
},

{
    "location": "#Basic-Usage-1",
    "page": "Home",
    "title": "Basic Usage",
    "category": "section",
    "text": "Do something here."
},

{
    "location": "manual/fields/#",
    "page": "Fields",
    "title": "Fields",
    "category": "page",
    "text": ""
},

{
    "location": "manual/fields/#Fields-1",
    "page": "Fields",
    "title": "Fields",
    "category": "section",
    "text": "DocTestSetup = quote\n  using ViscousFlow\n  using Random\n  Random.seed!(1)\nenddefddt1fracmathrmd1mathrmdt\n\nrenewcommandvecboldsymbol\nnewcommanduvec1vechat1\nnewcommandutangentuvectau\nnewcommandunormaluvecn\n\nrenewcommanddmathrmdusing ViscousFlow\nusing PlotsIn whirl, field data, such as velocity, vorticity and pressure, are stored on a staggered uniform grid. Such a grid is divided into cells, with edges (which, on a two-dimensional grid, are the same as faces) and nodes (cell centers). Nodes hold scalar-valued data. Edges, on the other hand, hold the components of vector-valued data that are normal to the respective edges; one component lies on the vertical edges, while the other is on the horizontal edges.Furthermore, there are two different cell types: primal and dual. On the physical grid, these cell types are offset with respect to each other by half a cell spacing in each direction. In other words, the four corners of the primal (resp. dual) cell are the nodes of four dual (resp. primary) cells.Thus, on a two-dimensional staggered grid, there are four distinct vector spaces, associated with where the data are held on the grid:dual nodes,\ndual edges,\nprimal nodes, and\nprimal edges.In whirl, these are each distinct data types. Furthermore, the relationships between these types are defined by an underlying grid shared by all. By convention, this grid is defined by the number of dual cells NX and NY in each direction; we will often refer to it as the dual grid. For example, Nodes{Dual,NX,NY} is the type for dual node data on this grid; Edges{Primal,NX,NY} is the type for edge data on the primal cells within this same NX by NY dual grid. Note that, even though this latter type is parameterized by NX and NY, these values do not correspond to the number of primal edges in each direction on this dual grid. These values always correspond to the number of dual cells on the grid, for any data type. This makes it clear the grid is shared by all data."
},

{
    "location": "manual/fields/#Setting-up-field-data-1",
    "page": "Fields",
    "title": "Setting up field data",
    "category": "section",
    "text": "Let\'s see an example of creating a blank set of dual node data and filling it with something:w = Nodes(Dual,(5,4))\nw .= reshape(1:20,5,4)Other data types on the same grid can be set up in similar fashion. To ensure that they have a size that is consistent with the dual node data w, we can use this in place of the size:q = Edges(Primal,w);\nq.u[2,3] = 1;\nq"
},

{
    "location": "manual/fields/#Field-differencing-operations-1",
    "page": "Fields",
    "title": "Field differencing operations",
    "category": "section",
    "text": "Field operations transform one data type to another. Some of these are differencing operations, analogous to differential counterparts in continuum calculus: curl, divergence, and gradient. For example, a curl operation can act upon dual nodal data (like streamfunction) and return primal edge data (i.e. velocity); a divergence operation acts on edge data (primal or dual) and returns nodal data of the same cell type. Note that these operations are mimetic: they maintain some of the same properties as the continuous counterparts. For example, the divergence of the curl of any dual nodal data is exactly zero. The curl of the gradient of primal nodal data is also zero.Let\'s take the curl of the dual nodal data we constructed:curl(w)We could also make this a little more cute by giving the curl operator a symbol and then acting upon the data as though it were a matrix-vector operation:C = Curl()\nC*wNote that C is not actually a matrix. Rather, it is simply another name for the curl operator, and * is defined in this context to apply curl to whatever is to the right of it. The other operators have similar constructs.Suppose we wish to apply the curl operation over and over. The curl() function allocates memory for the result whenever it is used; this would become expensive if it is done often. So it makes sense to preallocate space for this result and use the curl!() function, which simply fills in the elements:q = Edges(Primal,w)\ncurl!(q,w)Note that we used a convenience function for setting up primal edge data q of a size that corresponds with w.Let\'s check that divergence of the curl is indeed zero:D = Divergence()\nD*(C*w)"
},

{
    "location": "manual/fields/#The-Laplacian-and-its-inverse-1",
    "page": "Fields",
    "title": "The Laplacian and its inverse",
    "category": "section",
    "text": "ViscousFlow also makes heavy use of the discrete Laplacian operator, L. This mimics the continuous operator, nabla^2, and acts upon data of any type. Let\'s apply this to the original data:laplacian(w)As with the other operators, we can also construct a shorthand of the discrete Laplacian operator,L = plan_laplacian(size(w))\nL*wAn important part of whirl is the inverse of this operator. That is, we need the ability to solve the discrete Poisson systemLs = wfor s, for given data w. We achieve this in whirl with the lattice Green\'s function. To outfit the operator with its inverse, we simply set the optional flag:L = plan_laplacian(size(w),with_inverse=true)Then, the Poisson system is solved with the backslash (\\),s = L\\w\nL*sIt should be observed that the cells on the perimeter have not recovered the original values of w. These are the ghost cells, and the Laplacian operation does not apply to these.It is also important to note that, although it looks as though we\'ve constructed a matrix L and performed various matrix-vector operations with it, this is not actually the case. In fact, the \\ operation associated with L is significantly faster than a matrix inversion. Internally, it carries out a fast convolution between the data in w and the lattice Green\'s function, via fast Fourier transform. The lattice Green\'s function (LGF) table is pre-computed and pre-transformed in the original construction of L. (In fact, because this table is not dependent on the size of the grid, it is actually computed once for all time and stored in a file; subsequent applications of it just load it in and use the portion of it necessary for a certain grid.)The lattice Green\'s function has the advantage that it is independent of the grid size. Let\'s solve the Poisson system when w is a unit field, i.e. a field of zeros, except for a single 1 entry at one node. The solution s represents the influence of this point on all nodes. To see that the LGF does not depend on the grid size, let\'s use a grid that is long and skinny and plot the solution on itw = Nodes(Dual,(50,10));\nw[20,5] = 1.0\nL = plan_laplacian(w,with_inverse=true)\nplot(L\\w)\nsavefig(\"Linvw.svg\"); nothing # hide(Image: )The influence is not affected by the narrow grid dimensions."
},

{
    "location": "manual/fields/#The-integrating-factor-1",
    "page": "Fields",
    "title": "The integrating factor",
    "category": "section",
    "text": "An operator related to the lattice Green\'s function is the integrating factor. Suppose we have the system of ODEsddt u = L u + f(ut) quad u(0) = u_0where L is the discrete Laplacian (on an infinite uniform grid), and u are nodal data (and f is a nodal-valued function acting on this nodal data). The exact solution of this problem isu(t) = E(t)u_0 + int_0^t E(t-tau) f(u(tau)tau)mathrmdtauwhere E(t) is the integrating factor (or matrix exponential) for the system. The easiest way to understand the role of E(t) is to consider its behavior when f is zero and u_0 contains a field of zeros except for a single 1 entry at one cell. Let\'s set up this initial data:u0 = Nodes(Dual,(100,100));\nu0[40,50] = 1.0\nplot(u0)\nsavefig(\"w1.svg\"); nothing # hide(Image: )Then, E(t)u_0 diffuses this initial unit perturbation in each direction. Here, we apply it with t = 5:E = plan_intfact(5,u0)\nplot(E*u0)\nsavefig(\"Ew1.svg\"); nothing # hide(Image: )Note that E(0) = I, where I is the identity. Also, the integrating factor has the useful property that E(t+tau) = E(t)E(tau). From these properties, it follows that E^-1(t) = E(-t). Let\'s suppose we wish to advance u from time t = tau-h to time t = tau. Then we can define an auxiliary quantity, v(ttau) = E(tau-t)u(t), and this new quantity satisfies the modified set of ODEsddt v = E(tau-t) fleft E(t-tau) v(ttau)trightquad v(tau-htau) = E(h)u(tau-h)The result of integrating this set of ODEs to t = tau is v(tautau) = u(tau). In other words, the integrating factor allows us to solve a somewhat reduced set of ODEs."
},

{
    "location": "manual/fields/#Other-field-operations-1",
    "page": "Fields",
    "title": "Other field operations",
    "category": "section",
    "text": "Other field operations shift the data, by local averaging, from one data type to another. These operations are all called cellshift!, and they require that the target data be preallocated. For example, to shift dual node data to the dual edges,w = Nodes(Dual,(5,4));\nw .= reshape(1:20,5,4)\nWw = Edges(Dual,w);\ncellshift!(Ww,w)Note that the edges in the ghost cells are 0; these edges are not assigned any values in the shift operation.We can then shift this to primal edges:q = Edges(Primal,w);\ncellshift!(q,Ww)We can also compute the Hadamard (i.e. element by element) product of any data of the same type, e.g.,q∘q"
},

{
    "location": "manual/fields/#The-grid-in-physical-space-1",
    "page": "Fields",
    "title": "The grid in physical space",
    "category": "section",
    "text": "Thus far, we have not had to consider the relationship between the grid\'s index space and some physical space. All of the operations thus far have acted on the entries in the discrete fields, based only on their relative indices, and not on their physical coordinates. In this section, we will discuss the relationship between the grid\'s index space and physical space, and then in the next section we\'ll discuss how we can transfer data between these spaces.Generically, we can write the relationship between the physical coordinates x and y, and the indices i and j of any grid point asx(i) = (i - Delta i - i_0)Delta x quad y(j) = (j - Delta j - j_0)Delta xThe scaling between these spaces is controlled by Delta x, which represents the uniform size of each grid cell; note that grid cells are presumed to be square in whirl. The indices I_0 = (i_0j_0) represent the location of the origin in the index space for primal nodes. Why primal nodes? Since the underlying grid is composed of dual cells, then primal nodes sit at the corners of the domain, so it is the most convenient for anchoring the grid to a specific point. But, since some field data of the same index are shifted by half a cell in one or both directions, then Delta i and Delta j are included for such purposes; these are either 0 or 12, depending on the field type. For example, for a primal node, Delta i = 0, so that x(i_0) = 0; for a dual node, Delta i = 12, so that x(i_0) = -Delta x2.In particular, for our four different data types and their componentsPrimal nodes: Delta i = 0, Delta j = 0\nDual nodes: Delta i = 12, Delta j = 12\nPrimal edges u: Delta i = 12, Delta j = 0\nPrimal edges v: Delta i = 0, Delta j = 12\nDual edges u: Delta i = 0, Delta j = 12\nDual edges v: Delta i = 12, Delta j = 0"
},

{
    "location": "manual/fields/#Regularization-and-interpolation-1",
    "page": "Fields",
    "title": "Regularization and interpolation",
    "category": "section",
    "text": "Based on this relationship between the physical space and the index space, we can now construct a means of transferring data between a point (xy) in the physical space and the grid points in its immediate vicinity. We say that such a point is immersed in the grid. The process of transferring from the point to the  grid is called regularization, since we are effectively smearing this data over some extended neighborhood; the opposite operation, transferring grid field data to an arbitrary point, is  interpolation. In whirl, both operations are carried out with the discrete  delta function (DDF), which is a discrete analog of the Dirac delta function. The  DDF generally has compact support, so that  it only interacts with a small number of grid points in the vicinity of a  given physical location. Since each of the different field types reside at  slightly different locations, the range of indices invoked in this interaction  will be different for each field type.Regularization can actually take different forms. It can be a simple point-wise interpolation, the discrete analog of simply multiplying by the Dirac delta function:f_i delta(mathbfx - mathbfx_i)to immerse a value f_i based at point mathbfx_i = (x_iy_i).Alternatively, regularization can be carried out over a curve mathbfX(s), the analog ofint f(s) delta(mathbfx - mathbfX(s))mathrmdsor it can be performed volumetrically, corresponding toint f(mathbfy) delta(mathbfx - mathbfy)mathrmdmathbfyIn this case, the function f is distributed over some region of space. In each of these cases, the discrete version is simply a sum over data at a finite number of discrete points, and the type of regularization is specified by providing an optional argument specifying the arclength, area or volume associated with each discrete point. These arguments are used to weight the sum.Let\'s see the regularization and interpolation in action. We will set up a ring  of 100 points on a circle of radius 14 centered at (1212). This curve-  type regularization will be weighted by the arclength, ds, associated with each  of the 100 points.  On these points, we will  set vector-valued data in which the x component is uniformly equal to 1.0,  while the y component is set equal to the vertical position relative to the  circle center. We will regularize these vector data to a primal  edge field on the grid in which these points are immersed.using ViscousFlow\nusing Plots\npyplot()n = 100;\nθ = range(0,stop=2π,length=n+1);\nx = 0.5 .+ 0.25*cos.(θ[1:n]);\ny = 0.5 .+ 0.25*sin.(θ[1:n]);\nds = 2π/n*0.25;\nX = VectorData(x,y);The variable X now holds the coordinates of the immersed points. Now we will set up the vector-valued data on these pointsf = VectorData(X);\nfill!(f.u,1.0);\nf.v .= X.v.-0.5;Note that we have ensured that f has the correct dimensions by supplying the coordinate data X. This first step also initializes the data to zeros.Now, let\'s set up the grid. The physical domain will be of size 10 times 10, and we will use 100 dual grid cells in each direction. Allowing a single layer of ghost cells surrounding the domain, we use 102 cells, and set the cell size to 0.01. Also, we will set the (xy) origin to coincide with the lower left corner of the domain.nx = 102; ny = 102;\nq = Edges(Primal,(nx,ny));\nLx = 1.0;\ndx = Lx/(nx-2)Now we set up the regularization operator. To set it up, it needs to know the coordinate data of the set of immersed points, the grid cell size, and the weight to apply to each immersed point. Since this is a regularization of a curve, this weight is the differential arc length ds associated with each point. (This last argument is supplied as a scalar, since it is uniform.)H = Regularize(X,dx,weights=ds)We have omitted some optional arguments. For example, it chooses a default DDF kernel (the Roma kernel); this can be changed with the ddftype argument. Also, the lower left corner, where we\'ve set the origin, is the location of the (11) primal node; this is the default choice for I0 (the tuple I_0 of coordinates in index space discussed in the previous section).Now we can apply the regularization operator. We supply the target field q as the first argument and the source data f as the second argument.H(q,f);\nplot(q)\nsavefig(\"regq.svg\"); nothing # hide(Image: )We could also regularize this to a field of dual edges.p = Edges(Dual,(nx,ny));\nH(p,f);\nplot(p)\nsavefig(\"regp.svg\"); nothing # hide(Image: )Scalar-valued data on the immersed points can only be regularized to nodal fields; the syntax is similar, and the regularization operator does not need to be reconstructed:g = ScalarData(X);\nfill!(g,1.0);\nw = Nodes(Dual,(nx,ny));\nH(w,g);\nplot(w)\nsavefig(\"regw.svg\"); nothing # hide(Image: )For a given regularization operator, H, there is a companion interpolation operator, E. In whirl, this interpolation is also carried out with the same constructed operator, but with the arguments reversed: the grid field data are the source and the immersed points are the target. Note that interpolation is always a volumetric operation, so the weights assigned during the construction of the operator are not used in interpolation. Let\'s interpolate our regularized field back onto the immersed points.f2 = VectorData(X);\nH(f2,q);\nplot(f2.u,lab=\"u\")\nplot!(f2.v,lab=\"v\")\nsavefig(\"interpf.svg\"); nothing # hide(Image: )Note that interpolation is not the inverse of regularization; we don\'t recover the original data when we regularize and then interpolate. However, there is generally a way to scale the quantities on the immersed points and on the grid so that H = E^T. If we want to force these operations to be transposes of each other, we can supply the issymmetric=true flag. But here, we will exclude it so that it defaults to the asymmetric form.H = Regularize(X,dx)This flag will override any supplied weights.If we expect to carry out the regularization and interpolation a lot, then it is often sensible to construct matrix versions of these operators. This construction is sometimes a bit slow, but the resulting operators perform their operations much faster than the matrix-free operators described above. To generate these matrix operators, we have to supply the data types of the source and target of the operation. For example, for regularization from scalar field data to dual node data,g = ScalarData(X);\nw = Nodes(Dual,(nx,ny));\nHmat = RegularizationMatrix(H,g,w);\nfill!(g,1.0);\nw .= Hmat*g;In general, the interpolation matrix is separately constructed, and the source and target are reversed:Emat = InterpolationMatrix(H,w,g);\ng .= Emat*w;Alternatively, if the regularization and interpolation are symmetric, then we can get them both when we call for the regularization matrix:H = Regularize(X,dx,issymmetric=true)\nHmat, Emat = RegularizationMatrix(H,g,w)It might seem a bit funny to store them separately if they are just transposes of each other, but it is essential for the method dispatch that they are given separate types."
},

{
    "location": "manual/fields/#ViscousFlow.Fields.CircularConvolution",
    "page": "Fields",
    "title": "ViscousFlow.Fields.CircularConvolution",
    "category": "type",
    "text": "CircularConvolution{M, N}\n\nA preplanned, circular convolution operator on an M × N matrix.\n\nFields\n\nĜ: DFT coefficients of the convolution kernel\nF: preplanned rFFT operator\nF⁻¹: preplanned irFFT operator\npaddedSpace: scratch space to zero-pad the input matrix\nÂ: scratch space to store the DFT coefficients of the zero-padded input matrix\n\nConstructors:\n\nCircularConvolution(G::Matrix{Float64})\n\nExample:\n\njulia> G = repeat(1.0:3,1,4)\n3×4 Array{Float64,2}:\n 1.0  1.0  1.0  1.0\n 2.0  2.0  2.0  2.0\n 3.0  3.0  3.0  3.0\n\njulia> C = CircularConvolution(G)\nCircular convolution on a 3 × 4 matrix\n\njulia> C*reshape(1:12, 3, 4)\n3×4 Array{Int64,2}:\n 164  164  164  164\n 130  130  130  130\n 148  148  148  148\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.DDF-Tuple{}",
    "page": "Fields",
    "title": "ViscousFlow.Fields.DDF",
    "category": "method",
    "text": "DDF([ddftype=Fields.Roma],[dx=1.0])\n\nConstruct a discrete delta function operator. This is generally only needed internally by the Regularize operator, so the user doesn\'t have much need for accessing this directly. The default DDF is the Roma function, which has a support of 3 grid cells. Other choices are the Goza operator, which is a truncated Gaussian with 28 cells support, and the Witchhat, which has 2 cells support. The resulting operator is evaluated with one, two or three coordinate arguments, producing, respectively, 1-d, 2-d, or 3-d smeared delta functions. It can also be called with the usual Julia vectorized dot notation with arrays of arguments. The optional cell spacing argument dx rescales the coordinates by this spacing, and the result is also rescaled by this spacing (raised to the number of dimensions). This spacing argument defaults to 1.0.\n\njulia> ddf = DDF(ddftype=Fields.Roma)\nDiscrete delta function operator of type ViscousFlow.Fields.Roma, with spacing 1.0\n\njulia> ddf(1)\n0.16666666666666666\n\njulia> ddf(-1)\n0.16666666666666666\n\njulia> ddf.([-1,0,1])\n3-element Array{Float64,1}:\n 0.16666666666666666\n 0.6666666666666666\n 0.16666666666666666\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.Edges",
    "page": "Fields",
    "title": "ViscousFlow.Fields.Edges",
    "category": "type",
    "text": "Edges{Dual/Primal}\n\nEdges is a wrapper for vector-valued data that lie at the faces of either dual cells or primary cells. Edges type data have fields u and v for the components of the vector field. These are the normal components of the vector field on the vertical and horizontal faces of the corresponding cell.\n\nConstructors\n\nEdges(C,dims) creates a vector field of zeros in cells of type C (where C is either Dual or Primal), on a grid of dimensions dims. Note that dims represent the number of dual cells on the grid.\nEdges(C,w) performs the same construction, but uses existing field data w of Nodes type to determine the size of the grid.\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.InterpolationMatrix",
    "page": "Fields",
    "title": "ViscousFlow.Fields.InterpolationMatrix",
    "category": "type",
    "text": "InterpolationMatrix(H::Regularize,u::Nodes/Edges,f::ScalarData/VectorData) -> Emat\n\nConstruct and store a matrix representation of interpolation associated with H for data of type u to data of type f. The resulting matrix Emat can then be used to apply on grid data of type u to interpolate it to point data of type f, using mul!(f,Emat,u). It can also be used as just Emat*u.\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.Nodes",
    "page": "Fields",
    "title": "ViscousFlow.Fields.Nodes",
    "category": "type",
    "text": "Nodes{Dual/Primal}\n\nNodes is a wrapper for scalar-valued data that lie at the centers of either dual cells or primary cells. A Nodes type can be accessed by indexing like any other array, and allows the use of [size].\n\nConstructors\n\nNodes(C,dims) creates a field of zeros in cells of type C (where C is either Dual or Primal), on a grid of dimensions dims. Note that dims represent the number of dual cells on the grid, even if C is Primal.\nNodes(C,w) performs the same construction, but uses existing field data w of Nodes type to determine the size of the grid.\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.RegularizationMatrix",
    "page": "Fields",
    "title": "ViscousFlow.Fields.RegularizationMatrix",
    "category": "type",
    "text": "RegularizationMatrix(H::Regularize,f::ScalarData/VectorData,u::Nodes/Edges) -> Hmat\n\nConstruct and store a matrix representation of regularization associated with H for data of type f to data of type u. The resulting matrix Hmat can then be used to apply on point data of type f to regularize it to grid data of type u, using mul!(u,Hmat,f). It can also be used as just Hmat*f.\n\nIf H is a symmetric regularization and interpolation operator, then this actually returns a tuple Hmat, Emat, where Emat is the interpolation matrix.\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.Regularize-Union{Tuple{T}, Tuple{Array{T,1},Array{T,1},T}} where T<:Real",
    "page": "Fields",
    "title": "ViscousFlow.Fields.Regularize",
    "category": "method",
    "text": "Regularize(x,y,dx,[ddftype=Roma],[I0=(1,1)], [weights=1.0], [filter=false],\n                   [issymmetric=false])\n\nConstructor to set up an operator for regularizing and interpolating data from/to points immersed in the grid to/from fields on the grid itself. The supplied x and y represent physical coordinates of the immersed points, and dx denotes a uniform physical cell size of the grid. The separate arguments x and y can be replaced by a single argument X of type VectorData holding the coordinates.\n\nThe operations of regularization and interpolation are carried out with a discrete delta function (ddf), which defaults to the type Roma. Others are also possible, such as Goza. The optional tuple I0 represents the indices of the primary node that coincides with (x,y) = (0,0). This defaults to (1,1), which leaves one layer of ghost (dual) cells and sets the physical origin in the lower left corner of the grid of interior dual cells.\n\nAnother optional parameter, weights, sets the weight of each point in the regularization. This would generally be set with, say, the differential arc length for regularization of data on a curve. It can be a vector (of the same length as x and y) or a scalar if uniform. It defaults to 1.0.\n\nThe optional Boolean parameter filter can be set to true if it is desired to apply filtering (see Goza et al, J Comput Phys 2016) to the grid data before interpolating. This is generally only used in the context of preconditioning the solution for forces on the immersed points.\n\nIf the optional Boolean parameter issymmetric is set to true, then the regularization and interpolation are constructed to be transposes of each other. Note that this option overrides any supplied weights. The default of this parameter is false.\n\nThe resulting operator can be used in either direction, regularization and interpolation, with the first argument representing the target (the entity to regularize/interpolate to), and the second argument the source (the entity to regularize/interpolate from). The regularization does not use the filtering option.\n\nExample\n\nIn the example below, we set up a 12 x 12 grid. Using the default value for I0 and setting dx = 0.1, the physical dimensions of the non-ghost part of the grid are 1.0 x 1.0. Three points are set up in the interior, and a vector field is assigned to them, with the x component of each of them set to 1.0. These data are regularized to a field of primal edges on the grid.\n\njulia> x = [0.25,0.75,0.25]; y = [0.75,0.25,0.25];\n\njulia> X = VectorData(x,y);\n\njulia> q = Edges(Primal,(12,12));\n\njulia> dx = 0.1;\n\njulia> H = Regularize(x,y,dx)\nRegularization/interpolation operator with non-filtered interpolation\n  3 points in grid with cell area 0.01\n\njulia> f = VectorData(X);\n\njulia> fill!(f.u,1.0);\n\njulia> H(q,f)\nEdges{Primal,12,12} data\nu (in grid orientation)\n11×12 Array{Float64,2}:\n 0.0  0.0  0.0       0.0     0.0      …  0.0       0.0     0.0      0.0  0.0\n 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0\n 0.0  0.0  8.33333  33.3333  8.33333     0.0       0.0     0.0      0.0  0.0\n 0.0  0.0  8.33333  33.3333  8.33333     0.0       0.0     0.0      0.0  0.0\n 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0\n 0.0  0.0  0.0       0.0     0.0      …  0.0       0.0     0.0      0.0  0.0\n 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0\n 0.0  0.0  8.33333  33.3333  8.33333     8.33333  33.3333  8.33333  0.0  0.0\n 0.0  0.0  8.33333  33.3333  8.33333     8.33333  33.3333  8.33333  0.0  0.0\n 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0\n 0.0  0.0  0.0       0.0     0.0      …  0.0       0.0     0.0      0.0  0.0\nv (in grid orientation)\n12×11 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.ScalarData",
    "page": "Fields",
    "title": "ViscousFlow.Fields.ScalarData",
    "category": "type",
    "text": "ScalarData\n\nA wrapper for a one-dimensional array of scalar-valued data. The resulting wrapper can be indexed in the same way as the array itself.\n\nConstructors\n\nScalarData(d) constructs a wrapper for the one-dimensional array of data d\nScalarData(n::Int) constructs a wrapper for an array of zeros of length n.\nScalarData(x::ScalarData) constructs a wrapper for an array of zeros of the  same length as that wrapped by x.\nScalarData(x::VectorData) constructs a wrapper for an array of zeros of the   same length as that wrapped by x.\n\nExample\n\njulia> f = ScalarData(10);\n\njulia> f[5] = 1.0;\n\njulia> f\n10 points of scalar-valued data\n10-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n 1.0\n 0.0\n 0.0\n 0.0\n 0.0\n 0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.VectorData",
    "page": "Fields",
    "title": "ViscousFlow.Fields.VectorData",
    "category": "type",
    "text": "VectorData\n\nA wrapper for a one-dimensional array of two-component vector-valued data. The resulting wrapper can be indexed as though the first component and second component are stacked on top of each other.\n\nConstructors\n\nVectorData(u,v) constructs a wrapper for the vector components data u and v.\nVectorData(n::Int) constructs a wrapper with zeros of length n for both components.\nVectorData(x::ScalarData) constructs a wrapper for zero components of the  same length as that wrapped by x.\nVectorData(x::VectorData) constructs a wrapper for zero components of the   same length as that wrapped by x.\n\nExample\n\njulia> f = VectorData(10);\n\njulia> f.v[1:5] = 1:5;\n\njulia> f\n10 points of vector-valued data\n10×2 Array{Float64,2}:\n 0.0  1.0\n 0.0  2.0\n 0.0  3.0\n 0.0  4.0\n 0.0  5.0\n 0.0  0.0\n 0.0  0.0\n 0.0  0.0\n 0.0  0.0\n 0.0  0.0\n\njulia> f[7] = 1; f[18] = 0.2;\n\njulia> f\n10 points of vector-valued data\n10×2 Array{Float64,2}:\n 0.0  1.0\n 0.0  2.0\n 0.0  3.0\n 0.0  4.0\n 0.0  5.0\n 0.0  0.0\n 1.0  0.0\n 0.0  0.2\n 0.0  0.0\n 0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.cellshift!-Union{Tuple{NY}, Tuple{NX}, Tuple{Edges{Dual,NX,NY},Edges{Primal,NX,NY}}} where NY where NX",
    "page": "Fields",
    "title": "ViscousFlow.Fields.cellshift!",
    "category": "method",
    "text": "cellshift!(v::Edges{Dual/Primal},q::Edges{Primal/Dual})\n\nShift (by linear interpolation) the primal (resp. dual) edge data q to the edges of the dual (resp. primal) cells, and return the result in v.\n\nExample\n\njulia> q = Edges(Primal,(8,6));\n\njulia> q.u[3,2] = 1.0;\n\njulia> v = Edges(Dual,(8,6));\n\njulia> Fields.cellshift!(v,q)\nEdges{Dual,8,6} data\nu (in grid orientation)\n6×7 Array{Float64,2}:\n 0.0  0.0   0.0   0.0  0.0  0.0  0.0\n 0.0  0.0   0.0   0.0  0.0  0.0  0.0\n 0.0  0.0   0.0   0.0  0.0  0.0  0.0\n 0.0  0.25  0.25  0.0  0.0  0.0  0.0\n 0.0  0.25  0.25  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0   0.0  0.0  0.0  0.0\nv (in grid orientation)\n5×8 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.cellshift!-Union{Tuple{NY}, Tuple{NX}, Tuple{Edges{Dual,NX,NY},Nodes{Dual,NX,NY}}} where NY where NX",
    "page": "Fields",
    "title": "ViscousFlow.Fields.cellshift!",
    "category": "method",
    "text": "cellshift!(q::Edges{Dual},w::Nodes{Dual})\n\nShift (by linear interpolation) the dual nodal data w to the edges of the dual cells, and return the result in q.\n\nExample\n\njulia> w = Nodes(Dual,(8,6));\n\njulia> w[3,4] = 1.0;\n\njulia> q = Edges(Dual,w);\n\njulia> cellshift!(q,w)\nEdges{Dual,8,6} data\nu (in grid orientation)\n6×7 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.5  0.5  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\nv (in grid orientation)\n5×8 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.cellshift!-Union{Tuple{NY}, Tuple{NX}, Tuple{Tuple{Nodes{Dual,NX,NY},Nodes{Dual,NX,NY}},Edges{Primal,NX,NY}}} where NY where NX",
    "page": "Fields",
    "title": "ViscousFlow.Fields.cellshift!",
    "category": "method",
    "text": "cellshift!((wx::Nodes,wy::Nodes),q::Edges)\n\nShift (by linear interpolation) the edge data q (of either dual or primal type) to the dual or primal nodes, and return the result in wx and wy. wx holds the shifted q.u data and wy the shifted q.v data.\n\nExample\n\njulia> q = Edges(Primal,(8,6));\n\njulia> q.u[3,2] = 1.0;\n\njulia> wx = Nodes(Dual,(8,6)); wy = deepcopy(wx);\n\njulia> Fields.cellshift!((wx,wy),q);\n\njulia> wx\nNodes{Dual,8,6} data\nPrinting in grid orientation (lower left is (1,1))\n6×8 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\njulia> wy\nNodes{Dual,8,6} data\nPrinting in grid orientation (lower left is (1,1))\n6×8 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.coordinates",
    "page": "Fields",
    "title": "ViscousFlow.Fields.coordinates",
    "category": "function",
    "text": "cooordinates(w::Nodes/Edges;[dx=1.0],[I0=(1,1)])\n\nReturn a tuple of the ranges of the physical coordinates in each direction for grid data w. If w is of Nodes type, then it returns a tuple of the form xg,yg. If w is of Edges or NodePair type, then it returns a tuple of the form xgu,ygu,xgv,ygv.\n\nThe optional keyword argument dx sets the grid spacing; its default is 1.0. The optional keyword I0 accepts a tuple of integers to set the index pair of the primal nodes that coincide with the origin. The default is (1,1).\n\nExample\n\njulia> w = Nodes(Dual,(12,22));\n\njulia> xg, yg = coordinates(w,dx=0.1)\n(-0.05:0.1:1.05, -0.05:0.1:2.0500000000000003)\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.curl!-Union{Tuple{NY}, Tuple{NX}, Tuple{Edges{Primal,NX,NY},Nodes{Dual,NX,NY}}} where NY where NX",
    "page": "Fields",
    "title": "ViscousFlow.Fields.curl!",
    "category": "method",
    "text": "curl!(q::Edges{Primal},w::Nodes{Dual})\n\nEvaluate the discrete curl of w and return it as q.\n\nExample\n\njulia> w = Nodes(Dual,(8,6));\n\njulia> w[3,4] = 1.0;\n\njulia> q = Edges(Primal,w);\n\njulia> curl!(q,w)\nEdges{Primal,8,6} data\nu (in grid orientation)\n5×8 Array{Float64,2}:\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  -1.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   1.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0\nv (in grid orientation)\n6×7 Array{Float64,2}:\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  -1.0  1.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.curl-Union{Tuple{Nodes{Dual,NX,NY}}, Tuple{NY}, Tuple{NX}} where NY where NX",
    "page": "Fields",
    "title": "ViscousFlow.Fields.curl",
    "category": "method",
    "text": "curl(w::Nodes{Dual}) --> Edges{Primal}\n\nEvaluate the discrete curl of w. Another way to perform this operation is to construct a Curl object and apply it with *.\n\nExample\n\njulia> C = Curl();\n\njulia> w = Nodes(Dual,(8,6));\n\njulia> w[3,4] = 1.0;\n\njulia> C*w\nEdges{Primal,8,6} data\nu (in grid orientation)\n5×8 Array{Float64,2}:\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  -1.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   1.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0\nv (in grid orientation)\n6×7 Array{Float64,2}:\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  -1.0  1.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.divergence!-Union{Tuple{NY}, Tuple{NX}, Tuple{Nodes{Primal,NX,NY},Edges{Primal,NX,NY}}} where NY where NX",
    "page": "Fields",
    "title": "ViscousFlow.Fields.divergence!",
    "category": "method",
    "text": "divergence!(w::Nodes,q::Edges)\n\nEvaluate the discrete divergence of edge data q and return it as nodal data w. Note that q can be either primal or dual edge data, but w must be of the same cell type.\n\nExample\n\njulia> q = Edges(Primal,(8,6));\n\njulia> q.u[3,2] = 1.0;\n\njulia> w = Nodes(Primal,(8,6));\n\njulia> divergence!(w,q)\nNodes{Primal,8,6} data\nPrinting in grid orientation (lower left is (1,1))\n5×7 Array{Float64,2}:\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  1.0  -1.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.divergence-Union{Tuple{Edges{T,NX,NY}}, Tuple{NY}, Tuple{NX}, Tuple{T}} where NY where NX where T<:ViscousFlow.Fields.CellType",
    "page": "Fields",
    "title": "ViscousFlow.Fields.divergence",
    "category": "method",
    "text": "divergence(q::Edges) --> Nodes\n\nEvaluate the discrete divergence of edge data q. Can also perform this operation by creating an object of Divergence type and applying it with *.\n\nExample\n\njulia> D = Divergence();\n\njulia> q = Edges(Primal,(8,6));\n\njulia> q.u[3,2] = 1.0;\n\njulia> D*q\nNodes{Primal,8,6} data\nPrinting in grid orientation (lower left is (1,1))\n5×7 Array{Float64,2}:\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  1.0  -1.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.grad!-Union{Tuple{NY}, Tuple{NX}, Tuple{Edges{Primal,NX,NY},Nodes{Primal,NX,NY}}} where NY where NX",
    "page": "Fields",
    "title": "ViscousFlow.Fields.grad!",
    "category": "method",
    "text": "grad!(q::Edges{Primal},w::Nodes{Primal})\n\nEvaluate the discrete gradient of primal nodal data w and return it as primal edge data q.\n\nExample\n\njulia> w = Nodes(Primal,(8,6));\n\njulia> w[3,4] = 1.0;\n\njulia> q = Edges(Primal,(8,6));\n\njulia> grad!(q,w)\nEdges{Primal,8,6} data\nu (in grid orientation)\n5×8 Array{Float64,2}:\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  1.0  -1.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\nv (in grid orientation)\n6×7 Array{Float64,2}:\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  -1.0  0.0  0.0  0.0  0.0\n 0.0  0.0   1.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.grad-Union{Tuple{Nodes{Primal,NX,NY}}, Tuple{NY}, Tuple{NX}} where NY where NX",
    "page": "Fields",
    "title": "ViscousFlow.Fields.grad",
    "category": "method",
    "text": "grad(w::Nodes{Primal}) --> Edges{Primal}\n\nEvaluate the discrete gradient of primal nodal data w. Can also perform this operation by creating an object of Grad type and applying it with *.\n\nExample\n\njulia> w = Nodes(Primal,(8,6));\n\njulia> w[3,4] = 1.0;\n\njulia> G = Grad();\n\njulia> G*w\nEdges{Primal,8,6} data\nu (in grid orientation)\n5×8 Array{Float64,2}:\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  1.0  -1.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\nv (in grid orientation)\n6×7 Array{Float64,2}:\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  -1.0  0.0  0.0  0.0  0.0\n 0.0  0.0   1.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0   0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.laplacian!-Union{Tuple{NY}, Tuple{NX}, Tuple{Nodes{Dual,NX,NY},Nodes{Dual,NX,NY}}} where NY where NX",
    "page": "Fields",
    "title": "ViscousFlow.Fields.laplacian!",
    "category": "method",
    "text": "laplacian!(v,w)\n\nEvaluate the discrete Laplacian of w and return it as v. The data w can be of type dual/primary nodes or edges; v must be of the same type.\n\nExample\n\njulia> w = Nodes(Dual,(8,6));\n\njulia> v = deepcopy(w);\n\njulia> w[4,3] = 1.0;\n\njulia> laplacian!(v,w)\nNodes{Dual,8,6} data\nPrinting in grid orientation (lower left is (1,1))\n6×8 Array{Float64,2}:\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   1.0  0.0  0.0  0.0  0.0\n 0.0  0.0  1.0  -4.0  1.0  0.0  0.0  0.0\n 0.0  0.0  0.0   1.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.laplacian-Union{Tuple{Nodes{T,NX,NY}}, Tuple{NY}, Tuple{NX}, Tuple{T}} where NY where NX where T<:ViscousFlow.Fields.CellType",
    "page": "Fields",
    "title": "ViscousFlow.Fields.laplacian",
    "category": "method",
    "text": "laplacian(w)\n\nEvaluate the discrete Laplacian of w. The data w can be of type dual/primary nodes or edges. The returned result is of the same type as w.\n\nExample\n\njulia> q = Edges(Primal,(8,6));\n\njulia> q.u[2,2] = 1.0;\n\njulia> laplacian(q)\nEdges{Primal,8,6} data\nu (in grid orientation)\n5×8 Array{Float64,2}:\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0   1.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  -4.0  1.0  0.0  0.0  0.0  0.0  0.0\n 0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0\nv (in grid orientation)\n6×7 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.plan_intfact",
    "page": "Fields",
    "title": "ViscousFlow.Fields.plan_intfact",
    "category": "function",
    "text": "plan_intfact(a::Real,dims::Tuple,[fftw_flags=FFTW.ESTIMATE])\n\nConstructor to set up an operator for evaluating the integrating factor with real-valued parameter a. This can then be applied with the * operation on data of the appropriate size.\n\nThe dims argument can be replaced with w::Nodes to specify the size of the domain.\n\nExample\n\njulia> w = Nodes(Dual,(6,6));\n\njulia> w[4,4] = 1.0;\n\njulia> E = plan_intfact(1.0,(6,6))\nIntegrating factor with parameter 1.0 on a (nx = 6, ny = 6) grid\n\njulia> E*w\nNodes{Dual,6,6} data\nPrinting in grid orientation (lower left is (1,1))\n6×6 Array{Float64,2}:\n 0.00268447   0.00869352  0.0200715   0.028765    0.0200715   0.00869352\n 0.00619787   0.0200715   0.0463409   0.0664124   0.0463409   0.0200715\n 0.00888233   0.028765    0.0664124   0.0951774   0.0664124   0.028765\n 0.00619787   0.0200715   0.0463409   0.0664124   0.0463409   0.0200715\n 0.00268447   0.00869352  0.0200715   0.028765    0.0200715   0.00869352\n 0.000828935  0.00268447  0.00619787  0.00888233  0.00619787  0.00268447\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.plan_intfact!",
    "page": "Fields",
    "title": "ViscousFlow.Fields.plan_intfact!",
    "category": "function",
    "text": "plan_intfact!(a::Real,dims::Tuple,[fftw_flags=FFTW.ESTIMATE])\n\nSame as plan_intfact, but operates in-place on data.\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.plan_laplacian",
    "page": "Fields",
    "title": "ViscousFlow.Fields.plan_laplacian",
    "category": "function",
    "text": "plan_laplacian(dims::Tuple,[with_inverse=false],[fftw_flags=FFTW.ESTIMATE],\n                      [dx=1.0])\n\nConstructor to set up an operator for evaluating the discrete Laplacian on dual or primal nodal data of dimension dims. If the optional keyword with_inverse is set to true, then it also sets up the inverse Laplacian (the lattice Green\'s function, LGF). These can then be applied, respectively, with * and \\ operations on data of the appropriate size. The optional parameter dx is used in adjusting the uniform value of the LGF to match the behavior of the continuous analog at large distances; this is set to 1.0 by default.\n\nInstead of the first argument, one can also supply w::Nodes to specify the size of the domain.\n\nExample\n\njulia> w = Nodes(Dual,(5,5));\n\njulia> w[3,3] = 1.0;\n\njulia> L = plan_laplacian(5,5;with_inverse=true)\nDiscrete Laplacian (and inverse) on a (nx = 5, ny = 5) grid with spacing 1.0\n\njulia> s = L\\w\nNodes{Dual,5,5} data\nPrinting in grid orientation (lower left is (1,1))\n5×5 Array{Float64,2}:\n 0.16707    0.129276     0.106037     0.129276    0.16707\n 0.129276   0.0609665   -0.00734343   0.0609665   0.129276\n 0.106037  -0.00734343  -0.257343    -0.00734343  0.106037\n 0.129276   0.0609665   -0.00734343   0.0609665   0.129276\n 0.16707    0.129276     0.106037     0.129276    0.16707\n\njulia> L*s ≈ w\ntrue\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.plan_laplacian!",
    "page": "Fields",
    "title": "ViscousFlow.Fields.plan_laplacian!",
    "category": "function",
    "text": "plan_laplacian!(dims::Tuple,[with_inverse=false],[fftw_flags=FFTW.ESTIMATE],\n                      [dx=1.0])\n\nSame as plan_laplacian, but operates in-place on data.\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.product!-Union{Tuple{NY}, Tuple{NX}, Tuple{T}, Tuple{Edges{T,NX,NY},Edges{T,NX,NY},Edges{T,NX,NY}}} where NY where NX where T",
    "page": "Fields",
    "title": "ViscousFlow.Fields.product!",
    "category": "method",
    "text": "product!(out::Edges/Nodes,p::Edges/Nodes,q::Edges/Nodes)\n\nCompute the Hadamard (i.e. element by element) product of edge or nodal (primal or dual) data p and q and return the result in out.\n\nExample\n\njulia> q = Edges(Dual,(8,6));\n\njulia> out = p = deepcopy(q);\n\njulia> q.u[3,2] = 0.3;\n\njulia> p.u[3,2] = 0.2;\n\njulia> product!(out,p,q)\nEdges{Dual,8,6} data\nu (in grid orientation)\n6×7 Array{Float64,2}:\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\n 0.0  0.0  0.06  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\nv (in grid orientation)\n5×8 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#ViscousFlow.Fields.product-Union{Tuple{NY}, Tuple{NX}, Tuple{T}, Tuple{Edges{T,NX,NY},Edges{T,NX,NY}}} where NY where NX where T",
    "page": "Fields",
    "title": "ViscousFlow.Fields.product",
    "category": "method",
    "text": "product(p::Edges/Nodes,q::Edges/Nodes) --> Edges/Nodes\n\nCompute the Hadamard product of edge or nodal (primal or dual) data p and q and return the result. This operation can also be carried out with the ∘ operator:\n\nExample\n\njulia> q = Edges(Dual,(8,6));\n\njulia> p = deepcopy(q);\n\njulia> q.u[3,2] = 0.3;\n\njulia> p.u[3,2] = 0.2;\n\njulia> p∘q\nEdges{Dual,8,6} data\nu (in grid orientation)\n6×7 Array{Float64,2}:\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\n 0.0  0.0  0.06  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0   0.0  0.0  0.0  0.0\nv (in grid orientation)\n5×8 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "manual/fields/#Methods-1",
    "page": "Fields",
    "title": "Methods",
    "category": "section",
    "text": "Modules = [Fields]\nOrder   = [:type, :function]"
},

{
    "location": "manual/fields/#Index-1",
    "page": "Fields",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"fields.md\"]"
},

{
    "location": "manual/bodies/#",
    "page": "Bodies",
    "title": "Bodies",
    "category": "page",
    "text": ""
},

{
    "location": "manual/bodies/#Bodies-1",
    "page": "Bodies",
    "title": "Bodies",
    "category": "section",
    "text": "DocTestSetup = quote\nusing ViscousFlow\nusing Random\nRandom.seed!(1)\nendusing ViscousFlow\nusing Plots"
},

{
    "location": "manual/bodies/#ViscousFlow.Bodies.Ellipse",
    "page": "Bodies",
    "title": "ViscousFlow.Bodies.Ellipse",
    "category": "type",
    "text": "Ellipse(a,b,n)\n\nConstruct an elliptical body with semi-major axis a and semi-minor axis b, with n points distributed on the body perimeter.\n\nThe constructor Ellipse(a,n) creates a circle of radius a.\n\n\n\n\n\n"
},

{
    "location": "manual/bodies/#ViscousFlow.Bodies.Plate",
    "page": "Bodies",
    "title": "ViscousFlow.Bodies.Plate",
    "category": "type",
    "text": "Plate(length,thick,n,[λ=1.0])\n\nConstruct a flat plate with length length and thickness thick, with n points distributed on the body perimeter.\n\nThe optional parameter λ distributes the points differently. Values between 0.0 and 1.0 are accepted.\n\nThe constructor Plate(length,n,[λ=1.0]) creates a plate of zero thickness.\n\n\n\n\n\n"
},

{
    "location": "manual/bodies/#ViscousFlow.Bodies.RigidTransform",
    "page": "Bodies",
    "title": "ViscousFlow.Bodies.RigidTransform",
    "category": "type",
    "text": "RigidTransform(x::Tuple{Float64,Float64},α::Float64)\n\nConstruct a rigid-body transform operator, with rotation by angle α and translation specified by x. The translation coordinates are specified in the target coordinate system.\n\nThe resulting transform can be used as an operator on pairs of coordinate vectors, x and y, or on bodies. For transformation of bodies, it only overwrites the x and y fields of the body, but leaves the x̃ and ỹ (body coordinates) intact.\n\nThe translation can be provided as either a tuple (x,y) or as a complex number.\n\nExample\n\njulia> body = Bodies.Ellipse(0.5,0.1,100)\nElliptical body with 100 points and semi-axes (0.5,0.1)\n   Current position: (0.0,0.0)\n   Current angle (rad): 0.0\n\njulia> T = RigidTransform((1.0,1.0),π/4)\nRigid-body transform\n  Translation: (1.0,1.0)\n  Rotation angle (rad): 0.7853981633974483\n\njulia> T(body)\nElliptical body with 100 points and semi-axes (0.5,0.1)\n   Current position: (1.0,1.0)\n   Current angle (rad): 0.7853981633974483\n\n\n\n\n\n"
},

{
    "location": "manual/bodies/#ViscousFlow.Bodies.NACA4",
    "page": "Bodies",
    "title": "ViscousFlow.Bodies.NACA4",
    "category": "type",
    "text": "NACA4(cam,pos,thick[;np=20][,Zc=0.0+0.0im][,len=1.0]) -> Vector{ComplexF64}\n\nGenerates the vertices of a NACA 4-digit airfoil of chord length 1. The relative camber is specified by cam, the position of maximum camber (as fraction of chord) by pos, and the relative thickness by thick.\n\nThe optional parameter np specifies the number of points on the upper or lower surface. The optional parameter Zc specifies the mean position of the vertices (which is set to the origin by default). The optional parameter len specifies the chord length.\n\nExample\n\njulia> w = Bodies.NACA4(0.0,0.0,0.12);\n\n\n\n\n\n"
},

{
    "location": "manual/bodies/#Base.diff-Union{Tuple{Body{N}}, Tuple{N}} where N",
    "page": "Bodies",
    "title": "Base.diff",
    "category": "method",
    "text": "diff(body::Body) -> Tuple{Vector{Float64},Vector{Float64}}\n\nCompute the x and y differences of the faces on the perimeter of body body, whose centers are at the current x and y coordinates (in inertial space) of the body.\n\n\n\n\n\n"
},

{
    "location": "manual/bodies/#Base.length-Union{Tuple{Body{N}}, Tuple{N}} where N",
    "page": "Bodies",
    "title": "Base.length",
    "category": "method",
    "text": "length(body::Body)\n\nReturn the number of points on the body perimeter\n\n\n\n\n\n"
},

{
    "location": "manual/bodies/#ViscousFlow.Bodies.dlength-Union{Tuple{Body{N}}, Tuple{N}} where N",
    "page": "Bodies",
    "title": "ViscousFlow.Bodies.dlength",
    "category": "method",
    "text": "dlength(body::Body) -> Vector{Float64}\n\nCompute the lengths of the faces on the perimeter of body body, whose centers are at the current x and y coordinates (in inertial space) of the body.\n\n\n\n\n\n"
},

{
    "location": "manual/bodies/#ViscousFlow.Bodies.normal-Union{Tuple{Body{N}}, Tuple{N}} where N",
    "page": "Bodies",
    "title": "ViscousFlow.Bodies.normal",
    "category": "method",
    "text": "normal(body::Body) -> Tuple{Vector{Float64},Vector{Float64}}\n\nCompute the current normals (in inertial components) of the faces on the perimeter of body body, whose centers are at the current x and y coordinates (in inertial space) of the body.\n\n\n\n\n\n"
},

{
    "location": "manual/bodies/#Methods-1",
    "page": "Bodies",
    "title": "Methods",
    "category": "section",
    "text": "Modules = [Bodies]\nOrder   = [:type, :function]"
},

{
    "location": "manual/bodies/#Index-1",
    "page": "Bodies",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"bodies.md\"]"
},

{
    "location": "manual/saddlesystems/#",
    "page": "Saddle point systems",
    "title": "Saddle point systems",
    "category": "page",
    "text": ""
},

{
    "location": "manual/saddlesystems/#Saddle-point-systems-1",
    "page": "Saddle point systems",
    "title": "Saddle point systems",
    "category": "section",
    "text": "DocTestSetup = quote\nusing ViscousFlow\nenddefddt1fracmathrmd1mathrmdt\n\nrenewcommandvecboldsymbol\nnewcommanduvec1vechat1\nnewcommandutangentuvectau\nnewcommandunormaluvecn\n\nrenewcommanddmathrmdusing ViscousFlow\nusing PlotsSaddle systems comprise an important part of solving mechanics problems with constraints. In such problems, there is an underlying system to solve, and the addition of constraints requires that the system is subjected to additional forces (constraint forces, or Lagrange multipliers) that enforce these constraints in the system. Examples of such constrained systems are the divergence-free velocity constraint in incompressible flow (for which pressure is the associated Lagrange multiplier field), the no-slip and/or no-flow-through condition in general fluid systems adjacent to impenetrable bodies, and joint constraints in rigid-body mechanics.A general saddle-point system has the formleft beginarraycc A  B_1^T  B_2  0endarrayright left(beginarraycuf endarrayright) = left(beginarraycr_1r_2 endarrayright)We are primarily interested in cases when the operator A is symmetric and positive definite, which is fairly typical. It is also fairly common for B_1 = B_2, so that the whole system is symmetric.whirl allows us to solve such systems for u and f in a fairly easy way. We need only to provide rules for how to evaluate the actions of the various operators in the system. Let us use an example to show how this can be done."
},

{
    "location": "manual/saddlesystems/#Translating-cylinder-in-potential-flow-1",
    "page": "Saddle point systems",
    "title": "Translating cylinder in potential flow",
    "category": "section",
    "text": "In irrotational, incompressible flow, the streamfunction psi satisfies Laplace\'s equation,nabla^2 psi = 0On the surface of an impenetrable body, the streamfunction must obey the constraintpsi = psi_bwhere psi_b is the streamfunction associated with the body\'s motion. Let us suppose the body is moving vertically with velocity 1. Then psi_b = -x for all points inside or on the surface of the body. Thus, the streamfunction field outside this body is governed by Laplace\'s equation subject to the constraint.Let us solve this problem on a staggered grid, using the tools discussed in the Fields section, including the regularization and interpolation methods to immerse the body shape on the grid. Then our saddle-point system has the formleft beginarraycc L  H  E  0endarrayright left(beginarraycpsif endarrayright) = left(beginarrayc0psi_b endarrayright)where L is the discrete Laplacian, H is the regularization operator, and E is the interpolation operator.Physically, f isn\'t really a force here, but rather, represents the strengths of distributed singularities on the surface. In fact, this strength represents the jump in normal derivative of psi across the surface. Since this normal derivative is equivalent to the tangential velocity, f is the strength of the bound vortex sheet on the surface. This will be useful to know when we check the value of f obtained in our solution.First, let us set up the body, centered at (11) and of radius 12. We will also initialize a data structure for the force:using ViscousFlow\nusing Plots\npyplot()n = 128; θ = range(0,stop=2π,length=n+1);\nxb = 1.0 .+ 0.5*cos.(θ[1:n]); yb = 1.0 .+ 0.5*sin.(θ[1:n]);\nX = VectorData(xb,yb);\nf = ScalarData(X);Now let\'s set up a grid of size 102times 102 (including the usual layer of ghost cells) and physical dimensions 2times 2.nx = 102; ny = 102; Lx = 2.0; dx = Lx/(nx-2);\nw = Nodes(Dual,(nx,ny));We need to set up the operators now. First, the Laplacian:L = plan_laplacian(size(w),with_inverse=true)\nL⁻¹(w::T) where {T} = L\\wThe last line just defines another operator for computing the inverse of L. We have called it L⁻¹ for useful shorthand. This operator acts upon dual nodal data and returns data of the same type, e.g. ψ = L⁻¹(w). The saddle point system structure requires operators that have this sort of form.Now we need to set up the regularization H and interpolation E operators.regop = Regularize(X,dx;issymmetric=true)\nHmat, Emat = RegularizationMatrix(regop,f,w);Now we are ready to set up the system.S = SaddleSystem((w,f),(L⁻¹,Hmat,Emat),issymmetric=true,isposdef=true)Note that we have provided a tuple of the types of data, w and f, that we want the solver to work with, along with a tuple of the definitions of the three operators. The operators can be in the form of a function acting on its data (as for L⁻¹) or in the form of a matrix (or matrix-like) operator (as for Hmat and Emat); the constructor sorts it out. However, the order is important: we must supply A^-1, B_1^T, and B_2, in that order.We have also set two optional flags, to specify that the system is symmetric and positive definite. This instructs on which solver to use. (This is actually not quite true for the Laplacian: it is only positive semi-definite, since this operator has a null space. It is   adequate criteria for using the conjugate gradient method, but we will have to   be careful of some aspects of the solution, we we will see below.)Let\'s solve the system. We need to supply the right-hand side.w = Nodes(Dual,(nx,ny));\nψb = ScalarData(X);\nψb .= -(xb.-1);The right-hand side of the Laplace equation is zero. The right-hand side of the constraint is the specified streamfunction on the body. Note that we have subtracted the circle center from the x positions on the body. The reason for this will be discussed in a moment.We solve the system with the convenient shorthand of the backslash:@time ψ,f = S\\(w,ψb)Just to point out how fast it can be, we have also timed it. It\'s pretty fast.Now, let\'s plot the solution in physical space. We\'ll plot the body shape for reference, also.xg, yg = coordinates(ψ,dx=dx)\nplot(xg,yg,ψ)\nplot!(xb,yb,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black)\nsavefig(\"sfunc.svg\"); nothing # hide(Image: )The solution shows the streamlines for a circle in vertical motion, as expected. All of the streamlines inside the circle are vertical."
},

{
    "location": "manual/saddlesystems/#Core.Type-Union{Tuple{FP}, Tuple{FB2}, Tuple{FB1}, Tuple{FA}, Tuple{TF}, Tuple{TU}, Tuple{Tuple{TU,TF},Tuple{FA,FB1,FB2}}} where FP where FB2 where FB1 where FA where TF where TU",
    "page": "Saddle point systems",
    "title": "Core.Type",
    "category": "method",
    "text": "SaddleSystem((u,f),(A⁻¹,B₁ᵀ,B₂);[tol=1e-3],[issymmetric=false],[isposdef=false],[conditioner=Identity],[store=false])\n\nConstruct the computational operators for a saddle-point system of the form A B₁ᵀ B₂ 0uf. Note that the constituent operators are passed in as a tuple in the order seen here. Each of these operators could act on its corresponding data type in a function-like way, e.g. A⁻¹(u), or in a matrix-like way, e.g., A⁻¹*u.\n\nThe optional argument tol sets the tolerance for iterative solution (if   applicable). Its default is 1e-3.\n\nThe optional argument conditioner can be used to supply a function that acts upon the result f to \'condition\' it (e.g. filter it). It is, by default, set to the identity.\n\nThe optional Boolean argument store will compute and store the Schur complement matrix\'s factorization. This makes the inversion faster, though it comes at the expense of memory and overhead time for pre-computing it. The resulting solution is somewhat noiser, too.\n\nArguments\n\nu : example of state vector data.\nf : example of constraint force vector data. This data must be of       AbstractVector supertype.\nA⁻¹ : operator evaluating the inverse of A on data of type u, return type u\nB₁ᵀ : operator evaluating the influence of constraint force,           acting on f and returning type u\nB₂ : operator evaluating the influence of state vector on constraints,           acting on u and returning type f\n\n\n\n\n\n"
},

{
    "location": "manual/saddlesystems/#LinearAlgebra.ldiv!-Union{Tuple{N}, Tuple{FP}, Tuple{FBA}, Tuple{FAB}, Tuple{FA}, Tuple{TF}, Tuple{TU}, Tuple{Tuple{TU,TF},SaddleSystem{FA,FAB,FBA,FP,TU,TF,N,false},Tuple{TU,TF}}} where N where FP where FBA where FAB where FA where TF where TU",
    "page": "Saddle point systems",
    "title": "LinearAlgebra.ldiv!",
    "category": "method",
    "text": "ldiv!(state,sys::SaddleSystem,rhs)\n\nSolve a saddle-point system. rhs is a tuple of the right-hand side (ru,rf). Output state, a tuple (u,f), is updated. Note that sys is also mutated: its scratch space sys.B₂A⁻¹r₁ and sys.A⁻¹B₁ᵀf hold the intermediate results of the solution.\n\nA shorthand can be used for this operation: state = sys hs\n\n\n\n\n\n"
},

{
    "location": "manual/saddlesystems/#Methods-1",
    "page": "Saddle point systems",
    "title": "Methods",
    "category": "section",
    "text": "Modules = [SaddlePointSystems]\nOrder   = [:type, :function]"
},

{
    "location": "manual/saddlesystems/#Index-1",
    "page": "Saddle point systems",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"saddlesystems.md\"]"
},

{
    "location": "manual/timemarching/#",
    "page": "Time marching",
    "title": "Time marching",
    "category": "page",
    "text": ""
},

{
    "location": "manual/timemarching/#Time-marching-1",
    "page": "Time marching",
    "title": "Time marching",
    "category": "section",
    "text": "DocTestSetup = quote\nusing ViscousFlow\nenddefddt1fracmathrmd1mathrmdt\n\nrenewcommandvecboldsymbol\nnewcommanduvec1vechat1\nnewcommandutangentuvectau\nnewcommandunormaluvecn\n\nrenewcommanddmathrmdusing ViscousFlow\nusing Plotswhirl is equipped with a few classes of time marching schemes for advancing time-dependent equations."
},

{
    "location": "manual/timemarching/#Integrating-factor-systems-1",
    "page": "Time marching",
    "title": "Integrating factor systems",
    "category": "section",
    "text": "Integrating factor systems that we encounter in whirl are of the formddt u = A u + r_1(ut) quad u(0) = u_0The operator A may be a matrix or a scalar, but is generally independent of time. (The   method of integrating factors can deal with time-dependent A, but we don\'t encounter   such systems in the whirl context so we won\'t discuss them.) For this purpose, we use the IFRK class of solver, which stands for Integrating Factor Runge-Kutta. This method solves   the part associated with A exactly, via the integrating factor, and advances a modified   equation by Runge-Kutta method to account for the remaining part f.We discussed the construction   of the integrating factor in the context of fields in Fields. But first, let\'s   give an example of how we can solve a simpler problem with just a single scalar-valued   u. The example we will solve isddt u = -alpha u + cos(omega t)quad u(0) = u_0The exact solution is easily obtained:u(t) = u_0 e^-alpha t + frac1alpha^2+omega^2 left alpha(cos(omega t) - e^-alpha t) + omega sin (omega t)rightLet\'s solve it numerically, so we can evaluate the accuracy of the solver. We should note that the integrating factor for this system is e^-alpha t.For demonstration, we will set alpha = 1, omega = 4, and u_0 = 1.using ViscousFlow\nusing Plots\npyplot()α = 1; ω = 4; u₀ = 1.0;Here is the exact solution for later comparisonuex(t) = u₀*exp(-α*t) + (α*(cos(ω*t)-exp(-α*t))+ω*sin(ω*t))/(α^2+ω^2)The first steps are to define operators that provide the integrating factor and the right-hand side of the equations. For the integrating factor, we extend the definition of plan_intfact from Fields.ViscousFlow.plan_intfact(t::Float64,u::Vector{Float64}) = exp(-α*t);Note that we have defined this extended form of plan_intfact to adhere to the standard form, accepting arguments for time t and the state vector u, even though the state vector isn\'t strictly needed here. The state \'vector\' in this problem is actually only a scalar, of course. But the time marching method does not accept scalar-type states currently, so we will make u a 1-element vector to use the whirl tools.Now let us define the right-hand side function. This function should also adhere to the standard form, which requires the state vector u and the time t as arguments.r₁(u::Vector{Float64},t::Float64) = cos(ω*t);We also need to set the time-step size (001) and the initial condition. For the latter, we set up the state vector as a 1-element vector, as discussed earlier:Δt = 0.01;\nu = [u₀];We can now construct the integrator. We supply a form of the state vector (for use as a template   for pre-allocating space for internal storage variables), the time-step size, and the   definitions of the integrating factor and the right-hand side function:ifrk = IFRK(u,Δt,plan_intfact,r₁,rk=TimeMarching.RK31)We have set the time step size to 001. We have also specified that the Runge-Kutta method to be used is a third-order method, RK31, specially designed for storing as few different versions of the integrating factor as necessary. This is actually the default method, so we could have omitted this keyword argument. There are other choices, as well, such as TimeMarching.Euler for the forward Euler method.Now we can solve the system. The integrator has a simple form, accepting as arguments the current time and state, and returning the updated versions of these at the end of the step. We place this integrator inside of a loop and store the results. (Since u is set up   as a 1-element vector, then we will store only the element of this vector.)uhist = Float64[]; # for storing the solution\nT = 0:Δt:10;\nt = 0.0;\nfor ti in T\n  push!(uhist,u[1]) # storage\n  global t, u = ifrk(t,u) # advancement by one step by the integrator\nendNow we can plot the result and compare it with the exact solution.plot(T,uhist,label=\"numerical\",xlabel=\"t\",ylabel=\"u(t)\")\nplot!(T,uex.(T),label=\"exact soln\")\nsavefig(\"ifrk.svg\"); nothing # hide(Image: )As we can see, the results are nearly indistinguishable."
},

{
    "location": "manual/timemarching/#Constrained-systems-1",
    "page": "Time marching",
    "title": "Constrained systems",
    "category": "section",
    "text": ""
},

{
    "location": "manual/timemarching/#Constrained-integrating-factor-systems-1",
    "page": "Time marching",
    "title": "Constrained integrating factor systems",
    "category": "section",
    "text": "Constrained integrating factor systems that we encounter in whirl are of the formddt u = A u - B_1^T f + r_1(ut) quad B_2 u = r_2(ut) quad u(0) = u_0where f is again the Lagrange multiplier for enforcing the constraints on u. Now, we combine the ideas of the last two sections into a single integrator.Let\'s demonstrate this on the example of heat diffusion from a circular ring whose temperature is held constant. In this case, A is the discrete Laplace operator, L, times the heat diffusivity, r_1 is zero (in the absence of volumetric heating sources), and r_2 is the temperature of the ring. The operators B_1^T and B_2 will be the regularization and interpolation operators between discrete point-wise data on the ring and the field data.The ring will have radius 12 and fixed temperature 1, and the heat diffusivity is 1. (In other words, the problem has been non-dimensionalized by the diameter of the circle, the dimensional ring temperature, and the dimensional diffusivity.)First, we will construct a field to accept the temperature onnx = 129; ny = 129; Lx = 2.0; Δx = Lx/(nx-2);\nu₀ = Nodes(Dual,(nx,ny)); # field initial conditionNow set up a ring of points on the circle at center (11).n = 128; θ = range(0,stop=2π,length=n+1);\nR = 0.5; xb = 1.0 .+ R*cos.(θ); yb = 1.0 .+ R*sin.(θ);\nX = VectorData(xb[1:n],yb[1:n]);\nf = ScalarData(X); # to be used as the Lagrange multiplierFrom this, construct the regularization and interpolation operators in their usual symmetric form, and then set up a routine that will provide these operators inside the integrator:reg = Regularize(X,Δx;issymmetric=true)\nHmat, Emat = RegularizationMatrix(reg,f,u₀);\nplan_constraints(u::Nodes{Dual,nx,ny},t::Float64) = Hmat, EmatNow set up the right-hand side operators. Both must take the standard form, with arguments of the types of u and t. For r_1, we will simply set it to a field of zeros in the same type as u. For r_2, we set the result uniformly to 1.r₁(u::Nodes{T,NX,NY},t::Float64) where {T,NX,NY} = Nodes(T,u); # sets to zeros\nr₂(u::Nodes{T,NX,NY},t::Float64) where {T,NX,NY} = 1.0; # sets uniformly to 1.0We will set the time-step size to a large value (10) for demonstration purposes. The method remains stable for any choice. We also initialize time t and the state u:Δt = 1.0;\nt = 0.0;\nu = deepcopy(u₀);Now we can construct the integrator. We supply examples for the state u and the Lagrange multiplier data f, the time-step size, the constructor for the integrating factor, a tuple of the operators for computing the actions of B_1^T and B_2 on data of type f and u, respectively (which, in this case, are matrices Hmat and Emat), and a tuple of the right-hand side functions.ifherk = IFHERK(u,f,Δt,plan_intfact,plan_constraints,(r₁,r₂),rk=TimeMarching.Euler)Here we\'ve set the method to forward Euler. The resulting integrator accepts as arguments the current time t and the current state u, and returns the time, state, and Lagrange multiplier data at the end of the time step.Now, let\'s advance the system. We\'ll also time it.@time for i = 1:20\n  global t, u, f = ifherk(t,u)\nendNow let\'s plot itxg, yg = coordinates(u,dx=Δx);\nplot(xg,yg,u)\nplot!(xb,yb,linecolor=:black,linewidth=1.5)\nsavefig(\"ifherk.svg\"); nothing # hide(Image: )From a side view, we can see that it enforces the boundary condition:plot(xg,u[65,:],xlabel=\"x\",ylabel=\"u(x,1)\")\nsavefig(\"ifherk-side.svg\"); nothing # hide(Image: )"
},

{
    "location": "manual/timemarching/#ViscousFlow.TimeMarching.IFHERK",
    "page": "Time marching",
    "title": "ViscousFlow.TimeMarching.IFHERK",
    "category": "type",
    "text": "IFHERK(u,f,Δt,plan_intfact,B₁ᵀ,B₂,r₁,r₂;[tol=1e-3],[issymmetric=false],[rk::RKParams=RK31])\n\nConstruct an integrator to advance a system of the form\n\ndu/dt - Au = -B₁ᵀf + r₁(u,t) B₂u = r₂(u,t)\n\nThe resulting integrator will advance the system (u,f) by one time step, Δt. The optional argument tol sets the tolerance of iterative saddle-point solution, if applicable.\n\nArguments\n\nu : example of state vector data\nf : example of constraint force vector data\nΔt : time-step size\nplan_intfact : constructor to set up integrating factor operator for A that             will act on type u (by left multiplication) and return same type as u\nplan_constraints : constructor to set up the\nB₁ᵀ : operator acting on type f and returning type u\nB₂ : operator acting on type u and returning type f\nr₁ : operator acting on type u and t and returning u\nr₂ : operator acting on type u and t and returning type f\n\n\n\n\n\n"
},

{
    "location": "manual/timemarching/#ViscousFlow.TimeMarching.IFRK",
    "page": "Time marching",
    "title": "ViscousFlow.TimeMarching.IFRK",
    "category": "type",
    "text": "IFRK(u,Δt,plan_intfact,r₁;[rk::RKParams=RK31])\n\nConstruct an integrator to advance a system of the form\n\ndu/dt - Au = r₁(u,t)\n\nThe resulting integrator will advance the state u by one time step, Δt.\n\nArguments\n\nu : example of state vector data\nΔt : time-step size\nplan_intfact : constructor to set up integrating factor operator for A that             will act on type u (by left multiplication) and return same type as u\nr₁ : operator acting on type u and t and returning u\n\n\n\n\n\n"
},

{
    "location": "manual/timemarching/#ViscousFlow.TimeMarching.RK",
    "page": "Time marching",
    "title": "ViscousFlow.TimeMarching.RK",
    "category": "type",
    "text": "RK(u,Δt,r₁;[rk::RKParams=RK31])\n\nConstruct an integrator to advance a system of the form\n\ndu/dt = r₁(u,t)\n\nThe resulting integrator will advance the state u by one time step, Δt.\n\nArguments\n\nu : example of state vector data\nΔt : time-step size\nr₁ : operator acting on type u and t and returning u\n\n\n\n\n\n"
},

{
    "location": "manual/timemarching/#ViscousFlow.TimeMarching.System",
    "page": "Time marching",
    "title": "ViscousFlow.TimeMarching.System",
    "category": "type",
    "text": "Abstract type for a system of ODEs\n\n\n\n\n\n"
},

{
    "location": "manual/timemarching/#Methods-1",
    "page": "Time marching",
    "title": "Methods",
    "category": "section",
    "text": "Modules = [TimeMarching]\nOrder   = [:type, :function]"
},

{
    "location": "manual/timemarching/#Index-1",
    "page": "Time marching",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"timemarching.md\"]"
},

]}
