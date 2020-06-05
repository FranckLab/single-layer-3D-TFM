#~~~~~~~~~~~~~ Single-layer TPT-based Traction Force Microscopy ~~~~~~~~~~~~~~~~
# Run script using the FEniCS package in Python3 to compute displacment and stress
# from a single top layer of beads used to track 3D displacements
#
# June, 2019; Alex Landauer, Mohak Patel, and Kaushik Vijaykumar
# Franck Lab and Kasari Lab, Brown Univerisity and University of Wisc - Madison
#
#imports
from dolfin import *
from fenics import *
import numpy as np
import scipy as sp
import scipy.io as spio

#function definition in which the entire solution process is conducted
def sl_tfm_solve(data_name_in,data_name_out,load_steps,E,nu,thickness,x_center_norm,y_center_norm):
    print("Running FEniCS...")

    tag = 'quad7'

    #---------------------SET OPTIMIZATION PARAMETERS FOR THE COMPILER----------

    tol = 1e-12;
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]["quadrature_degree"] = 7
    #parameters["form_compiler"]["representation"] = 'quadrature'
    parameters["allow_extrapolation"] = True
    ffc_options = {"optimize": True, \
                   "eliminate_zeros": True, \
                   "precompute_basis_const": True, \
                   "precompute_ip_const": True}

    #---------------------SET UP THE BASELINE FIXED MESH------------------------
    print("Meshing...")
    #Generate predefined mesh
    N = 48
    L = 18
    mesh = UnitCubeMesh.create(N, N, L, CellType.Type.hexahedron)

    # refine mesh near top surface and centered on the cell
    mesh.coordinates()[:,2] = mesh.coordinates()[:,2]**0.60 # refine z

    # refine x and y
    x_dist_ = np.abs(mesh.coordinates()[:,0] - y_center_norm)
    y_dist_ = np.abs(mesh.coordinates()[:,1] - x_center_norm)

    sigma_ = 0.15
    x_gauss_dist = np.exp(-(x_dist_)/(2*sigma_**2))
    y_gauss_dist = np.exp(-(y_dist_)/(2*sigma_**2))

    # x_gauss_dist = (x_dist_)
    # y_gauss_dist = (y_dist_)

    x_all = mesh.coordinates()[:,0] - y_center_norm
    y_all = mesh.coordinates()[:,1] - x_center_norm

    x_shifted_ = x_all*(1.0-x_gauss_dist)
    y_shifted_ = y_all*(1.0-y_gauss_dist)

    x_shifted = np.zeros(x_shifted_.shape)
    y_shifted = np.zeros(y_shifted_.shape)

    right_side = x_shifted_[x_shifted_>=0.0]
    x_right_side = (1.0-y_center_norm)*right_side/(np.max(np.abs(right_side)))
    left_side = x_shifted_[x_shifted_<=0.0]
    x_left_side = x_center_norm*left_side/(np.max(np.abs(left_side)))
    x_shifted[x_shifted_>=0.0] = x_right_side
    x_shifted[x_shifted_<=0.0] = x_left_side

    top_part = y_shifted_[y_shifted_>=0.0]
    y_top_part = (1-y_center_norm)*top_part/(np.max(np.abs(top_part)))
    bottom_part = y_shifted_[y_shifted_<=0]
    y_bottom_part = y_center_norm*bottom_part/(np.max(np.abs(bottom_part)))
    y_shifted[y_shifted_>=0] = y_top_part
    y_shifted[y_shifted_<=0] = y_bottom_part

    mesh.coordinates()[:,0] = x_shifted+y_center_norm
    mesh.coordinates()[:,1] = y_shifted+x_center_norm


    # # mesh size is smaller near z=0 and mapped to a [-1;0] domain along z
    # mesh.coordinates()[:, 2] = -mesh.coordinates()[:, 2]**1.8
    # def denser(x,y,z):
    #     return [a*x, a*y, b*z]
    # x_bar, y_bar, z_bar = denser(x, y, z)
    # xyz_bar_coor = np.array([x_bar, y_bar, z_bar]).transpose()
    # mesh.coordinates()[:] = xyz_bar_coor
    # coor = mesh.coordinates()

    # vtkfile = File("ex_mesh_densified"+tag+".pvd")
    # vtkfile << mesh

    #---------------------IMPORT MATLAB DATA------------------------------------
    print("Import data...")
    #load experimental data (output from the regularized TPT field)
    data = spio.loadmat(data_name_in)
    dispdata = data['dispdata']
    coor = data['coor']
    um2px = data['um2px']


    #---------------------WARP THE MESH TO THE IMAGING DIMENSIONS---------------
    print("Mesh warping...")
    #starting mesh coordinates
    x_ = mesh.coordinates()[:,0]
    y_ = mesh.coordinates()[:,1]
    z_ = mesh.coordinates()[:,2]

    #keep everything in um-scale, so that disps and mesh remain O(1) or O(10)
    #since this is in the microscope. From Matlab, the mesh starts at (0,0,0)
    #and max from the image size (in practice the max location of the beads)

    #warp the mesh to the size of the imagem
    min_x = np.min(x_)
    min_y = np.min(y_)
    min_z = np.min(z_)
    x_ = x_ - min_x
    y_ = y_ - min_y
    z_ = z_ - min_z
    max_x = np.max(coor[:,0])
    max_y = np.max(coor[:,1])

    if thickness > tol:
        max_z = thickness
    else:
        max_uz = np.max(np.abs(dispdata[:,2]))
        max_z = 10*max_uz #for "thick" substrate assumption, make thickness
                          #10x the max z disp

    def rescale(x,y,z):
        return [max_x/np.max(x_)*x, max_y/np.max(y_)*y, max_z/np.max(z_)*z]

    x_bar, y_bar, z_bar = rescale(x_, y_, z_) #rescale to a matching size
    z_bar = z_bar - np.max(z_bar) + np.max(coor[:,2]) #shift to match the top planes
    xyz_bar_coor = np.array([x_bar, y_bar, z_bar]).transpose()
    mesh.coordinates()[:] = xyz_bar_coor
    coor_final = mesh.coordinates()

    #make a function space on the mesh
    V = VectorFunctionSpace(mesh, "Lagrange", 1)

    # vtkfile = File("ex_mesh_warped"+tag+".pvd")
    # vtkfile << mesh

    #--------------------- DEFINE THE TOP SURFACE DISPLACEMENT------------------
    print("Imposing BCs...")
    # Dataset grid
    xi,yi,zi = coor[:,0],coor[:,1],coor[:,2]
    ux,uy,uz = dispdata[:,0],dispdata[:,1],dispdata[:,2]

    #def at x,y,z return only x y -> i.e. turn "slightly" 3D data (from the
    #TPT result) into 2D data so that it can be imposed as a BC
    pts = coor[:,[0,1]]
    fux_ = sp.interpolate.CloughTocher2DInterpolator(pts, ux, fill_value=0)
    fuy_ = sp.interpolate.CloughTocher2DInterpolator(pts, uy, fill_value=0)
    fuz_ = sp.interpolate.CloughTocher2DInterpolator(pts, uz, fill_value=0)

    # Could also use a linear interpolant, but Clough seems to do well
    # fux_ = sp.interpolate.LinearNDInterpolator(pts, dispdata[:,0])
    # fuy_ = sp.interpolate.LinearNDInterpolator(pts, dispdata[:,1])
    # fuz_ = sp.interpolate.LinearNDInterpolator(pts, dispdata[:,2])
    def fux(x,y,z):
        return fux_(x,y)
    def fuy(x,y,z):
        return fuy_(x,y)
    def fuz(x,y,z):
        return fuz_(x,y)

    #define a class to determine the incremental displacment at "time" step t
    class DispIncrementFunc(UserExpression):
        #indented disp
        def __init__(self, t, **kwargs):
            super(DispIncrementFunc, self).__init__(**kwargs)
            self.t = t
        def eval(self, value, x):
            value[0] = fux(x[0],x[1],x[2])*self.t
            value[1] = fuy(x[0],x[1],x[2])*self.t
            value[2] = fuz(x[0],x[1],x[2])*self.t
        def value_shape(self):
            return (3,)

    #define a function space and expression for the imposed displacement
    temp_disp_inc = Function(V)
    temp_disp_inc = DispIncrementFunc(t=0.0, degree=5)


    #---------------------SUBDOMAINS AND BOUNDARY CONDITIONS--------------------

    #set up each boundary subdomain
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1e-10 # tolerance for coordinate comparisons
            return on_boundary and abs(x[0] - 0.0) < tol

    class Right(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1e-10 # tolerance for coordinate comparisons
            return on_boundary and abs(x[1] - 0.0) < tol

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1e-10 # tolerance for coordinate comparisons
            # print(abs(x[2] - np.min(z_bar)))
            return on_boundary and abs(x[2] - np.min(z_bar)) < tol

    class Front(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1e-10 # tolerance for coordinate comparisons
            return on_boundary and abs(x[0] - max_x) < tol

    class Back(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1e-10 # tolerance for coordinate comparisons
            return on_boundary and abs(x[1] - max_y) < tol

    class Top(SubDomain):
        def inside(self, x, on_boundary):
            tol = 1e-10 # tolerance for coordinate comparisons
            return on_boundary and abs(x[2] - np.max(z_bar)) < tol

    # Initialize sub-domain instances
    left   = Left()
    right  = Right()
    top    = Top()
    front = Front()
    back = Back()
    bottom = Bottom()

    # define meshfunction to identify boundaries by numbers
    # boundaries = FacetFunction("size_t", mesh) #"old" FEniCS way of doing things
    boundaries = MeshFunction("size_t", mesh, 2)
    boundaries.set_all(0)
    left.mark(boundaries, 1)   # mark left as 1
    right.mark(boundaries, 2)  # mark right as 2
    bottom.mark(boundaries, 3) # mark bottom as 3
    front.mark(boundaries, 4)   # mark front as 4
    back.mark(boundaries, 5)  # mark back as 5
    top.mark(boundaries, 6)    # mark top as 6

    # Define new measure including boundary naming
    # left: ds(1), right: ds(2), top: ds(3), bottom: ds(4) #"old" way
    dss = ds(subdomain_data=boundaries) #"new" way

    #define BCs (zero disp), can be modified as approperiate
    bcleft = DirichletBC(V.sub(0), Constant(0.0), boundaries, 1)
    bcright = DirichletBC(V.sub(1), Constant(0.0), boundaries, 2)
    bcbottom0 = DirichletBC(V.sub(0), Constant(0.0), boundaries, 3)
    bcbottom1 = DirichletBC(V.sub(1), Constant(0.0), boundaries, 3)
    bcbottom2 = DirichletBC(V.sub(2), Constant(0.0), boundaries, 3)
    bcfront = DirichletBC(V.sub(0), Constant(0.0), boundaries, 4)
    bcback = DirichletBC(V.sub(1), Constant(0.0), boundaries, 5)
    bctop = DirichletBC(V, temp_disp_inc, boundaries, 6)

    bcs = [bcleft,bcback,bcbottom0,bcbottom1,bcbottom2,bctop]

    # other reasonable ways to set bcs, depending on config
    # bcs = [bcleft,bcright,bcfront,bcback,bcbottom0,bcbottom1,bcbottom2,bctop]
    # bcs = [bcbottom0,bcbottom1,bcbottom2,bctop]

    #save mesh to check in Paraview (or similar) if needed
    vtkfile = File("ex_mesh_boundaries_N"+str(N)+tag+".pvd")
    vtkfile << boundaries

    # vtkfile = File("ex_u_post_interp"+str(N)+".pvd")
    # vtkfile << bctop

    #---------------------FUNCTION SPACES TO SOLVE ON---------------------------

    du = TrialFunction(V)            # Incremental displacement
    v  = TestFunction(V)             # Test function
    u  = Function(V)                 # Displacement from previous iteration

    #---------------------KINEMATICS--------------------------------------------

    #the standard finite deformation continuum definitions
    d = u.geometric_dimension() # 2D OR 3D
    I = Identity(d)             # Identity tensor
    F = I + grad(u)             # Deformation gradient
    C = F.T*F                   # Right Cauchy-Green tensor
    B = F*F.T

    #invariants of deformation tensors
    Ic = tr(C)
    J  = det(F)


    #---------------------ELASTICITY AND STRAIN ENERGY--------------------------

    # Elasticity parameters in Lame moduli
    mu, lmbda = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))

    #stress from the compressible Neo-hookean model
    stress = (mu/pow(J,5/3))*(B - (1./3.)*tr(B)*I) + (lmbda + (2.*mu/3.))*(J - 1.)*I

    #stored strain energy density (compressible neo-Hookean model)
    psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

    #total potential energy - no body forces or fancy BCs (add as needed)
    Pi = psi*dx

    # first variation of Pi (directional derivative about u in direction of v)
    F = derivative(Pi, u, v)

    #compute Jacobian of F
    J = derivative(F, u, du)


    #---------------------SET UP THE SOLVERS FOR THE VARIATIONAL PROBLEM--------

    print("Solving for u...")
    # basis vectors and traction, to be used later in the script
    # (old - traction computed in Matlab now)
    #e1 = Constant([1.,0.,0.])
    #e2 = Constant([0.,1.,0.])
    #e3 = Constant([0.,0.,1.])
    #trac1 = stress*e3

    # Non linear variantional problem
    problem_u_nl = NonlinearVariationalProblem(F, u, bcs, J)

    # Non-linear solver
    solver_u = NonlinearVariationalSolver(problem_u_nl)

    # Solver parameters - these seem okay, possible to optimize more though
    prm = solver_u.parameters
    prm['newton_solver']['linear_solver'] = 'bicgstab'
    prm['newton_solver']['preconditioner'] = 'petsc_amg'
    prm['newton_solver']['absolute_tolerance'] = 1E-5
    prm['newton_solver']['relative_tolerance'] = 1E-5
    prm['newton_solver']['maximum_iterations'] = 500
    prm['newton_solver']['relaxation_parameter'] = 1.0


    #---------------------SOLVE THE VARIATIONAL PROBLEM-------------------------

    #loading
    ut = 1 # reference value for the loading (imposed displacement)
    load_min = 0 # load multiplier min value
    load_max = 1 # load multiplier max value
    load_multipliers = np.linspace(load_min,load_max,load_steps)

    #set up solve variables
    var = 0;

    #loop through each "time" step (i.e. increment of displacement) and compute
    #the solve
    for (i_t, t) in enumerate(load_multipliers):
        print(t)
        temp_disp_inc.t = t*ut
        # Solve variational problem
        solver_u.solve()
        var = var + 1

    # save the vtk file for displacement only, helpful for debugging
    #vtkfile = File(data_name_out+"_u_"+tag+".pvd")
    #vtkfile << u

    #---------------------SOLVE FOR STRESSES (AND TRACTIONS)--------------------

    print("Projecting stresses...")

    #set up tensor space and project stresses onto it
    V_uT = TensorFunctionSpace(mesh,"Lagrange",1)
    bcs_stress = [bcleft,bcback,bcbottom0,bcbottom1,bcbottom2] # leave out "top"
    sigma_c = project(stress,V_uT,bcs_stress,solver_type="cg",preconditioner_type="petsc_amg")

    # vtk of stress for debugging
    # print("testing...")
    # vtkfile = File(data_name_out+"_stress_bcs"+tag+".pvd")
    # vtkfile << sigma_c

    #---------------------SAVE DATA TO VTK AND MAT FILES------------------------

    print("Saving output files...")

    #Displacement
    Npoints = np.shape(coor_final)
    vertex_values_u = u.compute_vertex_values(mesh)
    U = np.reshape(vertex_values_u,[Npoints[0],3],'F')

    # vertex_val_T = surface_traction_top.compute_vertex_values(mesh)
    # T = np.reshape(vertex_val_T,[Npoints[0],3],'F')

    #save each component of mesh nodes and u to a .mat variable
    idx_dict = {'x':coor_final[:,0],'y':coor_final[:,1],'z':coor_final[:,2]}
    u_dict = {'u1':U[:,0],'u2':U[:,1],'u3':U[:,2]}
    spio.savemat(data_name_out+"_pts_fenics",idx_dict)
    spio.savemat(data_name_out+"_u_fenics",u_dict)
    #collect displacements from each point in the volume defined by mesh vertex
    #points (i.e., idx_dict)
    s_vals = sigma_c.compute_vertex_values(mesh)
    s_vals = np.reshape(s_vals,[Npoints[0],9],'F')

    #save each component to a .mat variable
    stress_dict = {'s11':s_vals[:,0],'s12':s_vals[:,1],'s13':s_vals[:,2],
                   's21':s_vals[:,3],'s22':s_vals[:,4],'s23':s_vals[:,5],
                   's31':s_vals[:,6],'s32':s_vals[:,7],'s33':s_vals[:,8]}
    spio.savemat(data_name_out+"_stress",stress_dict)

    print("All done!")


#-------------------------END OF SOLUTION FUNCTION------------------------------
