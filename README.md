# Multigrid-FVM

## InputData file
InputData file has the following format:
### First line:
&nfx, &nfy, &ngrid, &ig_obs, &ibeg_obs, &ig_inlt, &ig_outlt, &nvcycle, &nswpm, &nswpp, &isolve, &omega, &convergence
### Second line:
&xl, &yl, &amuf, &prsc, &relxm, &relxp
### Third line:
&ibc_xm, &ibc_xp, &ibc_ym, &ibc_yp
### Fourth through to seventh line:
&u_xm, &v_xm, &p_xm, &sc_xm
&u_xp, &v_xp, &p_xp, &sc_xp
&u_ym, &v_ym, &p_ym, &sc_ym
&u_yp, &v_yp, &p_yp, &sc_yp
### Eighth line:
&u_guess, &v_guess, &p_guess, &sc_guess
### Ninth line (and k more)
&ninlt
&jfinlt[k], &jlinlt[k],&uinlt[k], &vinlt[k], &scinlt[k]
### Tenth line (and k more)
&noutlt
&jfoutlt[k], &jloutlt[k],&uoutlt[k], &voutlt[k], &scoutlt[k]

nfx      : number of cells in x\
nfy      : number of cells in y\
ngrid    : number of grids for multigriding\
ig_obs   : the grid in which the obstacles are defined (usually lowest and same value as ngrid)\
ibeg_obs : index of the beginning cell of obstacles in x based on the lowest grid (Origin for the obstacle grid) \
ig_inlt  : the grid in which the inlets are defined (unless multiple inlets are used set it as zero)\
ig_outlt : the grid in which the outlets are defined (unless multiple inlets are used set it as zero)\
nvcycle  : number of vcycles\
nswpm    : number of sweeps on momentum equation \
nswpp    : number of sweeps on pressure equation\
isolve   : Choice of solver, 1 for Jacobi, and 2 for Gauss-Seidel\
omega    : (UNUSED)\
convergence: Convergence limit, when to stop the simulation \
xl       : non dimensionalised domain length in x\
yl       : non dimensionalised domain length in y\
amuf     : viscosity (1/Re)\
prsc     : Prandtl number for the scalar variable\
relxm    : over/under relaxation for momentum\
relxp    : over/under relaxation for pressure\
ibc_xm   : Type of boundary condition at left boundary\
ibc_xp   : Type of boundary condition at right boundary\
ibc_ym   : Type of boundary condition at bottom boundary\
ibc_yp   : Type of boundary condition at top boundary\

u_xm, v_xm, p_xm, sc_xm : u, v, p and scalar values at the left boundary
u_xp, v_xp, p_xp, sc_xp : u, v, p and scalar values at the right boundary
u_ym, v_ym, p_ym, sc_ym : u, v, p and scalar values at the bottom boundary
u_yp, v_yp, p_yp, sc_yp : u, v, p and scalar values at the top boundary

u_guess, v_guess, p_guess, sc_guess : u, v, p and scalar guess values for the domain (initial condition)

ninlt    : number of inlets
jfinlt[k], jlinlt[k], uinlt[k], vinlt[k], scinlt[k] : y location begining, y location ending, u, v, scalar value at the inlet

noutlt   : number of outlets
jfoutlt[k], jloutlt[k], uoutlt[k], voutlt[k], scoutlt[k] : y location begining, y location ending, u, v, scalar value at the outlet
