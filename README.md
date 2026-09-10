# Multigrid-FVM

## InputData file
InputData file has the following format:
### First line:
&nfx, &nfy, &ngrid, &ig_obs, &ibeg_obs, &ig_inlt, &ig_outlt, &nvcycle, &nswpm, &nswpp, &isolve, &omega, &convergence
### Second line:
&xl, &yl, &amuf, &prsc, &relxm, &relxp
### Third line:
&ibc_xm, &ibc_xp, &ibc_ym, &ibc_yp
### Fourth line:
&u_xm, &v_xm, &p_xm, &sc_xm
### Fifth through to eighth line:
&u_xp, &v_xp, &p_xp, &sc_xp
&u_ym, &v_ym, &p_ym, &sc_ym
&u_yp, &v_yp, &p_yp, &sc_yp
&u_guess, &v_guess, &p_guess, &sc_guess
### Ninth line (and k more)
&ninlt
&jfinlt[k], &jlinlt[k],&uinlt[k], &vinlt[k], &scinlt[k]
### Tenth line (and k more)
&noutlt
&jfoutlt[k], &jloutlt[k],&uoutlt[k], &voutlt[k], &scoutlt[k]

nfx      : number of cells in x
nfy      : number of cells in y
ngrid    : number of grids for multigriding
ig_obs   :  
ibeg_obs : 
ig_inlt  : 
ig_outlt :
nvcycle  : 
nswpm    : number of sweeps on momentum equation 
nswpp    : number of sweeps on pressure equation
isolve   : 
omega    : 
convergence: Convergence limit, when to stop the simulation 
