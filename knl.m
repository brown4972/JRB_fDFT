function [sol, it_hist, ierr, x_hist] = knl(x,f, nloptions,static_data)
% KNL v0.05 and it'll be v0.05 for a while.
%
% C. T. Kelley, May 2014
%
% This code comes with no guarantee or warranty of any kind.
%
% function [sol, it_hist, ierr, x_hist] = knl(x,f, nloptions,static_data)
%
% Kelley's latest nonlinear solver. This will replace the codes from
% my other books when it's done.
%
% Input:
%       initial iterate = x
%        function = f
%           Calling sequence is either fout = f(x)
%           or fout=f(x, static_data)
%           If f uses precomputed data, add that data to the
%           input arguments of knl as the final (fourth) argument.
%
%        nloptions = options structure. See knl_optset for
%                    documentation.
% optional input:
%        static_data = precomputed data for function evaluation,
%                     jacobian/preconditioner-vector product
% output:
%        sol = solution
%        it_hist(maxit,3) = l2 norms of nonlinear residuals
%            for the iteration, number of function evaluations,
%            and number of steplength reductions
%        ierr = 0 upon successful termination
%        ierr = 1 if after maxit iterations
%             the termination criterion is not satsified
%        ierr = 2 failure in the line search. The iteration
%             is terminated if too many steplength reductions
%             are taken.
%
%    x_hist = matrix of the entire iteration history.
%             The columns are the nonlinear iterates. This
%             is useful for making movies, for example, but
%             can consume way too much storage. This is an
%             OPTIONAL argument. Storage is only allocated
%             if x_hist is in the output argument list.
%
%
if nargin == 4
   nloptions=knl_optset('static_data',static_data,nloptions);
end
%
% Process options and ship to knl_core
%
fdlocal=struct('fdflag',nloptions.static_data_flag,...
               'fddata',nloptions.static_data,...
               'fxdata',nloptions.fx_data,'flocal',f);

if nargout == 4
[sol, it_hist, ierr, x_hist] = knl_core(x,@f_internal, nloptions,fdlocal);
else
[sol, it_hist, ierr] = knl_core(x,@f_internal, nloptions,fdlocal);
end


function [fval, pj_data] = f_internal(x, fdlocal)
% F_INTERNAL
%
% Figure out if the function is/is not using optional arguments, and
% package it for the core code.
%
fdflag=fdlocal.fdflag;
f_data=fdlocal.fddata;
fxdata=fdlocal.fxdata;
fun=fdlocal.flocal;
fevaltype = fdflag + 10*fxdata*(nargout-1);
switch fevaltype
    case 0
       fval = feval(fun,x);
    case 1
       fval = feval(fun,x,f_data);
    case 10
       [fval, pj_data] = feval(fun,x);
    case 11
       [fval, pj_data] = feval(fun,x,f_data);
end 
function [sol, it_hist, ierr, x_hist] = knl_core(x,f, nloptions,static_data) 
% KNL_CORE Newton-Krylov solver, globally convergent 
%        solver for f(x) = 0
%
% This is the old nsoli.m code which has been modified 
% to work inside knl. You should never use this directly.
%
% Inexact-Newton-Armijo iteration
%
% Eisenstat-Walker forcing term
%
% Parabolic line search via three point interpolation.
%
% C. T. Kelley, March 2011
%
% This code comes with no guarantee or warranty of any kind.
%
% function [sol, it_hist, ierr, x_hist] = knl_core(x,f, nloptions,static_data) 
%
% input:
%        initial iterate = x
%        function = f
%        nloptions = options structure. See knl_optset for 
%                    documentation.
% NON-optional input (not optional in the core code)
%        static_data = precomputed data for function evaluation,
%                     jacobian/preconditioner-vector product
% output:
%        sol = solution
%        it_hist(maxit,3) = l2 norms of nonlinear residuals
%            for the iteration, number of function evaluations,
%            and number of steplength reductions
%        ierr = 0 upon successful termination
%        ierr = 1 if after maxit iterations
%             the termination criterion is not satsified
%        ierr = 2 failure in the line search. The iteration
%             is terminated if too many steplength reductions
%             are taken.
%
%    x_hist = matrix of the entire interation history.
%             The columns are the nonlinear iterates. This
%             is useful for making movies, for example, but
%             can consume way too much storage. This is an
%             OPTIONAL argument. Storage is only allocated
%             if x_hist is in the output argument list.
%
%
%
% internal parameters:
%
%       alpha = 1.d-4, parameter to measure sufficient decrease
%
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
%
%
% Set internal parameters.
%
alpha = 1.d-4; sigma0 = .1; sigma1 = .5; gamma = .9;
%
%
%
% Read the options structure
%
atol=nloptions.atol;
debug=nloptions.debug;
rtol=nloptions.rtol;
fx_data=nloptions.fx_data;
maxarm=nloptions.maxarm;
maxitl=nloptions.maxitl;
maxit=nloptions.maxit;
etamax=nloptions.etamax;
nloptions=knl_optset('static_data',static_data,nloptions);
eta=abs(etamax);
it_histx = zeros(maxit+1,3);
ierr = 0; 
x_hist=[];
%
n = length(x); fnrm = 1; itc = 0;
%
% Evaluate f at the initial iterate,and
% compute the stop tolerance.
%
if fx_data == 0
    f0 = feval(f,x,static_data);
else
    [f0, fx_data_value]=feval(f,x,static_data);
end
fnrm = norm(f0); s0=zeros(size(x));
if debug == 1
   disp([itc fnrm 0 0 0])
end
it_histx(itc+1,1) = fnrm; it_histx(itc+1,2) = 1; it_histx(itc+1,3) = 0;
fnrmo = 1;
stop_tol = atol + rtol*fnrm;
outstat(itc+1, :) = [itc fnrm 0 0 0];
%
%
% main iteration loop
%
while(fnrm > stop_tol & itc < maxit)
%
% Keep track of the ratio (rat = fnrm/frnmo)
% of successive residual norms and 
% the iteration counter (itc).
%
    rat = fnrm/fnrmo;
    fnrmo = fnrm; 
    itc = itc+1;
    if fx_data ==  1
       nloptions=knl_optset('fx_data_value',fx_data_value,nloptions);
    end
    [step, errstep, inner_it_count,inner_f_evals] = ...
                     knl_linear_solve(s0,f0,x,f,eta,nloptions);
%
%   The line search starts here.
%
    xold = x;
    lambda = 1; lamm = 1; lamc = lambda; iarm = 0;
    xt = x + lambda*step;
    if fx_data == 0
       ft = feval(f,xt,static_data);
    else
       [ft, fx_data_value]=feval(f,xt,static_data);
    end
    nft = norm(ft); nf0 = norm(f0); ff0 = nf0*nf0; ffc = nft*nft; ffm = nft*nft;
    while nft >= (1 - alpha*lambda) * nf0;
%
%   Apply the three point parabolic model.
%
        if iarm == 0
            lambda = sigma1*lambda; 
        else
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm); 
        end
%
% Update x; keep the books on lambda.
%
        xt = x+lambda*step;
        lamm = lamc;
        lamc = lambda;
%
% Keep the books on the function norms.
%
        if fx_data == 0
           ft = feval(f,xt,static_data);
        else
           [ft, fx_data_value]=feval(f,xt,static_data);
        end
        nft = norm(ft);
        ffm = ffc;
        ffc = nft*nft;
        iarm = iarm+1;
        if iarm > maxarm
            disp(' Armijo failure, too many reductions ');
            ierr = 2;
            disp(outstat)
            it_hist = it_histx(1:itc+1,:);
        if nargout == 4, x_hist = [x_hist,x]; end
            sol = xold;
            return;
        end
    end
    x = xt;
    f0 = ft;
%
%   End of line search.
%
    if nargout == 4, x_hist = [x_hist,x]; end
    fnrm = norm(f0);
    it_histx(itc+1,1) = fnrm; 
%
%   How many function evaluations did this iteration require?
%
    it_histx(itc+1,2) = it_histx(itc,2)+inner_f_evals+iarm+1;
%    if itc == 1, it_histx(itc+1,2) = it_histx(itc+1,2)+1; end;
    it_histx(itc+1,3) = iarm;
%
    rat = fnrm/fnrmo;
%
%   Adjust eta as per Eisenstat-Walker.
%
    if etamax > 0
        etaold = eta;
        etanew = gamma*rat*rat;
        if gamma*etaold*etaold > .1
            etanew = max(etanew,gamma*etaold*etaold);
        end
        eta = min([etanew,etamax]);
        eta = max(eta,.5*stop_tol/fnrm);
    end
%
    outstat(itc+1, :) = [itc fnrm inner_it_count rat iarm];
    if debug == 1
       disp([itc fnrm inner_it_count rat iarm])
    end
%
end
sol = x;
it_hist = it_histx(1:itc+1,:);
%if debug == 1
%    disp(outstat)
%    it_hist = it_histx(1:itc+1,:);
%end
%
% on failure, set the error flag
%
if fnrm > stop_tol, ierr = 1; end
%
%
function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
% Apply three-point safeguarded parabolic model for a line search.
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
%
% input:
%       lambdac = current steplength
%       lambdam = previous steplength
%       ff0 = value of \| F(x_c) \|^2
%       ffc = value of \| F(x_c + \lambdac d) \|^2
%       ffm = value of \| F(x_c + \lambdam d) \|^2
%
% output:
%       lambdap = new value of lambda given parabolic model
%
% internal parameters:
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
%

%
% set internal parameters
%
sigma0 = .1; sigma1 = .5;
%
% compute coefficients of interpolation polynomial
%
% p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
%
% d1 = (lambdac - lambdam)*lambdac*lambdam < 0
%      so if c2 > 0 we have negative curvature and default to
%      lambdap = sigam1 * lambda
%
c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
if c2 >= 0
    lambdap = sigma1*lambdac; return
end
c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0);
lambdap = -c1*.5/c2;
if lambdap < sigma0*lambdac, lambdap = sigma0*lambdac; end
if lambdap > sigma1*lambdac, lambdap = sigma1*lambdac; end
%
function [step, errstep, inner_it_count,inner_f_evals] = ...
                 knl_linear_solve(s0,f0,x,f,eta,nloptions)
% KNL_LINEAR_SOLVE
%
% Sort out the options and call kl to solve for the step.
% Do not modify or hack this function unless you are Tim Kelley.
%
% [step, errstep, inner_it_count,inner_f_evals] = ...
%                 knl_linear_solve(s0,f0,x,f,nloptions)
%
% Read the jtv and ptv options.
%
fx_data=nloptions.fx_data;
fx_data_value=nloptions.fx_data_value;
xjacf=nloptions.x_dependent_jtv_flag;
xpcf=nloptions.x_dependent_ptv_flag;
x_dep_pc=nloptions.x_dependent_ptv;
x_dep_jac=nloptions.x_dependent_jtv;
pj_data=nloptions.pj_data;
pjdataf=nloptions.pj_data_flag;
%
% Sanity check.
%
if pjdataf == 1 && fx_data == 1
   disp('Warning fx_data and pj_data are incompatible.');
   disp('Results can be unpredictable.');
end
%
% This core code uses f_internal. This function has as its
% static_data the structure
% static_data_internal = struct('fdflag',nloptions.static_data_flag,...
%               'fddata',nloptions.static_data,'flocal',f);
%
% Be careful with this. If you are talking to f_internal, then
% fdlocal is the static data to use. If you are talking to your
% own Jac-vec or Preconditioner, then you probably want
% static_data=nloptions.static_data.fddata. I have taken care of
% this for you.
%
static_data_internal=nloptions.static_data;
%
% This is the setup for a finite-difference Jac-vec.
%
jac_data=struct('fx',f0,'x',x,'f',f,'static_data',static_data_internal);
loptions=kl_optset('matvec_data',jac_data,'ltol',eta,nloptions);
%
% Are you doing anything yourself? If so, put the original static_data
% where it's supposed to be.
%
if xjacf == 1 || xpcf == 1
        static_data_orig=nloptions.static_data.fddata;
        pjac_data=struct('fx',f0,'x',x,'static_data',static_data_orig);
        if pjdataf == 1
            pjac_data = feval(pj_data, x, f0, static_data_orig);
        end
        if fx_data == 1
        pjac_data=struct('fx',f0,'x',x,'static_data',static_data_orig,...
              'fx_data_value',fx_data_value);
        end
end
%
if xjacf == 1 
%        pjac_data=struct('fx',f0,'x',x,'static_data',static_data_orig);
%        if pjdatf == 1
%            pjac_data = feval(pj_data, x, f0, static_data_orig);
%        end
        loptions=kl_optset('matvec_data',pjac_data,loptions);
end
%
if xpcf == 1
        loptions=kl_optset('p_data',pjac_data,loptions);
        loptions=kl_optset('ptv',x_dep_pc,loptions);
end
%
% Use the finite-difference Jac-vec or your own.
%
switch xjacf 
    case 0
       [step, errstep, inner_it_count,inner_f_evals] = ...
             kl(s0,-f0,@knl_jacvec,loptions);
    case 1
       [step, errstep, inner_it_count,inner_f_evals] = ...
             kl(s0,-f0,x_dep_jac,loptions);
end
function z = knl_jacvec(w,jacvec_data)
% JACVEC 
% Compute a Jacobian-vector product with finite differences
% This is for use in knl.m.
%
% z = f'(x) w
%
static_data=jacvec_data.static_data;
x = jacvec_data.x;
f = jacvec_data.f;
fx = jacvec_data.fx;
z=dirder(x,w,f,fx,static_data);

function z = dirder(x,w,f,f0,static_data)
% Finite difference directional derivative
% Approximate f'(x) w
%
% C. T. Kelley, November 25, 1993
%
% This code comes with no guarantee or warranty of any kind.
%
% function z = dirder(x,w,f,f0,static_data)
%
% inputs:
%           x, w = point and direction
%           f = function
%           f0 = f(x), in nonlinear iterations
%                f(x) has usually been computed
%                before the call to dirder
%  static_data =  precomputed data from calling routine
%
%
% Hardwired difference increment.
epsnew=1.d-7;
%
n=length(x);
%
% scale the step
%
if norm(w) == 0
    z=zeros(n,1);
return
end
epsnew = epsnew/norm(w);
if norm(x) > 0
    epsnew=epsnew*norm(x);
end
%
% del and f1 could share the same space if storage
% is more important than clarity
%
del=x+epsnew*w;
f1=feval(f,del,static_data);
z = (f1 - f0)/epsnew;

