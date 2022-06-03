function optout = knl_optset(varargin)
% KNL_OPTSET
%
% Set the options in KNL.
% C. T. Kelley, May 2011
%
% function optout = knl_optset(varargin)
%
% THIS IS A PROTOTYPE HACK JOB. 
%
%
% OPTIONS = KNL_OPTSET('name1',val1,'name2',val2, ... );
%           creates the options structure with the options named
%           in the list set to something other than the defaults.
%
% OPTOUT = KNL_OPTSET('name1',val1,'name2',val2, ... , OPTIN);
%           modifies the structure OPTIN and creates a new structure
%           OPTOUT.
%
% I am putting options into this code as needed. If you are one of
% my beta-testers and want something, bother me.
%
% And now for the long list of options. knl uses the linear solver
% kl.m, so some of the options go directly to kl.m. I have grouped
% those under the heading ``Linear Solver Options'', which is below.
%
%                    Nonlinear Solver Algorithmic Options
%
% atol = absolute residual termination criterion; default = 1.d-12
% rtol = relative residual termination criterion; default = 1.d-6;
%        The iteration terminates when 
%                || F(x_n) || <= atol + rtol * ||F(x_0)||
%
% debug = flag for intermediate output, default = 0 (off)
%         debug = 1 prints the histout file as it's built.
%
% etamax = tolerance for the linear iteration. Setting etamax > 0 will
%          use the Eisenstat-Walker method. Setting etamax < 0 means
%          the forcing term will be | etamax | for the entire solve.
%
% maxarm = limit on linesearch reductions. After maxarm reductions without
%          sufficient decrease, knl reports a failure.
%          default=20
%
% maxit  = limit on number of nonlinear iterations
%          default=40
%
%                    Nonlinear solver data options
%
% static_data = structure with any precomputed data you need to send
%              to the evaluation of f, the preconditioner, or the
%              Jacobian. This is indended to eliminate
%              the need for global variables. You can also set this
%              option by using the optional argument to knl. If you 
%              intend to pass this to f, it should be as a second
%              argument, so you'd use f(x, static_data).
%
% fx_data and pj_data only make sense if you are using one or both
% of the x_dependent_jtv or  x_dependent_ptv options
%
% fx_data    = a flag to tell knl that your function can provide
%              preconditioning/analytic Jacobian data as an optional
%              second output argument. If you set this flag to 1, your
%              function will be called with the second argument during
%              the ENTIRE LINE SEARCH to be sure it will be ready 
%              for the linear solve for the inexact Newton direction.
%              You should use this option if the data you need incur
%              a low incremental cost once you've evaluated f. 
%              Otherwise, using the pj_data option might be more efficient.
%
%              If fx_data is on, then your function should look like
%              [fval, fx_data_value] = f(x)  or
%              [fval, fx_data_value] = f(x, static_data) 
%
% WARNING!! fx_data and pj_data are, for the present, incompatible. You
%           may call one or the other, but not both. I intend to fix 
%           this later.
%
% pj_data = a function handle to compute preconditioning and/or analytic
%            Jacobian-vector product data that depends
%            on the nonlinear iteration x. This function has a rigid 
%            calling sequence that knl will use to set up the linear solve.
%            The calling sequence is
%            
%            pj_data = pj_data(pjac_data)
%           
%            where jac_data is a subset of the structure knl uses 
%            for finite-difference Jacobian-vector products:
%
%            pjac_data=struct('fx',fx,'x',x,'static_data',static_data);
%
%            Here, x = current point, fx = f(x), f = handle to f,
%            static_data = your data for the function set with the 
%            static_data option. You must include as much of jac_data as
%            you need within pj_data. 
%
%            The default, if pj_data is not set, is to simply pass 
%            pj_data to your Jacobian-vector and/or preconditioner-vector
%            product.
%
%            You give knl_optset a handle to your preconditioner-vector
%            product function with the ptv option below.
%
% x_dependent_jtv and x_dependent_ptv are handles to Jacobian-vector (jtv)
%            and/or preconditioner-vector (ptv) products which look for
%            x-dependent data from either fx_data, pj_data, or the default
%            structure of 
%
%            struct('fx',fx,'x',x,'static_data',static_data);
%
%            The calling sequences for preconditioner or Jacobian 
%            matrix-vector products must look like
%            ptv=my_ptv(v, my_data)
%
%                    Linear Solver Options
%
% maxitl = maximum number of linear iterations. Default = 40
%
% lmethod = Krylov method. Your choices are 'gmres' (default),
%           'bicgstab', or 'tfqmr'
%
% orthog = orthogonalization methods. Your choices are
%          mgs = modified Gram-Schmidt (default)
%          cgs = classical Gram-Schmidt (orthogonalize twice)
%
%          If you have > 2 cores, you should use cgs for large problems
%          which need many gmres iterations to converge.
%
% ptv = matlab function handle for the preconditioner-vector product. The
%       default is nothing, which means that you have no preconditioner.
%
% p_data = precomputed data you create for the preconditioner. Don't use
%          this if you use a preconditioner which depends on the nonlinear
%          iteration. Your preconditioner set-up function should give me
%          this right before the call to kl. On the other hand, if your
%          preconditioner does not depend on the nonlinear iteration, this
%          is the place to put the data.
%
% p_side = left or right? The default is 'left'.
%
% reorth = 1 (mgs) if you want to test for loss of orthogonality and fix it
%            (the default) or be crazy and do nothing ( = 0).
%            (cgs) reorthogonalize with every linear iteration (good)
%            Never use reorth = 0 and cgs.
%
%
nv=length(varargin);
%
% Even number of arguments means you're setting options from scratch.
%
if 2*floor(.5*nv) == nv
   optin = knl_set_defaults;
   mv=nv/2;
else
%
% Odd number of arguments means you're changing options.
%
   mv=(nv-1)/2;
   optin=varargin{nv};
end
optout=optin;
for i=1:mv
   optout=knl_optset_base(varargin{2*i-1},varargin{2*i},optout);
end

function optout=knl_optset_base(str,val,optin)
%
% This is the internal function to manage the knl options structure.
%
optout=knl_set_defaults;
if nargin == 3
   optout=optin;
end
if nargin > 0
   str=lower(str);
   parms=fieldnames(optout);
   n=length(parms);
   okflag=0;
%
% Test for verbal shortcuts and do the right thing.
%
%   if ischar(val)
%   [okflag,optout]=verbal_shortcut(str,val,optin);
%      if okflag == 0
%         disp('error in knl_optset');
%         disp([val,' is not a legal verbal short cut for parameter ', str]);
%      end
%   else
%
% Update options.
%
       for p=1:n
           if strcmp(str,parms(p))
              optout=setfield(optout,str,val);
              okflag=1;
            end
       end
       if okflag == 0
           disp('error in knl_optset');
           disp([str,' is not a known parameter name']);
       end
%   end
%
% Turn on the flags for the data and preconditioning options.
% Also watch out for thing you are not supposed to set with knl_optset.
%
    if okflag == 1
       switch str
          case 'fx_data'
             if optout.pj_data_flag == 1
                disp('Warning fx_data and pj_data are incompatible');
             end
             optout=setfield(optout,'pj_data_flag',0);
             optout=setfield(optout,'pj_data',[]);
          case 'static_data'
             optout=setfield(optout,'static_data_flag',1);
          case 'ltol'
             disp('Warning. You should not set ltol with knl_optset.');
             disp('The linear solver tolerance is controled with etamax.');
          case 'matvec_data'
             disp('Warning. You should not set the matvec_data option'); 
             disp('with knl_optset. knl does that for you and will.');
             disp('overwrite your settings with its own.');
             optout=setfield(optout,'matvec_data_flag',1);
          case 'ptv'
             optout=setfield(optout,'p_flag',1);
          case 'p_data'
             optout=setfield(optout,'p_data_flag',1);
          case 'pj_data'
             if optout.fx_data == 1
                disp('Warning fx_data and pj_data are incompatible');
             end
             optout=setfield(optout,'fx_data',0);
             optout=setfield(optout,'pj_data_flag',1);
          case 'x_dependent_jtv'
             optout=setfield(optout,'x_dependent_jtv_flag',1);
          case 'x_dependent_ptv'
             optout=setfield(optout,'x_dependent_ptv_flag',1);
       end
    end
end


function knl_defaults = knl_set_defaults
%
% Default options
%
atol = 1.d-12  ;     % termination criteria for absolute nonlinear residual
debug = 0;           % shut up
etamax = .9;         % forcing term for inexact Newton iteration as per
                     % Eisenstat-Walker
fx_data= 0;          % nonlinear function does not speak to preconditioner
                     % or Jacobian
%
fx_data_value = [];  % f will return this as the second argument if fx_data
                     % is turned on
%
static_data = [];     % nonlinear function evaluation needs no data
static_data_flag = 0;
%
matvec_data = [];     % matvec needs no data
matvec_data_flag = 0;
%
maxarm = 20;     % upper limit on steplength reductions
maxit = 40;      % upper limit on nonlinear iterations
maxitl = 40;     % upper limit on linear iterations
%
lmethod = 'gmres'; % use gmres
ltol = 1.d-3;      % linear solver tolerance. DON'T set in knl_optset.
%
pj_data = [];      % No x-dependent preconditioner/Jacobian data.
pj_data_flag = 0;
%
orthog = 'mgs'; % Do the usual thing.
%
ptv = [];        % no preconditioner
p_flag = 0;     
p_side = 'left'; % left is the default
p_data = [];     % preconditioner needs no data
p_data_flag = 0;
%
reorth = 1;    % Test for loss of orthogonality with modified
                 % Gram-Schmidt and reorthogonalize if needed.
%
rtol = 1.d-6;     % termination criteria for relative nonlinear residual
%
x_dependent_ptv = [];  % No x-dependent preconditioning.
x_dependent_ptv_flag = 0;
%                      
x_dependent_jtv = [];  % Take the default finite-difference Jacobian-vector
%                        product, please.
x_dependent_jtv_flag = 0;
%
% Fill the structure
%
knl_defaults=struct(...
'atol', atol,...
'debug', debug,...
'etamax', etamax,...
'fx_data', fx_data,...
'fx_data_value', fx_data_value,...
'static_data', static_data,...
'static_data_flag', static_data_flag,...
'matvec_data', matvec_data,...
'matvec_data_flag', matvec_data_flag,...
'maxarm', maxarm,...
'maxit', maxit,...
'maxitl', maxitl,...
'lmethod', lmethod,...
'orthog',orthog,...
'ptv',ptv,...
'p_flag', p_flag,...
'p_data', p_data,...
'p_side', p_side,...
'p_data_flag', p_data_flag,...
'pj_data', pj_data,...
'pj_data_flag', pj_data_flag,...
'reorth', reorth, ...
'rtol', rtol,...
'ltol',ltol,...
'x_dependent_jtv', x_dependent_jtv,...
'x_dependent_jtv_flag', x_dependent_jtv_flag,...
'x_dependent_ptv', x_dependent_ptv,...
'x_dependent_ptv_flag', x_dependent_ptv_flag);
