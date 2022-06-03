function optout = kl_optset(varargin)
% KL_OPTSET
%
% Set the options in KL.
% C. T. Kelley, October 1, 2014
%
% THIS IS A PROTOTYPE HACK JOB. 
%
%
% OPTIONS = KL_OPTSET('name1',val1,'name2',val2, ... );
%           creates the options structure with the options named
%           in the list set to something other than the defaults.
%
% OPTOUT = KL_OPTSET('name1',val1,'name2',val2, ... , OPTIN);
%           modifies the structure OPTIN and creates a new structure
%           OPTOUT.
%
% The options are:
%
% matvec_data = data structure you precompute for the matrix-vector
%               product. The default is nothing.
%
% maxitl = maximum number of linear iterations. Default = 40
%
% lmethod = Krylov method. Your choices are 'gmres' (default), 'cg',
%           'bicgstab', or 'tfqmr'
%
% ltol = termination criterion. Iteration stops when resid <= ltol*resid_0
%        The default is 1.d-3.
%
% orthog = orthogonalization methods. Modified Gram-Schmidt 'mgs' or
%          classical 'cgs'. 'cgs' will reorthogonalize with every step
%          if reorthog = 1. Default: cgs + orthogonalize twice.
%
% ptv = matlab function handle for the preconditioner-vector product. The
%       default is nothing, assuming you have no preconditioner.
%
% p_data = precomputed data you create for the preconditioner
%
% p_side = left or right? The default is 'left'.
%
% reorth = 1 if you want to test for loss of orthogonality and fix it 
%            or reorthogonalize every cgs step
%            (the default) or be crazy and do nothing ( = 0).
%
% I am putting options into this code as needed. If you are one of
% my beta-testers and want something, bother me.
%
%
nv=length(varargin);
%
% Even number of arguments means you're setting options from scratch.
%
if 2*floor(.5*nv) == nv
   optin = kl_set_defaults;
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
   optout=kl_optset_base(varargin{2*i-1},varargin{2*i},optout);
end

function optout=kl_optset_base(str,val,optin)
%
% This is the internal function to manage the kl options structure.
%
optout=kl_set_defaults;
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
%         disp('error in kl_optset');
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
           disp('error in kl_optset');
           disp([str,' is not a known parameter name']);
       end
%   end
%
% Turn on the flags for the data and preconditioning options.
%
    if okflag == 1
       switch str
          case 'matvec_data'
             optout=setfield(optout,'matvec_data_flag',1);
          case 'ptv'
             optout=setfield(optout,'p_flag',1);
          case 'p_data'
             optout=setfield(optout,'p_data_flag',1);
       end
    end
end


function kl_defaults = kl_set_defaults
%
% Default options
%
matvec_data = [];     % matvec needs no data
matvec_data_flag = 0;
%
maxitl = 40;     % upper limit on linear iterations
%
lmethod = 'gmres'; % use gmres
%
orthog = 'cgs'; % CGS twice is the default.
%
ptv = [];        % no preconditioner
p_flag = 0;     
p_side = 'left'; % left is the default
p_data = [];     % preconditioner needs no data
p_data_flag = 0;
%
reorth = 1;    % Test for loss of orthogonality with modified
               % Gram-Schmidt and reorthogonalize if needed.
               % Orthogonalize twice with CGS (default)
%
ltol = 1.d-3;     % termination criteria for relative linear residual
%
% Fill the structure
%
kl_defaults=struct(...
'matvec_data', matvec_data,...
'matvec_data_flag', matvec_data_flag,...
'maxitl', maxitl,...
'lmethod', lmethod,...
'ltol',ltol,...
'orthog',orthog,...
'ptv',ptv,...
'p_flag', p_flag,...
'p_data', p_data,...
'p_side', p_side,...
'p_data_flag', p_data_flag,...
'reorth', reorth);
