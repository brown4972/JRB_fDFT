function [x, reshist, total_iters, total_matvecs, it_hist] ...
             = kl(x0, b, atv_in, loptions_in)
% KL v0.05
% prototype of the wrapper for my old linear solver codes
%
% C. T. Kelley, Oct 1, 2014
% This code comes with no guarantee or warranty of any kind.
%
% function [x, reshist, total_iters, total_matvecs, it_hist] ...
%             = kl(x0, b, atv_in, loptions_in)
%
% Input: x0 = initial iterate
%         b = right side of linear system
%    atv_in = function handle matrix-vector product. 
%             The calling sequence for this function should be
%             av = atv_in(x) or
%             av = atv_in(x, matvec_data) if you have data to pass to
%                  the matvec function. In this case, you put that data
%                  in the matvec_data field of the loptions_in structure
%                  with the kl_optset function.
%
%  loptions_in = options data structure for kl. Type "help kl_optset"
%                for the details.
%
% Output: x = solution
%         reshist = vector of residual norms for the iteration
%         total_iters = iteration count
%         total_matvecs = matvec count 
%                This is the same as the iteration count except
%                for BiCGStab. 
%         it_hist = the entire iteration. Useful for movies. Eats up
%                   lots of storage.
%         -->>          WORKS ONLY FOR PCG in this version.
%
if nargin == 3
   loptions=kl_optset;
else
   loptions=loptions_in;
end
if nargout == 5
   it_hist_flag = 1;
else
   it_hist_flag = 0;
end
atv = atv_in;
mv_struct=struct('loptions',loptions,'atv',atv);
ptv=loptions.ptv;
lmethod=loptions.lmethod;
p_flag = loptions.p_flag;
p_side = loptions.p_side;
p_data = loptions.p_data;
p_data_flag = loptions.p_data_flag;
if strcmp(lmethod,'cg')
   p_side='none';
else
   it_hist=[];
end
if p_flag == 1 && strcmp(p_side,'left')
   bt = mat_vec_0(b,p_data,ptv,p_data_flag);
else
   bt = b;
end
switch lmethod
   case 'gmres'
    goptions=struct('ltol',loptions.ltol,'maxitl',loptions.maxitl,...
           'orthog',loptions.orthog,'reorth',loptions.reorth,...
           'mdata',mv_struct);
    [xt, reshist, total_iters] = gmres_core(x0, bt, @mat_vec, goptions);
   case 'bicgstab'    
    boptions=struct('ltol',loptions.ltol,'maxitl',loptions.maxitl,...
           'mdata',mv_struct);
    [xt, reshist, total_iters] = bicgstab_core(x0, bt, @mat_vec, boptions);
   case 'tfqmr'    
    toptions=struct('ltol',loptions.ltol,'maxitl',loptions.maxitl,...
           'mdata',mv_struct);
    [xt, reshist, total_iters] = tfqmr_core(x0, bt, @mat_vec, toptions);
   case 'cg'
    toptions=struct('ltol',loptions.ltol,'maxitl',loptions.maxitl,...
           'mdata',mv_struct,'it_hist_flag',it_hist_flag);
    switch p_flag
       case 0
         [xt, reshist, total_iters, it_hist] = ...
                        pcg_core(x0, bt, @mat_vec, toptions);
       case 1
         [xt, reshist, total_iters, it_hist] = ...
                        pcg_core(x0, bt, @mat_vec, toptions, @pc_vec);
       end
otherwise
    error('illegal method for kl');
end
total_matvecs=total_iters;
if p_flag == 1 && strcmp(p_side,'right') 
   x = mat_vec_0(xt,p_data,ptv,p_data_flag);
else
   x = xt;
end

function ptv = pc_vec(w,mv_struct)
%
% Preconditioner-vector product for CG
%
loptions=mv_struct.loptions;
ptv = loptions.ptv;
p_data = loptions.p_data;
p_data_flag = loptions.p_data_flag;
ptv=mat_vec_0(w,p_data,ptv,p_data_flag);

function mtv = mat_vec(w,mv_struct)
atv=mv_struct.atv;
loptions=mv_struct.loptions;
ptv = loptions.ptv;
matvec_data_flag = loptions.matvec_data_flag;
matvec_data = loptions.matvec_data;
ptv=loptions.ptv;
lmethod=loptions.lmethod;
p_flag = loptions.p_flag * (1 - strcmp(lmethod,'cg'));
p_side = loptions.p_side;
p_data = loptions.p_data;
p_data_flag = loptions.p_data_flag;
p_side = loptions.p_side;
switch p_flag
   case 0  
      mtv=mat_vec_0(w,matvec_data,atv,matvec_data_flag);
   case 1
      switch p_side
         case 'left'
            mtv=mat_vec_0(w,matvec_data,atv,matvec_data_flag);
            mtv=mat_vec_0(mtv,p_data,ptv,p_data_flag);
         case 'right'
            mtv=mat_vec_0(w,p_data,ptv,p_data_flag);
            mtv=mat_vec_0(mtv,matvec_data,atv,matvec_data_flag);
      end
end     


function mtv0 = mat_vec_0(w,m_data,mtv,dflag)
switch dflag
   case 0
      mtv0=feval(mtv,w);
   case 1
      mtv0=feval(mtv,w,m_data);
end
function [x, reshist, total_iters] = gmres_core(x0, b, atv, goptions)
% GMRES_CORE linear equation solver
% Implementation following Saad-Schultz
%
% New version for KNL/KL. This is the core code with no preconditioning. 
% Use kl.m to manage preconditioning, restarts, and for a simpler interface.
%
% C. T. Kelley, July 2014
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, reshist, total_iters] = gmres_core(x0, b, atv, goptions)
%
%
% Input: x0 = initial iterate
%        b = right hand side
%        atv, a matrix-vector product routine
%            atv must return Ax when x is input
%            the format for atv is
%               function ax = atv(x, a_data)
%            where a_data is auxiliary data that you'll send to the
%            matrix-vector product. The core code requires this format.
%            
%        goptions = structure containing algorithmic parameters.
%
% Output: x=solution
%         reshist = vector of residual norms for the history of
%            the iteration
%         total_iters = number of iterations
%

%
% initialization
%
n=length(b);
errtol=goptions.ltol;
kmax=goptions.maxitl;
orthog=goptions.orthog;
reorth=goptions.reorth;
mdata=goptions.mdata;
x=x0;
%
%
h=zeros(kmax);
v=zeros(n,kmax);
c=zeros(kmax+1,1);
s=zeros(kmax+1,1);
if norm(x) ~=0
   r = b-feval(atv,x,mdata);
else
   r = b;
end
rho=norm(r);
g=rho*eye(kmax+1,1);
errtol=errtol*norm(b);
reshist=[];
%
% test for termination on entry
%
reshist=[reshist,rho];
total_iters=0;
if(rho < errtol) 
    return
end
%
v(:,1)=r/rho;
beta=rho;
k=0;
%
% GMRES iteration
%
while((rho > errtol) && (k < kmax))
    k=k+1;
    v(:,k+1)=feval(atv,v(:,k),mdata);
    normav=norm(v(:,k+1));
    switch orthog
       case 'mgs'
%
% Modified Gram-Schmidt
%
       for j=1:k
           h(j,k)=v(:,j)'*v(:,k+1);
           v(:,k+1)=v(:,k+1)-h(j,k)*v(:,j);
       end
       h(k+1,k)=norm(v(:,k+1));
       normav2=h(k+1,k);
%
% Reorthogonalize?
%
       if  (reorth == 1 && normav + .001*normav2 == normav) 
          for j=1:k
             hr=v(:,j)'*v(:,k+1);
             h(j,k)=h(j,k)+hr;
             v(:,k+1)=v(:,k+1)-hr*v(:,j);
          end
       end
       case 'cgs'
%
% Classical Gram-Schmidt twice.
% This is far more rococo than it needs to be because matlab accesses
% parts of large arrays very slowly. Simply asking for v(:,1:k) eats
% up a lot of time, so I only do that 1/4 of the time. If you wrote
% this loop like you would in C or fortran, you'd use v(:,1:k) directly
% instead of copying into another matrix. In matlab, the performance 
% penalty for doing it the normal way was severe when I did the experiment
% in late 2011.
%
      qt=zeros(k,1);
      kmod=mod(k,4);
      if kmod == 1
         wt=v(:,1:k);
         qt=wt'*v(:,k+1); 
         del=wt*qt;
         kbase = k;
      else
%         qt(1:k-1)=wt'*v(:,k+1);
         qt(1:kbase)=wt'*v(:,k+1);
         del=wt*qt(1:kbase);
         for l=kbase+1:k
           qt(l)=v(:,l)'*v(:,k+1);
           del=del + v(:,l)*qt(l);
         end
%         qt(k)=v(:,k)'*v(:,k+1);
%         del=wt*qt(1:k-1) + v(:,k)*qt(k);
      end
%      v(:,k+1) = v(:,k+1) - wt*qt;
      v(:,k+1) = v(:,k+1) - del;
      h(1:k,k)=qt;
      if reorth == 1      
         if kmod == 1
           qt=wt'*v(:,k+1); 
           del=wt*qt;
         else
           qt(1:kbase)=wt'*v(:,k+1);
           del=wt*qt(1:kbase);
           for l=kbase+1:k
             qt(l)=v(:,l)'*v(:,k+1);
             del=del + v(:,l)*qt(l);
           end
%         qt(1:k-1)=wt'*v(:,k+1);
%         qt(k)=v(:,k)'*v(:,k+1);
%         del=wt*qt(1:k-1) + v(:,l)*qt(k);
        end
%            qt=wt'*v(:,k+1);
%            v(:,k+1) = v(:,k+1) - wt*qt;
            v(:,k+1) = v(:,k+1) - del;
            h(1:k,k)=qt+h(1:k,k);
      end
      otherwise
          disp('Illegal orthogonalization method in gmres_core.');
      end
   h(k+1,k)=norm(v(:,k+1));
%
%   watch out for happy breakdown 
%
    if(h(k+1,k) ~= 0)
         v(:,k+1)=v(:,k+1)/h(k+1,k);
    end
%
%   Form and store the information for the new Givens rotation
%
    if k > 1
        h(1:k,k)=givapp(c(1:k-1),s(1:k-1),h(1:k,k),k-1);
    end
    nu=norm(h(k:k+1,k));
    if nu~=0
%        c(k)=h(k,k)/nu;
        c(k)=conj(h(k,k)/nu);
        s(k)=-h(k+1,k)/nu;
        h(k,k)=c(k)*h(k,k)-s(k)*h(k+1,k);
        h(k+1,k)=0;
        g(k:k+1)=givapp(c(k),s(k),g(k:k+1),1);
    end
%
% Update the residual norm
%
    rho=abs(g(k+1));
    reshist=[reshist,rho];
end
%
% At this point either k > kmax or rho < errtol.
% It's time to compute x and leave.
%
y=h(1:k,1:k)\g(1:k);
total_iters=k;
x = x0 + v(1:n,1:k)*y;

function vrot=givapp(c,s,vin,k)
%  Apply a sequence of k Givens rotations, used within gmres codes
% 
%  C. T. Kelley, July 10, 1994
%
% This code comes with no guarantee or warranty of any kind.
%
%  function vrot=givapp(c, s, vin, k)
%
vrot=vin;
for i=1:k
    w1=c(i)*vrot(i)-s(i)*vrot(i+1);
%    w2=s(i)*vrot(i)+c(i)*vrot(i+1);
    w2=s(i)*vrot(i)+conj(c(i))*vrot(i+1);
    vrot(i:i+1)=[w1,w2];
end
function [x, reshist, total_iters] = ...
                     bicgstab_core(x0, b, atv, boptions)
% BICGSTAB_CORE
% Bi-CGSTAB solver for linear systems. New version for KNL/KL.
%
% C. T. Kelley, April 7, 2011
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, reshist, total_iters]
%                    = bicgstab_core(x0, b, atv, boptions)
%
%
% Input:        x0=initial iterate
%               b=right hand side
%               atv, a matrix-vector product routine
%			atv must return Ax when x is input
%			the format for atv is
%                       function ax = atv(x)
%               boptions = structure containing algorithmic parameters.
%
% Output:       x=solution
%               reshist = vector of iteration residual norms
%               total_iters = number of iterations
%

%
% initialization
%
n=length(b); 
total_iters=0;
errtol=boptions.ltol;
kmax=boptions.maxitl;
mdata=boptions.mdata;
errtol = errtol*norm(b); 
x=x0;
rho=zeros(kmax+1,1);
%
if norm(x)~=0
    r=b-feval(atv,x,mdata);
else
    r=b;
end
reshist=[];
%
hatr0=r;
k=0; rho(1)=1; alpha=1; omega=1;
v=zeros(n,1); p=zeros(n,1); rho(2)=hatr0'*r;
zeta=norm(r); reshist=[reshist,zeta];
%
% Bi-CGSTAB iteration
%
while((zeta > errtol) & (k < kmax))
    k=k+1;
    if omega==0
       error('Bi-CGSTAB breakdown, omega=0');
    end
    beta=(rho(k+1)/rho(k))*(alpha/omega);
    p=r+beta*(p - omega*v);
    v=feval(atv,p,mdata);
    tau=hatr0'*v;
    if tau==0
        error('Bi-CGSTAB breakdown, tau=0');
    end 
    alpha=rho(k+1)/tau;
    s=r-alpha*v; 
    t=feval(atv,s,mdata);
    tau=t'*t;
    if tau==0
        error('Bi-CGSTAB breakdown, tau=0');
    end
    omega=t'*s/tau; 
    rho(k+2)=-omega*(hatr0'*t);
    x=x+alpha*p+omega*s;
    r=s-omega*t;
    zeta=norm(r);
    total_iters=k;
    reshist=[reshist, zeta];
end
function [x, reshist, total_iters, it_hist] = ...
                     pcg_core(x0, b, atv, cgoptions, pcv)
% PCG_CORE linear solver
% Preconditioned Conjugate-Gradient solver. This s a new version for KNL.
% Use knl.m to manage preconditioning and for a simpler interface.
%
% C. T. Kelley, April 2011
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, reshist, total_iters, it_hist] 
%                    = pcg_core(x0, b, atv, cgoptions, pcv)
%
%
% Input:        x0=initial iterate
%               b=right hand side
%               atv, a matrix-vector product routine
%			atv must return Ax when x is input
%			the format for atv is
%                       function ax = atv(x)
%               cgoptions = structure containing algorithmic parameters.
%               pcv, a routine to apply the preconditioner
%                       if omitted, the identity is used.
%                       The format for pcv is
%                       function px = pcv(x).
%
% Output:       x=solution
%               reshist = vector of iteration residual norms
%               total_iters = number of iterations
%               it_hist (optional) = matrix of all iterations
%			useful for movies
%
% 
%

%
% initialization
%
it_hist=[]; 
total_iters=0;
it_hist_flag=cgoptions.it_hist_flag;
n=length(b); reshist=[]; x=x0;
errtol=cgoptions.ltol;
maxiters=cgoptions.maxitl;
mdata=cgoptions.mdata;
if it_hist_flag == 1; it_hist=[it_hist, x]; end
r=b - feval(atv, x, mdata);
if nargin == 4
    z=r;
else
    z = feval(pcv, r, mdata);
end
rho=z'*r;
tst=norm(r);
terminate=errtol*norm(b);
reshist=[reshist,tst];
it=1;
while((tst > terminate) & (it <= maxiters))
%
%
%
if(it==1) 
	p = z;
else
	beta=rho/rhoold;
	p = z + beta*p;
%
% end if
%
end
w = feval(atv, p, mdata);
alpha=p'*w;
%
% Test here to make sure the linear transformation is positive definite.
% A non-positive value of alpha is a very bad sign.
%
if(alpha <= 0)
    [alpha, rho, it]
    error(' negative curvature ') 
end
alpha=rho/alpha;
x=x+alpha*p;
if it_hist_flag == 1; it_hist=[it_hist, x]; end
r = r - alpha*w;
tst=norm(r);
rhoold=rho;
if nargin == 4
    z=r;
else
    z = feval(pcv, r, mdata);
end
rho=z'*r;
it=it+1;
reshist=[reshist,tst];
%
% end while
%
total_iters=it-1;
end
function [x, reshist, total_iters] = ...
                     tfqmr_core(x0, b, atv, toptions)
% TFQMR_CORE
% TFQMR solver for linear systems. New version for KNL/KL.
%
% C. T. Kelley, April 8, 2011
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, reshist, total_iters]
%                    = tfqmr_core(x0, b, atv, toptions)
%
%
% Input:        x0=initial iterate
%               b=right hand side
%               atv, a matrix-vector product routine
%			atv must return Ax when x is input
%			the format for atv is
%                       function ax = atv(x)
%               toptions = structure containing algorithmic parameters.
%
% Output:       x=solution
%               reshist = vector of iteration residual norms
%               total_iters = number of iterations
%
% 
%

%
% initialization
%
n=length(b); 
errtol=toptions.ltol;
kmax=toptions.maxitl;
mdata=toptions.mdata;
errtol = errtol*norm(b);
x=x0;
%
if norm(x)~=0
    r=b-feval(atv,x,mdata);
else
    r=b;
end
reshist=[];
%
u=zeros(n,2); y=zeros(n,2); w = r; y(:,1) = r; 
k=0; d=zeros(n,1); v=feval(atv,y(:,1),mdata); u(:,1)=v;
theta=0; eta=0; tau=norm(r); reshist=[reshist,tau];
rho=tau*tau;
%
% TFQMR iteration
%
while( k < kmax)
    k=k+1;
    sigma=r'*v; 
%
    if sigma==0
        error('TFQMR breakdown, sigma=0')
    end
%
    alpha=rho/sigma;
%
% 
%
    for j=1:2
%
%   Compute y2 and u2 only if you have to
%
        if j==2 
            y(:,2)=y(:,1)-alpha*v;
            u(:,2)=feval(atv,y(:,2),mdata);
        end
        m=2*k-2+j;
        w=w-alpha*u(:,j);
        d=y(:,j)+(theta*theta*eta/alpha)*d;
        theta=norm(w)/tau; c=1/sqrt(1+theta*theta);
        tau=tau*theta*c; eta=c*c*alpha;
        x=x+eta*d;
%
%   Try to terminate the iteration at each pass through the loop
%
        if tau*sqrt(m+1) <= errtol
            reshist=[reshist, tau];
            total_iters=k;
            return
        end
    end
%
%
%
    if rho==0
        error('TFQMR breakdown, rho=0')
    end
%
    rhon=r'*w; beta=rhon/rho; rho=rhon;
    y(:,1)=w + beta*y(:,2);
    u(:,1)=feval(atv,y(:,1),mdata);
    v=u(:,1)+beta*(u(:,2)+beta*v);
    reshist=[reshist, tau];
    total_iters=k;
end
