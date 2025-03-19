%% Minimalistic DDE-Biftool demo Mackey-Glass Equation 
%
% <http://www.scholarpedia.org/article/Mackey-Glass_equation>
%
% The Mackey-Glass equation is given by
% 
% $$x'(t)=\beta \frac{x(t-\tau)}{1+x(t-\tau)^n}-\gamma x(t)$$
% 
% Parameters are (in this order) |beta|, |n|, |tau| (|gamma| is not part of
% parameter vector).
% See manual for extensive comments.
%% load DDE-Biftool into path
clear
ddebiftool_path([pwd(),'/ddebiftool_25_03_07']);
format compact
clr=lines();
%% Define problem
[ib,in,itau,ig]=deal(1,2,3,4); % define parameter indices
f= @(x,b,n,g)b.*x(2,:)./(1+x(2,:).^n)-g.*x(1,:); % r.h.s., all vectorized for speed-up
df1=@(x,n,dx,dn)1./(x.^n+1).^2.*(dx+dx.*x.^n-dx.*n.*x.^n-dn.*x.*x.^n.*log(x));
df=@(x,b,n,g,dx,db,dn,dg)db.*x(2,:)./(1+x(2,:).^n)+b.*df1(x(2,:),n,dx(2,:),dn)-dg.*x(1,:)-g.*dx(1,:);
%% Convert r.h.s. to DDE-Biftool f(x,p) format,
% note that the parameters are row vectors for historic reasons
funcs=set_funcs('sys_rhs',{1,[ib,in,ig],f},'sys_tau',@()itau,...
    'sys_dirderi',{df},'x_vectorized',true,'p_vectorized',true); % set up problem
%% Set initial guess for state and parameter
% set boundaries, max stepsize, and perform initial correction.
[b,n,tau,g]=deal(2, 10, 0,1);   % initial parameters
x0=(b/g-1)^(1/n);               % initial non-trivial equilibrium
bounds={'max_bound',[itau,2;ib,5],'max_step',[0,0.3]};
[eqbr,suc]=SetupStst(funcs,'x',x0,'parameter',[b,n,tau,g],...
    'step',0.1,'contpar',itau,bounds{:})
%% Compute, find stability and bifurcations of non-trivial equilibria 
disp('Nontrivial equilibria');
figure(1);clf;ax1=gca;xlabel(ax1,'tau');ylabel(ax1,'x (eq)');set(ax1,'FontSize',16);
eqbr=br_contn(funcs,eqbr,10,'ax',ax1);
[eqbr_wbif,nunst_eqs,eqbifs,ihopf]=MonitorChange(funcs,eqbr);
fprintf('Hopf bifurcation near point %d, tau=%g\n',ihopf(1),eqbr_wbif.point(ihopf(1)).parameter(itau));
%% Continue Hopf bifurcation in two parameters b and tau
% starting from point where stability change was detected
% (eqbr.point(ihopf)).
[hopfbr,suc]=SetupHopf(funcs,eqbr_wbif,ihopf(1),...
    'contpar',[ib,itau],'dir',itau,'step',1e-1)
figure(2);clf;ax2=gca; xlabel('b');ylabel('tau');set(ax2,'FontSize',16);
hopfbr=br_contn(funcs,hopfbr,30,'ax',ax2);
hopfbr=br_rvers(hopfbr);
hopfbr=br_contn(funcs,hopfbr,30,'ax',ax2);
[hopfbr,nunst_h]=br_stabl(funcs,hopfbr);
%% Re-plot Hopf bifurcation in 2 parameter plane
figure(3);clf;ax3=gca;
Plot2dBranch(hopfbr,'ax',ax3);
xlabel(ax3,'b');ylabel(ax3,'tau');set(ax3,'FontSize',16);
%% Branch off toward periodic orbits
[per_orb,suc]=SetupPsol(funcs,eqbr_wbif,ihopf,'intervals',20,'degree',4,'max_step',[itau,0.5],...
    'plot_measure',{@(p)p.parameter(itau),@(p)max(p.profile)});
%% Plot initial periodic orbit to illustrate discretization
figure(4);clf;ax4=gca;
plot(ax4,per_orb.point(2).mesh,per_orb.point(2).profile,'o-',...
    per_orb.point(1).mesh,per_orb.point(1).profile,'o-');
xlabel(ax4,'t/period');ylabel(ax4,'x(t)');set(ax4,'FontSize',16);
%% Continuation with "live plot" of periodic orbits
% Floquet multipliers are determined in call to br_stabl. The number of
% unstable Floquet multipliers, nunst_per, indicates that stability is lost
% for increasing tau.
figure(1);
per_orb=br_contn(funcs,per_orb,60,'ax',ax1);
[per_orb_wbif,nunst_per,perbifs,ipd]=MonitorChange(funcs,per_orb);
%% Check that first stability change is a period doubling
per_orb_wbif.point(ipd+1).stability.mu(1:3)
%% plot all profiles, highlighting first period doubling
figure(4);clf;ax4=gca;hold(ax4,'on');
for i=1:length(per_orb.point)
    pt=per_orb_wbif.point(i);
    plot(ax4,pt.mesh,pt.profile,'r-');
end
pt=per_orb_wbif.point(ipd+1);
plot(ax4,pt.mesh,pt.profile,'k-','linewidth',3);
xlabel(ax4,'t/period');ylabel(ax4,'x(t)');set(ax4,'FontSize',16);
%% Find period doubling bifurcations in two parameters
[pdfuncs,pdbr1,suc]=SetupPeriodDoubling(funcs,per_orb_wbif,ipd,...
    'contpar',[ib,itau],'dir',ib,'step',1e-1,bounds{:},'plot_measure',[])
%% Continuation
figure(2);
pdbr1=br_contn(pdfuncs,pdbr1,30,'ax',ax2);
pdbr1=br_rvers(pdbr1);
pdbr1=br_contn(pdfuncs,pdbr1,30,'ax',ax2);
[pdbr1,nunst_pd]=br_stabl(pdfuncs,pdbr1);
%% Branch off at period doubling 
% (Solutions at far end get inaccurate.)
[per2,suc]=DoublePsol(funcs,per_orb_wbif,ipd);
figure(1);
per2=br_contn(funcs,per2,60,'ax',ax1);
[per_orb2_wbif,nunst_p2,perbifs2,ipd2]=MonitorChange(funcs,per2);
%% Plot final one-parameter diagram
% Note that stability boundaries are not accurate as we have not refined
% the branch there. Use MonitorChange or LocateSpecialPoints for
% refinement.)
figure(5);clf;ax5=gca;hold(ax5,'on');
Plot2dBranch({eqbr_wbif,per_orb_wbif,per_orb2_wbif},'ax',ax5);
plot(per_orb_wbif.point(1).parameter(itau),mean(per_orb_wbif.point(1).profile),'ko',...
    'DisplayName','Hopf','MarkerFaceColor','k')
plot(per_orb_wbif.point(ipd).parameter(itau),max(per_orb_wbif.point(ipd).profile),'o','color',clr(1,:),...
    'DisplayName','PD1','MarkerFaceColor',clr(1,:));
plot(per_orb2_wbif.point(ipd2(1)).parameter(itau),max(per_orb2_wbif.point(ipd2(1)).profile),'o','color',clr(4,:),...
    'DisplayName','PD2','MarkerFaceColor',clr(4,:));
legend('Location','southeast')
xlabel(ax5,'tau');ylabel(ax5,'x (eq)');set(ax5,'FontSize',16);
%% Find secondary period doubling and track in two parameters
[pd2funcs,pdbr2,suc]=SetupPeriodDoubling(funcs,per_orb2_wbif,ipd2,...
    'contpar',[ib,itau],'dir',ib,'step',1e-1,'plot_measure',[],bounds{:});
figure(2);
pdbr2=br_contn(pd2funcs,pdbr2,30,'ax',ax2);
pdbr2=br_rvers(pdbr2);
pdbr2=br_contn(pd2funcs,pdbr2,30,'ax',ax2);
[pdbr2,nunst_pd2]=br_stabl(pd2funcs,pdbr2);
%% Two-parameter bifurcation plot cleaned-up
% pass on pdfuncs for period doublings to help identify type of bifurcation
% (without this hint, pdbr1 and pdbr2 look just like periodic orbits)
figure(3);hold(ax3,'on');
Plot2dBranch(pdbr1,'funcs',pd2funcs,'DisplayName','PD1');
Plot2dBranch(pdbr2,'funcs',pd2funcs,'DisplayName','PD2');
%% plot all profiles of period-2 orbits, highlighting period doubling
figure(4);clf;ax4=gca;hold(ax4,'on');
for i=1:length(per_orb2_wbif.point)
    pt=per_orb2_wbif.point(i);
    plot3(ax4,pt.mesh,pt.parameter(itau)+0*pt.mesh,pt.profile,'r-');
end
pt=per_orb2_wbif.point(ipd2);
plot3(ax4,pt.mesh,pt.parameter(itau)+0*pt.mesh,pt.profile,'k-','linewidth',3);
xlabel(ax4,'t/period');ylabel(ax4,'tau');zlabel(ax4,'x(t)');set(ax4,'FontSize',16);
grid(ax4,'on');view(ax4,[5,30])