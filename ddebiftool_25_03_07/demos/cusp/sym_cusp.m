function varargout=sym_cusp(action,varargin)
%% Automatically generated with matlabFunction
% 
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'ntau'
   varargout{1}=1;
   return
  case 'npar'
   varargout{1}=6;
   return
  case 'nf'
   varargout{1}=2;
   return
  case 'nx'
   varargout{1}=2;
   return
  case 'tp_del'
   varargout{1}=0;
   return
  case 'maxorder'
   varargout{1}=5;
   return
  case 'iscollected'
   varargout{1}=0;
   return
  case 'directional_derivative'
   varargout{1}=1;
   return
  case 'sys_tau_seq'
   varargout{1}={};
   return
  case 'sys_tau_in'
   varargout{1}=[];
   return
  case 'xpattern'
   varargout{1}=[1  1  2  1  1  2  2  2];
   return
end
ind=varargin{1};
order=varargin{2};
nout=varargin{3};
f=str2func(sprintf('sym_cusp_%s_%d_%d',action,ind,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{4:end});
end



function [out_1,out_2] = sym_cusp_rhs_1_0(in1,in2,in3,in4)
%SYM_CUSP_RHS_1_0
%    [OUT_1,OUT_2] = SYM_CUSP_RHS_1_0(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    06-Feb-2025 23:55:23

in_1_n1 = in1(:,1);
in_1_n2 = in1(:,2);
in_1_n3 = in1(:,3);
in_1_n4 = in1(:,4);
in_2_n1 = in2(:,1);
in_2_n2 = in2(:,2);
in_2_n3 = in2(:,3);
in_2_n4 = in2(:,4);
in_2_n5 = in2(:,5);
t2 = in_1_n3.*4.0;
t3 = -t2;
t4 = exp(t3);
t5 = t4+1.0;
t6 = 1.0./t5;
t7 = t6-1.0./2.0;
out_1 = -in_1_n1+in_2_n4-in_1_n4.*in_2_n2+in_2_n1.*t7;
if nargout > 1
    out_2 = -in_1_n2+in_2_n5+in_2_n3.*t7;
end
end


function [out_1,out_2] = sym_cusp_rhs_1_1(in1,in2,in3,in4)
%SYM_CUSP_RHS_1_1
%    [OUT_1,OUT_2] = SYM_CUSP_RHS_1_1(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    06-Feb-2025 23:55:23

in_1_n3 = in1(:,3);
in_1_n4 = in1(:,4);
in_2_n1 = in2(:,1);
in_2_n2 = in2(:,2);
in_2_n3 = in2(:,3);
in_3_n1 = in3(:,1);
in_3_n2 = in3(:,2);
in_3_n3 = in3(:,3);
in_3_n4 = in3(:,4);
in_4_n1 = in4(:,1);
in_4_n2 = in4(:,2);
in_4_n3 = in4(:,3);
in_4_n4 = in4(:,4);
in_4_n5 = in4(:,5);
t2 = in_1_n3.*4.0;
t3 = -t2;
t4 = exp(t3);
t5 = t4+1.0;
t6 = 1.0./t5;
t7 = t6.^2;
t8 = t6-1.0./2.0;
out_1 = -in_3_n1+in_4_n4-in_1_n4.*in_4_n2-in_2_n2.*in_3_n4+in_4_n1.*t8+in_2_n1.*in_3_n3.*t4.*t7.*4.0;
if nargout > 1
    out_2 = -in_3_n2+in_4_n5+in_4_n3.*t8+in_2_n3.*in_3_n3.*t4.*t7.*4.0;
end
end


function [out_1,out_2] = sym_cusp_rhs_1_2(in1,in2,in3,in4)
%SYM_CUSP_RHS_1_2
%    [OUT_1,OUT_2] = SYM_CUSP_RHS_1_2(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    06-Feb-2025 23:55:23

in_1_n3 = in1(:,3);
in_2_n1 = in2(:,1);
in_2_n3 = in2(:,3);
in_3_n3 = in3(:,3);
in_3_n4 = in3(:,4);
in_4_n1 = in4(:,1);
in_4_n2 = in4(:,2);
in_4_n3 = in4(:,3);
t2 = in_1_n3.*4.0;
t3 = in_3_n3.^2;
t5 = in_1_n3.*8.0;
t4 = -t2;
t7 = -t5;
t6 = exp(t4);
t8 = exp(t7);
t9 = t6+1.0;
t10 = 1.0./t9.^2;
t11 = 1.0./t9.^3;
out_1 = in_3_n4.*in_4_n2.*-2.0+in_3_n3.*in_4_n1.*t6.*t10.*8.0-in_2_n1.*t3.*t6.*t10.*1.6e+1+in_2_n1.*t3.*t8.*t11.*3.2e+1;
if nargout > 1
    out_2 = in_3_n3.*in_4_n3.*t6.*t10.*8.0-in_2_n3.*t3.*t6.*t10.*1.6e+1+in_2_n3.*t3.*t8.*t11.*3.2e+1;
end
end


function [out_1,out_2] = sym_cusp_rhs_1_3(in1,in2,in3,in4)
%SYM_CUSP_RHS_1_3
%    [OUT_1,OUT_2] = SYM_CUSP_RHS_1_3(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    06-Feb-2025 23:55:23

in_1_n3 = in1(:,3);
in_2_n1 = in2(:,1);
in_2_n3 = in2(:,3);
in_3_n3 = in3(:,3);
in_4_n1 = in4(:,1);
in_4_n3 = in4(:,3);
t2 = in_1_n3.*4.0;
t3 = in_3_n3.^2;
t4 = in_3_n3.^3;
t6 = in_1_n3.*8.0;
t7 = in_1_n3.*1.2e+1;
t5 = -t2;
t9 = -t6;
t10 = -t7;
t8 = exp(t5);
t11 = exp(t9);
t12 = exp(t10);
t13 = t8+1.0;
t14 = 1.0./t13.^2;
t15 = 1.0./t13.^3;
t16 = t14.^2;
out_1 = in_2_n1.*t4.*t8.*t14.*6.4e+1-in_2_n1.*t4.*t11.*t15.*3.84e+2+in_2_n1.*t4.*t12.*t16.*3.84e+2-in_4_n1.*t3.*t8.*t14.*4.8e+1+in_4_n1.*t3.*t11.*t15.*9.6e+1;
if nargout > 1
    out_2 = in_2_n3.*t4.*t8.*t14.*6.4e+1-in_2_n3.*t4.*t11.*t15.*3.84e+2+in_2_n3.*t4.*t12.*t16.*3.84e+2-in_4_n3.*t3.*t8.*t14.*4.8e+1+in_4_n3.*t3.*t11.*t15.*9.6e+1;
end
end


function [out_1,out_2] = sym_cusp_rhs_1_4(in1,in2,in3,in4)
%SYM_CUSP_RHS_1_4
%    [OUT_1,OUT_2] = SYM_CUSP_RHS_1_4(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    06-Feb-2025 23:55:23

in_1_n3 = in1(:,3);
in_2_n1 = in2(:,1);
in_2_n3 = in2(:,3);
in_3_n3 = in3(:,3);
in_4_n1 = in4(:,1);
in_4_n3 = in4(:,3);
t2 = in_1_n3.*4.0;
t3 = in_3_n3.^3;
t4 = in_3_n3.^4;
t6 = in_1_n3.*8.0;
t7 = in_1_n3.*1.2e+1;
t8 = in_1_n3.*1.6e+1;
t5 = -t2;
t10 = -t6;
t11 = -t7;
t12 = -t8;
t9 = exp(t5);
t13 = exp(t10);
t14 = exp(t11);
t15 = exp(t12);
t16 = t9+1.0;
t17 = 1.0./t16.^2;
t18 = 1.0./t16.^3;
t20 = 1.0./t16.^5;
t19 = t17.^2;
out_1 = in_2_n1.*t4.*t9.*t17.*-2.56e+2+in_2_n1.*t4.*t13.*t18.*3.584e+3-in_2_n1.*t4.*t14.*t19.*9.216e+3+in_2_n1.*t4.*t15.*t20.*6.144e+3+in_4_n1.*t3.*t9.*t17.*2.56e+2-in_4_n1.*t3.*t13.*t18.*1.536e+3+in_4_n1.*t3.*t14.*t19.*1.536e+3;
if nargout > 1
    out_2 = in_2_n3.*t4.*t9.*t17.*-2.56e+2+in_2_n3.*t4.*t13.*t18.*3.584e+3-in_2_n3.*t4.*t14.*t19.*9.216e+3+in_2_n3.*t4.*t15.*t20.*6.144e+3+in_4_n3.*t3.*t9.*t17.*2.56e+2-in_4_n3.*t3.*t13.*t18.*1.536e+3+in_4_n3.*t3.*t14.*t19.*1.536e+3;
end
end


function [out_1,out_2] = sym_cusp_rhs_1_5(in1,in2,in3,in4)
%SYM_CUSP_RHS_1_5
%    [OUT_1,OUT_2] = SYM_CUSP_RHS_1_5(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    06-Feb-2025 23:55:23

in_1_n3 = in1(:,3);
in_2_n1 = in2(:,1);
in_2_n3 = in2(:,3);
in_3_n3 = in3(:,3);
in_4_n1 = in4(:,1);
in_4_n3 = in4(:,3);
t2 = in_1_n3.*4.0;
t3 = in_3_n3.^4;
t4 = in_3_n3.^5;
t6 = in_1_n3.*8.0;
t7 = in_1_n3.*1.2e+1;
t8 = in_1_n3.*1.6e+1;
t9 = in_1_n3.*2.0e+1;
t5 = -t2;
t11 = -t6;
t12 = -t7;
t13 = -t8;
t14 = -t9;
t10 = exp(t5);
t15 = exp(t11);
t16 = exp(t12);
t17 = exp(t13);
t18 = exp(t14);
t19 = t10+1.0;
t20 = 1.0./t19.^2;
t21 = 1.0./t19.^3;
t23 = 1.0./t19.^5;
t22 = t20.^2;
t24 = t20.^3;
out_1 = in_2_n1.*t4.*t10.*t20.*1.024e+3-in_2_n1.*t4.*t15.*t21.*3.072e+4+in_2_n1.*t4.*t16.*t22.*1.536e+5-in_2_n1.*t4.*t17.*t23.*2.4576e+5+in_2_n1.*t4.*t18.*t24.*1.2288e+5-in_4_n1.*t3.*t10.*t20.*1.28e+3+in_4_n1.*t3.*t15.*t21.*1.792e+4-in_4_n1.*t3.*t16.*t22.*4.608e+4+in_4_n1.*t3.*t17.*t23.*3.072e+4;
if nargout > 1
    out_2 = in_2_n3.*t4.*t10.*t20.*1.024e+3-in_2_n3.*t4.*t15.*t21.*3.072e+4+in_2_n3.*t4.*t16.*t22.*1.536e+5-in_2_n3.*t4.*t17.*t23.*2.4576e+5+in_2_n3.*t4.*t18.*t24.*1.2288e+5-in_4_n3.*t3.*t10.*t20.*1.28e+3+in_4_n3.*t3.*t15.*t21.*1.792e+4-in_4_n3.*t3.*t16.*t22.*4.608e+4+in_4_n3.*t3.*t17.*t23.*3.072e+4;
end
end

