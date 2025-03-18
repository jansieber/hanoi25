function varargout=sym_ddae_demo(action,varargin)
%% Automatically generated with matlabFunction
% 
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'ntau'
   varargout{1}=1;
   return
  case 'npar'
   varargout{1}=4;
   return
  case 'nf'
   varargout{1}=4;
   return
  case 'nx'
   varargout{1}=4;
   return
  case 'tp_del'
   varargout{1}=0;
   return
  case 'maxorder'
   varargout{1}=2;
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
   varargout{1}=[1  1  2  1  3  1  4  1  1  2  2  2  3  2  4  2];
   return
end
ind=varargin{1};
order=varargin{2};
nout=varargin{3};
f=str2func(sprintf('sym_ddae_demo_%s_%d_%d',action,ind,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{4:end});
end



function [out_1,out_2,out_3,out_4] = sym_ddae_demo_rhs_1_0(in1,in2,in3,in4)
%SYM_DDAE_DEMO_RHS_1_0
%    [OUT_1,OUT_2,OUT_3,OUT_4] = SYM_DDAE_DEMO_RHS_1_0(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    07-Feb-2025 00:00:59

in_1_n1 = in1(:,1);
in_1_n2 = in1(:,2);
in_1_n3 = in1(:,3);
in_1_n4 = in1(:,4);
in_1_n5 = in1(:,5);
in_1_n6 = in1(:,6);
in_1_n7 = in1(:,7);
in_1_n8 = in1(:,8);
in_2_n1 = in2(:,1);
in_2_n2 = in2(:,2);
in_2_n3 = in2(:,3);
t2 = in_1_n1+in_1_n2;
t3 = -in_1_n2;
t4 = in_1_n1.*t2;
t5 = in_1_n2.*t2;
t6 = in_1_n1+t3;
t7 = -t4;
t8 = in_1_n1.*t6;
t9 = in_1_n2.*t6;
t10 = t5+t8;
t11 = in_2_n1+t7+t9;
out_1 = in_1_n3.*in_2_n3+in_1_n1.*t11+in_1_n2.*t10;
if nargout > 1
    out_2 = in_1_n4.*in_2_n3-in_1_n1.*t10+in_1_n2.*t11;
end
if nargout > 2
    out_3 = -in_1_n3+in_2_n2.*(in_1_n5-in_1_n7);
end
if nargout > 3
    out_4 = -in_1_n4+in_2_n2.*(in_1_n6-in_1_n8);
end
end


function [out_1,out_2,out_3,out_4] = sym_ddae_demo_rhs_1_1(in1,in2,in3,in4)
%SYM_DDAE_DEMO_RHS_1_1
%    [OUT_1,OUT_2,OUT_3,OUT_4] = SYM_DDAE_DEMO_RHS_1_1(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    07-Feb-2025 00:00:59

in_1_n1 = in1(:,1);
in_1_n2 = in1(:,2);
in_1_n3 = in1(:,3);
in_1_n4 = in1(:,4);
in_1_n5 = in1(:,5);
in_1_n6 = in1(:,6);
in_1_n7 = in1(:,7);
in_1_n8 = in1(:,8);
in_2_n1 = in2(:,1);
in_2_n2 = in2(:,2);
in_2_n3 = in2(:,3);
in_3_n1 = in3(:,1);
in_3_n2 = in3(:,2);
in_3_n3 = in3(:,3);
in_3_n4 = in3(:,4);
in_3_n5 = in3(:,5);
in_3_n6 = in3(:,6);
in_3_n7 = in3(:,7);
in_3_n8 = in3(:,8);
in_4_n1 = in4(:,1);
in_4_n2 = in4(:,2);
in_4_n3 = in4(:,3);
t2 = in_1_n1+in_1_n2;
t3 = in_3_n1+in_3_n2;
t4 = -in_1_n2;
t5 = -in_3_n2;
t6 = in_1_n1.*t2;
t7 = in_1_n2.*t2;
t8 = in_3_n1.*t2;
t9 = in_3_n2.*t2;
t10 = in_1_n1.*t3;
t11 = in_1_n2.*t3;
t12 = in_1_n1+t4;
t13 = in_3_n1+t5;
t14 = -t6;
t15 = -t8;
t16 = -t10;
t17 = in_1_n1.*t12;
t18 = in_1_n2.*t12;
t19 = in_3_n1.*t12;
t20 = in_3_n2.*t12;
t21 = in_1_n1.*t13;
t22 = in_1_n2.*t13;
t23 = t7+t17;
t24 = in_2_n1+t14+t18;
t25 = t9+t11+t19+t21;
t26 = in_4_n1+t15+t16+t20+t22;
out_1 = in_1_n3.*in_4_n3+in_2_n3.*in_3_n3+in_1_n1.*t26+in_1_n2.*t25+in_3_n1.*t24+in_3_n2.*t23;
if nargout > 1
    out_2 = in_1_n4.*in_4_n3+in_2_n3.*in_3_n4-in_1_n1.*t25+in_1_n2.*t26-in_3_n1.*t23+in_3_n2.*t24;
end
if nargout > 2
    out_3 = -in_3_n3+in_4_n2.*(in_1_n5-in_1_n7)+in_2_n2.*(in_3_n5-in_3_n7);
end
if nargout > 3
    out_4 = -in_3_n4+in_4_n2.*(in_1_n6-in_1_n8)+in_2_n2.*(in_3_n6-in_3_n8);
end
end


function [out_1,out_2,out_3,out_4] = sym_ddae_demo_rhs_1_2(in1,in2,in3,in4)
%SYM_DDAE_DEMO_RHS_1_2
%    [OUT_1,OUT_2,OUT_3,OUT_4] = SYM_DDAE_DEMO_RHS_1_2(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    07-Feb-2025 00:00:59

in_1_n1 = in1(:,1);
in_1_n2 = in1(:,2);
in_3_n1 = in3(:,1);
in_3_n2 = in3(:,2);
in_3_n3 = in3(:,3);
in_3_n4 = in3(:,4);
in_3_n5 = in3(:,5);
in_3_n6 = in3(:,6);
in_3_n7 = in3(:,7);
in_3_n8 = in3(:,8);
in_4_n1 = in4(:,1);
in_4_n2 = in4(:,2);
in_4_n3 = in4(:,3);
t2 = in_1_n1+in_1_n2;
t3 = in_3_n1+in_3_n2;
t4 = -in_1_n2;
t5 = -in_3_n2;
t6 = in_3_n1.*t2;
t7 = in_3_n2.*t2;
t8 = in_1_n1.*t3;
t9 = in_1_n2.*t3;
t10 = in_1_n1+t4;
t11 = in_3_n1+t5;
t12 = in_3_n1.*t3.*2.0;
t13 = in_3_n2.*t3.*2.0;
t14 = -t6;
t15 = -t8;
t17 = in_3_n1.*t10;
t18 = in_3_n2.*t10;
t19 = in_1_n1.*t11;
t20 = in_1_n2.*t11;
t21 = in_3_n1.*t11.*2.0;
t22 = in_3_n2.*t11.*2.0;
t23 = t13+t21;
t25 = t7+t9+t17+t19;
t26 = in_4_n1+t14+t15+t18+t20;
out_1 = in_3_n3.*in_4_n3.*2.0+in_1_n2.*t23+in_3_n1.*t26.*2.0+in_3_n2.*t25.*2.0-in_1_n1.*(t12-t22);
if nargout > 1
    out_2 = in_3_n4.*in_4_n3.*2.0-in_1_n1.*t23-in_3_n1.*t25.*2.0+in_3_n2.*t26.*2.0+t4.*(t12-t22);
end
if nargout > 2
    out_3 = in_4_n2.*(in_3_n5-in_3_n7).*2.0;
end
if nargout > 3
    out_4 = in_4_n2.*(in_3_n6-in_3_n8).*2.0;
end
end

