function T1D_dic=cal_T1d_dictionary(Rdosl, B1, dictionary)

all_inhom_b1=dictionary.all_inhom_b1;
all_jtmt=dictionary.all_jtmt;
var_T1d=dictionary.var_T1d;

% match B1
b1_check = abs(B1 - all_inhom_b1);
num_b1 = find(b1_check == min(b1_check));
% match T1D
if (Rdosl > 0) && (Rdosl < 5)
    yy_check = abs(Rdosl - all_jtmt(:,num_b1));
    num_T1d = find(yy_check == min(yy_check));
    T1D_dic = var_T1d(num_T1d);
end

end