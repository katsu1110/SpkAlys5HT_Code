function [l_enh, l_sup] = splitList

load('Z:\Corinna\SharedCode\Katsu\frdiff.mat')
load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stim_cond_2007.mat')

l_enh = incl_i(frdiff <= 0);
l_sup = incl_i(frdiff > 0);