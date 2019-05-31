function y = cal_translation_err(T_estm, T_true)

y= norm(T_estm-T_true)/norm(T_true)*100;

return
