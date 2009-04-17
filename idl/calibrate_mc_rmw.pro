function calibrate_mc_rmw, prob_obj, burn_in=burn_in, iter_lim=iter_lim, talkative=talkative, logplot=logplot, epsilon=epsilon, a=a
; calibrates the markov chain to an acceptance rate between 20 and 30 %
; for any object following the conventions
; prob_obj  :  probability object conforming to prob_obj_DEFINE.pro
; burn_in : number of burn-in iterations
; iter_lim : number of iterations for step width calibration
; mul : factor for adjusting the step width during calibration
; logplot : if set, plots are displayed with logarithmic scale
; rat_lim : if set, gives the average acceptance rates for individual parameters to be achieved
; adjust_step: if set, gives the factor with which to adjust the stepwidths after burn-in


print,'Beginning calibration of MCMC ...'
n_par = prob_obj->get_n_par()
iter = 0L
sub_iter = 0L
cont = 0
best_params = 0D

window,1

print,'Starting burn-in ...'
while iter lt burn_in do begin
    iter++    
    prob_obj = markov_chain(prob_obj)
    if keyword_set(talkative) then begin
       if iter mod 200 eq 0 then print,' Iteration :', iter 
    endif
    prob_obj = check_best(prob_obj, win_index=1, logplot=logplot)  
    if iter eq (burn_in/2.0) then begin
       print,'Re-initializing burn-in ...'
       prob_obj->set_param, prob_obj->get_param(/best), /all
       prob_obj->set_probability, -100000000D
    endif
endwhile     

prob_obj->set_param, prob_obj->get_param(/best), /all
prob_obj->set_ar, -1, /all, /clear
avgstep = prob_obj->get_steps(/all)
avgstep[*] = 0D
finstep = dblarr(n_par)
print,'Calibrating step widths (RMW)...(set cont=1 to abort)'

prob_obj->set_ar,/clear
while (iter le iter_lim) do begin
    iter++
    steps = prob_obj->get_steps(/all)
    if iter mod 200 eq 0 then begin
       print,'Report ---------------------------------------- iter :', iter
       print,'Steps: ', steps
       print,'Last known alpha: ', alpha
       print,'Global acceptance rate : ', prob_obj->get_ar(/global)         
       wset,0
       if not keyword_set(logplot) then begin
         plot,prob_obj->get_x(),prob_obj->get_data(), xs=1, ys=1
         oplot,prob_obj->get_x(),prob_obj->get_model(), color=255
       endif else begin
         plot, prob_obj->get_x(),prob_obj->get_data(), /ylog, /xlog, xs=1, ys=1
         oplot,prob_obj->get_x(),prob_obj->get_model(), color=255	  
       endelse
       wait,0.01D
       prob_obj->set_ar, -1, /all, /clear   
    endif
    
    p_old = prob_obj->get_probability()
    prob_obj = markov_chain(prob_obj)    
    p_last = prob_obj->get_probability(/last)   
    prob_obj = check_best(prob_obj, win_index=0, logplot=logplot)   	  
    seed = prob_obj->get_seed()
    alpha = exp(p_last - p_old)
    if alpha gt 1 then alpha=1.0D
    stepnew = steps + steps*randomu(seed, n_par)/sqrt(iter)*(alpha - 0.234D)
    use = where(stepnew lt epsilon)
    if use[0] ne -1 then stepnew[use] = epsilon
    use = where(stepnew gt a)
    if use[0] ne -1 then stepnew[use] = a
   
    prob_obj->set_steps, stepnew, /all              
    prob_obj->set_seed, seed    
    avgstep = [[avgstep], [stepnew]]
    ;wset,1
    ;if iter gt 0 && iter mod 2000 eq 0 then begin
    ;   dim = size(avgstep)
    ;   plot,iter - 2000 + dindgen(2000),avgstep[0,dim[2]-2000:dim[2]-1],psym=1
    ;endif   
endwhile

for p=0., n_par-1 do finstep[p] = median(avgstep[p,*])

print,'final steps : ', finstep
prob_obj->set_steps, finstep, /all


return, prob_obj

end
