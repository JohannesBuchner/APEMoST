function calibrate_markov_chain, prob_obj, burn_in=burn_in, iter_lim=iter_lim, mul=mul, talkative=talkative, logplot=logplot, rat_lim=rat_lim,adjust_step=adjust_step
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

;DEBUG!!!!
ar_old = 0D
st_old = 0D
flag = 0

;DEBUG ENDE !!!

if not keyword_set(mul) then mul = 0.85D
if not keyword_set(adjust_step) then adjust_step = 0.5D
if not keyword_set(rat_lim) then rat_lim = 0.25D^(1.0D/n_par)

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

prob_obj->set_steps, prob_obj->get_steps(/all)*adjust_step, /all
prob_obj->set_param, prob_obj->get_param(/best), /all

print,'Calibrating step widths ...(set cont=1 to abort)'

prob_obj->set_ar,/clear
while (iter le iter_lim) do begin
    iter++
    for i=0., n_par-1 do begin
       prob_obj = markov_chain(prob_obj, index=i+1, /calc_index)      
       prob_obj = check_best(prob_obj, win_index=1, logplot=logplot)  	  
    endfor
    
    if iter mod 200 eq 0 then begin
       steps = prob_obj->get_steps(/all)
       print,steps
       ar = prob_obj->get_ar(/all)
       descr = prob_obj->get_par_descr(/all)
       if flag eq 0 then begin
          flag = 1
	  ar_old = ar
	  st_old = steps
       endif 
       for i=0., n_elements(ar)-1 do begin
          print,'----------------------------------------- iteration : ', iter
          print,'Acceptance rates: '
          print,'Parameter :', descr[i], ' - ', ar[i], ' - changed by : ', (ar[i]-ar_old[i])/ar_old[i]*100D, ' percent'
	  print,'Steps width changed by ', (steps[i]-st_old[i])/st_old[i]*100D,' percent'
          print,'--------------------------------------------------------------'
       endfor	 
       
       st_old = steps
       ar_old = ar 
       
       use = where(ar gt (rat_lim + 0.05D))
       if use[0] ne -1 then steps[use] = steps[use]/mul
       use = where(ar lt (rat_lim - 0.05D))
       if use[0] ne -1 then steps[use] = steps[use]*mul
       print,steps
       
       if not keyword_set(logplot) then begin
         plot,prob_obj->get_x(),prob_obj->get_data(), xs=1, ys=1
         oplot,prob_obj->get_x(),prob_obj->get_model(), color=255
       endif else begin
         plot, prob_obj->get_x(),prob_obj->get_data(), /ylog, /xlog, xs=1, ys=1
         oplot,prob_obj->get_x(),prob_obj->get_model(), color=255	  
       endelse
       wait,0.01D
       prob_obj->set_steps, steps, /all              
       prob_obj->set_param, prob_obj->get_param(/best), /all
       prob_obj->set_probability, -100000000D
       prob_obj->set_ar,/clear             
       
       prob_obj->set_ar, -1, /all, /clear       
       sub_iter = 0L
       while (sub_iter lt 200) do begin
          sub_iter++
          prob_obj = markov_chain(prob_obj)       
          prob_obj = check_best(prob_obj, win_index=1, logplot=logplot)   	  
       endwhile	     
       print,'Global acceptance rate : ', prob_obj->get_ar(/global)     
       deltarat = prob_obj->get_ar(/global) - 0.23D
       if abs(deltarat) lt 0.01 then cont=1 else begin
          if deltarat lt 0 then rat_lim /= 0.99D $
	  else rat_lim *= 0.99D
       endelse          
    endif    
endwhile



return, prob_obj

end
