function markov_chain, prob_obj, index=index, talkative=talkative, calc_index=calc_index
; calculates a markov chain iteration for any object following the conventions
;
;
; receives: prob_obj ... object for calculation of Markov Chain
;           index  ... optional parameter specifying the index + 1 of a single parameter to change (for calibration)
;           calc_index ... if set, submits index to the model calculation
;           seed ... seed for the random number generator          

prob_old = prob_obj->get_probability()
seed = prob_obj->get_seed()
npar = prob_obj->get_n_par()
accept=0
old_mod = prob_obj->get_model()

if keyword_set(index) then begin
   step = prob_obj->get_steps(ind=index)
   value = prob_obj->get_param(ind=index)
   temp = value + randomn(seed)*step
   lim = prob_obj->get_minmax(ind=index)
   if (temp gt lim[1]) then temp = lim[1] - ((temp - lim[1]) mod (lim[1] - lim[0]))   $
   else if (temp lt lim[0]) then temp = lim[0] + ((lim[0] - temp) mod (lim[1] - lim[0]))
   prob_obj->set_param,temp,ind=index      
endif else begin
   steps = prob_obj->get_steps(/all)
   n_par = n_elements(steps)
   values = prob_obj->get_param(/all)
   temp =  values + randomn(seed, n_par)*steps
   lim = prob_obj->get_minmax(/all)
   use = where(temp lt lim[*,0] OR temp gt lim[*,1])
   if use[0] ne -1 then begin
      use2 = where(temp[use] lt lim[use,0])
      if use2[0] ne -1 then temp[use[use2]] = lim[use[use2],0] + ((lim[use[use2],0] - temp[use[use2]]) mod (lim[use[use2],1] - lim[use[use2],0]))
      use2 = where(temp[use] gt lim[use,1])
      if use2[0] ne -1 then temp[use[use2]] = lim[use[use2],1] - ((temp[use[use2]] - lim[use[use2],1]) mod (lim[use[use2],1] - lim[use[use2],0]))
   endif
   prob_obj->set_param, temp, /all 
endelse

if keyword_set(calc_index) then prob_obj->calc_model,ind=index $
else prob_obj->calc_model

prob_obj->calc_prob
prob_new = prob_obj->get_probability()

if prob_new gt prob_old then accept=1 $
else begin
   test = randomu(seed)   
   if alog(test) lt (prob_new-prob_old) then accept=1 $
   else begin
      prob_obj->set_model, old_mod
      prob_obj->set_probability, prob_old
      if keyword_set(index) then begin
	 prob_obj->set_param, value, ind=index
      endif else begin
         prob_obj->set_param, values, /all
      endelse      
   endelse      
endelse

prob_obj->set_seed, seed

if keyword_set(index) then begin
   prob_obj->set_ar,accept, ind=index   
endif else begin
   prob_obj->set_ar,accept, /all      
endelse

return,prob_obj

end

;----------------------------------------------------
