function check_best, prob_obj, win_index=win_index, logplot=logplot
; check if a new best value has been found


prob = prob_obj->get_probability()

if prob_obj->get_best() lt prob then begin
   ;new best value -> set the best value, and plot the new model!
   print,'New best prob. :', prob
   prob_obj->set_best, prob 
   wset,win_index
   if not keyword_set(logplot) then begin
      plot,prob_obj->get_x(),prob_obj->get_data(), xs=1, ys=1
      oplot,prob_obj->get_x(),prob_obj->get_model(), color=255
   endif else begin
      plot, prob_obj->get_x(),prob_obj->get_data(), /ylog, /xlog, xs=1, ys=1
      oplot,prob_obj->get_x(),prob_obj->get_model(), color=255	  
   endelse
   
   ;set the best value - parameters 
   prob_obj->set_param, prob_obj->get_param(/all), /best
endif

return, prob_obj

end
