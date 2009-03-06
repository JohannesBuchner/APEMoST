pro sinmod::calc_model, ind=ind 
; ATTENTION :  IND = actual index + 1 !!!!
; CALCULATES THE SIMPLE SIN MODEL

COMMON littleCommon, x_squared, sigma, fac

;pars = parameter array
pars = *self.params

;x = observed abscissa
x = *self.x_dat


*self.model = pars[0]*sin(2.0D*!DPI*pars[1]*x + pars[2])

end

;-----------------------------------------------------
; CALCULATES the model probability (Gaussian noise)
pro sinmod::calc_prob

COMMON littleCommon

fit = self.beta*total(-1.0D*(*self.model - *self.y_dat)^2.0D/(2.0D*sigma^2.0D))
 
self.prob = fit 

end

;-----------------------------------------------------
; SETS the tempering parameter
pro sinmod::set_beta, b

COMMON littleCommon

self.beta = b 

end

;-----------------------------------------------------
; returns the tempering parameter

function sinmod::get_beta

COMMON littleCommon

return, self.beta

end

;-----------------------------------------------------

pro simplesin, filen

COMMON littleCommon

sigma = 0.5D
struct = {sinmod, beta:0D, INHERITS prob_obj}
n_beta = 12.0
sinmod1 = objarr(n_beta)
n_swap = 30D

bet0 = 0.001D
dbet = (1.0D - bet0)/(n_beta-1)

;calculate the parallel tempering parameters
;and call the 'setup'-routine
if n_beta gt 0 then print,' Initializing parallel tempering ...'
for i=0., n_beta-1 do begin
   print, ' Chain ', i+1, ' - beta = ', 1.0D - i*dbet   
   sinmod1[i] = obj_new('sinmod')
   sinmod1[i]->set_beta, 1.0D - i*dbet
   sinmod1[i]->setup, filen
   if i eq 0 then param_descr = sinmod1[i]->get_par_descr(/all)
endfor

;reset all files that will contain analysis results
for i=0., n_elements(param_descr)-1 do begin
   openw,1,param_descr[i] + '_results.dat'
   close,1
endfor


;initialize the model for all tempering chains
x_squared = sinmod1[0]->get_x()
for i=0., n_beta-1.0D do begin
   sinmod1[i]->calc_model
   sinmod1[i]->calc_prob
endfor

n_p = n_elements(param_descr)

;start the markov chain calibration
for i=0., n_beta-1.0D do begin
   sinmod1[i] = calibrate_markov_chain(sinmod1[i], rat_lim=0.5D, burn_in=10000D, iter_lim=20000D, mul=0.85D, /talkative)
   if i eq 0 then begin
      ;set the starting point for the calibration of all hotter distributions to the 
      ; best parameter values of the (beta=1)-distribution
      for j=1., n_beta-1.0D do begin
         sinmod1[j]->set_param, sinmod1[0]->get_param(/best), /all 
	 sinmod1[j]->calc_model
	 sinmod1[j]->calc_prob
      endfor 	  
   endif
endfor

;some parameters for the program flow
end_MCMC = 0
iter = 0L
;setup the plotting windows
for w=0., n_p do window,w,xsize=300,ysize=250

;start the analysis -> run until end_MCMC != 0
while (end_MCMC eq 0) do begin
   ;perform one markov chain step
   for i=0., n_beta-1.0D do begin
      sinmod1[i] = markov_chain(sinmod1[i])
   endfor
    
   ;increase iter 
   iter++
   
   ;if iter is an integer multiple of 2000 -> show acceptance rate, plot distributions, etc.
   if iter mod 2000 eq 0 then begin
      print,'Iteration ' + strtrim(string(iter,FORMAT='(I)'),2)
      print,'Acceptance rate : ', sinmod1[0]->get_ar(/global)      
      for w=1., n_p do begin
         parhist = sinmod1[0]->get_hist(ind=w, nb=200)
         wset,w
         plot, parhist[*,0],parhist[*,1],xs=1,ys=1, xtitle=param_descr[w-1]
      endfor
      print,sinmod1[1]->get_param(/all)
   endif     
   
   ;check for the possibility of a parallel tempering parameter change
   ; and use the seed of the (beta=1)-chain for the random number generation
   seed = sinmod1[0]->get_seed()
   swap = randomu(seed)
   if swap le 1.0D/n_swap then begin
      a = floor((n_beta-1.0)*randomu(seed))
      b = (a + 1.0D)
      p1 = sinmod1[a]->get_probability()
      p1_beta = sinmod1[a]->get_beta()
      p2 = sinmod1[a]->get_probability()
      p2_beta = sinmod1[b]->get_beta()
      r = p1_beta*p2/p2_beta + p2_beta*p1/p1_beta - p1 - p2
      c = randomu(seed) 
      if r ge alog(c) then begin
         temp1 = sinmod1[a]->get_param(/all)
	 sinmod1[a]->set_param, sinmod1[b]->get_param(/all), /all
         sinmod1[b]->set_param, temp1, /all
	 print,'SWAP...', r, c, ' index :', a, b
      endif else print, 'NOSWAP ...' , r, alog(c), ' index :', a, b
   endif   
   ; don't forget to set the new seed !!!
   sinmod1[0]->set_seed, seed      
   
   ; check if a new best value has been found
   sinmod1[0] = check_best(sinmod1[0], win_index=0)	 
   
   ; add the current values to the parameter record
   sinmod1[0]->add_values, iter_lim=60000.0   
  
endwhile

end



