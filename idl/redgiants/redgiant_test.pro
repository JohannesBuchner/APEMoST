pro redg::calc_model, ind=ind 

COMMON littleCommon, n_prof, iter, noise, profs, powergauss, height_use, height_lim, jeffrey_mod, life_lim, life_use, noise_lim, noise_use
 
pars = *self.params
x = *self.x_dat

if keyword_set(ind) then begin
   if ind eq 1 then begin
      noise = pars[0]      
   endif
   if (ind gt 1 && ind lt  5) then begin      
      powergauss = pars[1]*exp(-1.0D*(x - pars[2])^2.0D/(2.0D*pars[3]^2.0D))        
   endif
   if (ind gt 4) then begin      
      index = floor((ind - 5)/2.0)
      indexp = (ind - 5) mod 2.0
      if indexp eq 1 then profs[*,index] = pars[ind-2]/(1.0D + (pars[ind-1]*x)^4.0D) $
      else profs[*,index] = pars[ind-1]/(1.0D + (pars[ind]*x)^4.0D)
   endif
endif else begin
   noise = pars[0]
   powergauss = pars[1]*exp(-1.0D*(x - pars[2])^2.0D/(2.0D*pars[3]^2.0D))            
   for k=0., n_prof-1 do begin
      index = 4 + 2*k
      profs[*,k] = pars[index]/(1.0D + (pars[index+1]*x)^4.0D)
   endfor   
endelse

*self.model =  noise + total(profs,2) + powergauss

end

;-----------------------------------------------------

pro redg::calc_prob

COMMON littleCommon
 
 
height_jeffrey =  1.0D/(alog((jeffrey_mod + height_lim)/jeffrey_mod))
;life_jeffrey = 1.0D / (alog((jeffrey_mod + life_lim)/jeffrey_mod))
noise_jeffrey = 1.0D / (alog((jeffrey_mod + noise_lim)/jeffrey_mod))

prior =  mean(alog(height_jeffrey/((*self.params)[height_use] + jeffrey_mod))) ;+ mean(alog(life_jeffrey/((*self.params)[life_use] + jeffrey_mod))) + alog(noise_jeffrey/((*self.params)[noise_use] + jeffrey_mod))
fit = total(-1.0D*(*self.y_dat)/(*self.model) - alog(*self.model)) + prior
 
self.prob = fit 

end

;-----------------------------------------------------

pro redgiant_test, filen

COMMON littleCommon

struct = {redg, INHERITS prob_obj}
giant1 = obj_new('redg')

giant1->setup, filen

param_descr = giant1->get_par_descr(/all)
n_p = n_elements(param_descr)
limits = giant1->get_minmax(/all)
height_use = where(strmid(param_descr,0,1) eq 'a')
height_lim = limits[height_use,1]
life_use = where(strmid(param_descr,0,1) eq 'b')
life_lim = limits[life_use,1]
noise_use = 0
noise_lim = limits[0,1]

;initialize the individual model contributions
noise = 0D
profs = 0D
powergauss = 0D

;initialize prior parameters
jeffrey_mod = 0.0001D

;reset all files that will contain analysis results
for i=0., n_elements(param_descr)-1 do begin
   openw,1,param_descr[i] + '_results.dat'
   close,1
endfor

;calculate the number of profiles: there are 2 parameters per profile
n_prof = (giant1->get_n_par() - 4.0D)/2.0D 
print,n_prof

n_x = n_elements(giant1->get_x())
profs = dblarr(n_x, n_prof)

;initialize the model
giant1->calc_model

giant1 = calibrate_markov_chain(giant1, rat_lim=0.55D, burn_in=3000D, iter_lim=6000D, mul=0.85D, /talkative, /logplot, adjust_step=0.01D)
end_MCMC = 0
iter = 0L

for w=1., n_p do window,w,xsize=300,ysize=250

while (end_MCMC eq 0) do begin
   giant1 = markov_chain(giant1)
   giant1 = check_best(giant1, win_index=5, /logplot)	          
   giant1->add_values, iter_lim=60000.0
   iter++
   if iter mod 2000 eq 0 then begin
      print,'Iteration ' + strtrim(string(iter,FORMAT='(I)'),2)
      print,'Acceptance rate : ', giant1->get_ar(/global)
      for w=1., n_p do begin
         parhist = giant1->get_hist(ind=w, nb=200)
         wset,w
         plot, parhist[*,0],parhist[*,1],xs=1,ys=1, xtitle=param_descr[w-1]
      endfor	  
   endif      
endwhile

end



