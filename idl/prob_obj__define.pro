;-------------------------------------------------------------------
;NEEDED FOR USAGE OF :
;
;markov_chain.pro 
;calibrate_markov_chain.pro
;-------------------------------------------------------------------
;
; defines the standard object for any MCMC analysis package:
; prob_obj_DEFINE - defines prob_obj for usage in other programs 
;------------------------------------------------------------------- 
;
; Description of prob_obj:
; 
; x) general class that contains all data and interfaces for MCMC analysis of 1-D data
; x) for implementation new classes need to be created which inherit prob_obj, but overwrite
;    the "calc_model"- as well as the "calc_prob"-procedure and, if necessary, the constructor and destructor  
;
;          -------------------------------------------------------------------
;          data: 


;n_par		    : number of parameters
;accept     	    : number of accepted steps for MCMC (after calibration)
;reject     	    : number of rejected steps for MCMC (after calibration)
;prob	    	    : probability of the most recently evaluated parameter values
;prob_best  	    : probability of best parameter values yets
;seed	    	    : pointer to 1D-array as seed for random number generators 
;params     	    : pointer to 1D-array containing the parameter values
;params_best        : pointer to 1D-array containing the best parameter values yet
;params_distr       : pointer to 2D-array containing the resulting values for each parameter (parameter record array)
;params_distr_buf   ; pointer to 2D-array which acts as a buffer for storing parameter values for quickly (parameter record buffer)
;params_descr	    : pointer to 1D-array containing description of parameters in string format
;params_priors	    : pointer to 1D-array containing the current value of the prior probability for each parameter
;params_ar  	    : pointer to 2D-array containing the number of accepted/rejected steps for individual parameters
;params_step 	    : pointer to 1D-array containing the current step width for individual parameters
;params_minmax 	    : pointer to 2D-array containing the lower and upper limits for each parameter
;x_dat      	    : pointer to 1D-array containing the abscissa of the observations (e.g., date, frequency, wavelength)
;y_dat	    	    : pointer to 1D-array containing the ordinate of the observations (e.g., intensity, power, radial vel.)
;model	    	    : pointer to 1D-array containing the model values corresponding to the observed abscissa values

;
;          -------------------------------------------------------------------
;          methods:
;
; 
; prob_obj::INIT  - constructor of prob_obj     
; prob_obj::CLEANUP - destructor of prob_obj    
; prob_obj::add_values, iter_lim=iter_lim 	
; prob_obj::calc_model, ind=ind, all=all 
; prob_obj::calc_prob 
; prob_obj::write2files
;
;---------------------------------------------------
; prob_obj::setup, filename 
;
; ;FILE FORMAT for input :
;-----------------
; 1.) first line contains the filename of the observations (two column ascii file)
; 2.) second to nth line is a multi-column format separated by white spaces (here denoted by | ) and containing :
;      | in. value | low. limit | up. limit | par. descr. | in. stepwidth |
;
; in. value: initial value for parameter (float or double)
; low. limit : lower boundary for parameter (float or double)
; up. limit : upper boundary for parameter (float or double)
; par. descr : description of parameter (string - should be filename compatible)
; in. stepwidth : stepwidth for each parameter as a fraction of the total range (e.g., 0.1 means  0.1*(up. limit - low.limit)) 
;
;-----------------------------------------------------
;
;    	    	    
; 
; prob_obj::get_ar(ind=ind, all=all, global=global)        
; prob_obj::get_best()	    	    	    
; prob_obj::get_data()	 
; prob_obj::get_hist(ind=ind,nb=nb)   	    	    
; prob_obj::get_minmax(ind=ind, all=all)    
; prob_obj::get_model()     	    	    
; prob_obj::get_n_par()     	    	    
; prob_obj::get_param(ind=ind, all=all, best=best)     
; prob_obj::get_par_descr(ind=ind, all=all) 
; prob_obj::get_probability()	    	    
; prob_obj::get_seed()	    	    	    
; prob_obj::get_steps(ind=ind, all=all)     
; prob_obj::get_x() 	    	    	    
;
; prob_obj::set_ar, new_ar, ind=ind, all=all, clear=clear   	    	
; prob_obj::set_best, new_best                	    	    	    	
; prob_obj::set_data, new_data                	    	    	    	
; prob_obj::set_minmax, new_minmax, ind=ind, all=all        	    	
; prob_obj::set_model, new_model, init=init                 	    	
; prob_obj::set_n_par, new_n_par                            	    	
; prob_obj::set_param, new_param, ind=ind, all=all, init=init, best=best       	
; prob_obj::set_par_descr, new_par_descr, ind=ind, all=all, init=init 	
; prob_obj::set_probability, new_prob
; prob_obj::set_seed, new_seed      	    	    	    	    	
; prob_obj::set_steps, new_steps, ind=ind, all=all, init=init       	
; prob_obj::set_x, new_x                                            	
;
;--------------------------------------------------------------

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::INIT
;constructor for prob_obj

;initialize seed for random number generators
temp = randomu(someseed)
self.seed = ptr_new(someseed)

;set initial values for probabilities (as low as possible)
self.prob = -100000000D
self.prob_best = -100000000D

;set initial values for acceptance/rejection parameters
self.accept = 0L
self.reject = 0L

;set initial values of parameters to 0
self.n_par = 0

return, 1

;!!!!!!!!!!!! OVERWRITE if necessary !!!!!!!!!!!!
;!!!!!! in which case you should call :  !!!!!!!!
;!!!!!!!!! someVar = prob_obj::INIT() !!!!!!!!!!!

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::CLEANUP
;destructor for prob_obj

;deallocate memory from standard pointers
ptr_free, params, params_descr, params_ar, x_dat, y_dat, model, params_step, params_min, params_max

;!!!!!!!!!!!! OVERWRITE if necessary !!!!!!!!!!!!
;!!!!!! in which case you should call :  !!!!!!!!
;!!!!!!!!!!!! prob_obj::CLEANUP() !!!!!!!!!!!!!!!
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::add_values, iter_lim=iter_lim
; adds values to the record of visited parameter configurations during runtime
; limit : if set, gives the maximum number of iterations to be stored before writing to file

; determine the number of elements that have been added to the buffer
buf_size = size(*self.params_distr_buf)

; determine the number of elements that have been added to the actual configuration array
dist_size = size(*self.params_distr)

;if the buffer has more than 200 entries per parameter, and has at least 2 dimensions (n_elements(buf_size) >= 5), 
;append the buffer to the main array and empty the buffer
if buf_size[2] gt 200 && n_elements(buf_size) eq 5 then begin

   ;check whether the main array already contains real values (sum(array) != 0)
   if n_elements(dist_size) eq 5 || total(*self.params_distr) ne 0 then begin
   
      ;if the number of elements in the actual array is larger than iter_lim -> write to files, empty the array and add elements of the buffer
      if dist_size[2] gt iter_lim then begin
      
         ;write to files
         self->write2files
	
	 ;reset the array to contain only the buffer elements 
	 *self.params_distr = *self.params_distr_buf

      ;else append the contents of the buffer to the main array	 
      endif else *self.params_distr = [[(*self.params_distr)], [(*self.params_distr_buf)]]
   endif else begin
   
      ;main array has not yet been initialized -> initialize it with the contents of the buffer      
      *self.params_distr = *self.params_distr_buf
      
   endelse   
   
   ;reset the buffer to contain the current parameter values
   *self.params_distr_buf = *self.params
   
endif else begin

   ;check if the buffer contains some real values (sum(array) != 0)
   if total(*self.params_distr_buf) ne 0 then begin
      
      ;buffer contains real values but is not full -> append current valus to buffer
      *self.params_distr_buf = [[(*self.params_distr_buf)],[(*self.params)]]
      
   endif else begin
   
      ;buffer does not contain real values -> initialize the buffer with the current values
      *self.params_distr_buf = *self.params     
      
   endelse
endelse

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::calc_model, ind=ind, all=all
;calculates the model
;ind : if set, gives the index + 1 of the parameter that was changed prior to recalculating 
;all : if set, signifies that all parameters are assumed to be changed - the whole model should be recalculated


;!!!!!!!!!!!! OVERWRITE !!!!!!!!!!!!!!!

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::calc_prob
;calculates the probability of the model

;!!!!!!!!!!!! OVERWRITE !!!!!!!!!!!!!!!

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::write2files
;writes the analysis-results to the model

;determine the dimensions of the parameter record
dim = size(*self.params_distr)

;print warning for making sure that the user does not disturb the writing proces
print,'ATTENTION - Writing to File - do not interrupt'

;for all parameters ...
for i=0., self.n_par-1 do begin

   ;open/append to a file for the parameter, using the parameter description as a part of the filename 
   
   openw,1,(*self.params_descr)[i] + '_results.dat', /append
   
   ;print the parameter record subarray for the current parameter using a 5-column format
   printf,1, (*self.params_distr)[i,*], FORMAT='(D,D,D,D,D)'
   
   ;close the file
   close,1
endfor

;print 'done'
print,'---> Writing done <---'


end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::setup, filename
;reads the initial setup of the model from a file and sets up the structure of the prob. object
;filename : path to the file containing the necessary information 
;
;FILE FORMAT:
;-----------------
; 1.) first line contains the filename of the observations (two column ascii file)
; 2.) second to nth line is a multi-column format separated by white spaces (here denoted by | ) and containing :
;      | in. value | low. limit | up. limit | par. descr. | in. stepwidth |
;
; in. value: initial value for parameter (float or double)
; low. limit : lower boundary for parameter (float or double)
; up. limit : upper boundary for parameter (float or double)
; par. descr : description of parameter (string - should be filename compatible)
; in. stepwidth : stepwidth for each parameter as a fraction of the total range (e.g., 0.1 means  0.1*(up. limit - low.limit)) 


;initialize local variables ...
error = 0
iofile = ''

;setup error handling for I/O - as long as no problem is encountered error remains 0
;the code will continue to run even if an error is encountered
catch, error


;error handling -> this code will be run if error != 0 anywhere in the procedure
if error ne 0 then 
   ;print an error message   
   print,'File I/O error - there was a problem with ', filename
   
   ;get user input to ensure that the error message was read
   read,dummy
endif

; open the setup file
openr,1,filename

; read first line which contains filename of the observed data file
readf,1,iofile

; close the setup file
close,1

; read the observed data
rdfloat,iofile,tempx,tempy,/double

; initialize the x_dat and y_dat fields from the data just read
self.x_dat = ptr_new(tempx)
self.y_dat = ptr_new(tempy)

; initialize the model field from the data just read
self.model = ptr_new(dblarr(n_elements(tempy)))


; read the parameter setup from the inital input file
readcol,filename,inval, lowlim, uplim, pardescr, instepw, FORMAT='(D,D,D,A,D)',skipline=1

; initialize the number of parameters by getting the number of elements in inval
self.n_par = n_elements(inval)

; print number of parameters
print,'Number of parameters : ', self.n_par

; initialize the parameter array with the array of the initial values read from the file
self.params = ptr_new(inval)

; initialize the array containing the hitherto best parameters with the initial values
self.params_best = ptr_new(inval)

; initialize the lower and upper boundaries of the parameters 
self.params_minmax = ptr_new(dblarr(self.n_par, 2))
(*self.params_minmax)[*,0] = lowlim
(*self.params_minmax)[*,1] = uplim

; initialize the step widths from the fractions given in the setup file
self.params_step = ptr_new(dblarr(self.n_par))
*self.params_step = (uplim - lowlim)*instepw 

; initialize the parameter descriptions
self.params_descr = ptr_new(pardescr)
self.params_ar = ptr_new(lonarr(self.n_par,2))

;initialize the parameter record arrays
self.params_distr = ptr_new(dblarr(self.n_par))
self.params_distr_buf = ptr_new(dblarr(self.n_par))


end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_ar, ind=ind, all=all, global=global
; returns acceptance rates for specific or all parameters
; ind : if set, gives the index + 1 of the parameter for which to return the acceptance rate
; all : if set, acceptance rates for all parameters are returned
; global : if set, the global acceptance rate is returned 

;check if ind was set
if keyword_set(ind) then begin
   ; ind was set -> calculate the acceptance rate for the parameter with the index !!! ind - 1 !!!
   temp = ((*self.params_ar)[ind-1,0])/(1.0D*total((*self.params_ar)[ind-1,*]))
endif else begin
   ; ind was not set 
   if keyword_set(all) then begin
      ; all was set -> calculate the acceptance rates for all parameters
      temp = ((*self.params_ar)[*,0])/(1.0D*total(*self.params_ar,2))
   endif else begin
      ; all was not set
      if keyword_set(global) then begin
         ; global was set -> calculate the overall acceptance rate using the fields accept and reject
         temp = 1.0D*self.accept/(self.accept + self.reject)
	 ; print ...
	 print,'ACCEPT : ', self.accept
	 print,'REJECT : ', self.reject
      endif else begin
         ;oops ....
         print,'Warning in prob_obj::get_ar - arguments missing'
         temp = -1D
      endelse
   endelse
endelse

; return the calculated value
return, temp

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_best
; returns probability of best parameter values yet
   
   return, self.prob_best
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_hist, ind=ind, nb=nb
; returns a 2D-array containing the current parameter distribution as approximated by a histogram
; ind : gives the index + 1 for the parameter for which to return the current distribution
; nb : sets the number of bins used in the histogram

; get the requested parameter record for the parameter with the index !!! ind - 1 !!!
dat = (*self.params_distr)[ind-1,*]
   
; calculate the binwidth from the number of bins nb and the upper and lower parameter values yet encountered
binwidth = 1.0D*(max(dat) - min(dat))/nb

; calculate the histogram
h = histogram(dat, nbins=nb)

; calculate the indices corresponding to the histogram values
indices = min(dat) + dindgen(nb)*binwidth

; normalize the histogram
density = h/total(h)

; make a 2-dim array 
hist = [[indices],[density]]
   
; return the calculated data   
return, hist

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_data
; returns the ordinate values of observations (e.g., intensity, power)

   return, *self.y_dat
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_minmax, ind=ind, all=all
; returns an array of the lower and upper limits of specific or all parameters
; ind : if set, gives the index + 1 of the parameter for which to return the limits (returns a 1D-array : [lower, upper])
; all : if set, limits for all parameters are returned (returns a 2D-array : [[parameter indices],[lower, upper]] )

;check if ind was set
if keyword_set(ind) then begin

   ;ind was set -> return parameter limits for the parameter with the index !!! ind - 1 !!!
   temp = reform((*self.params_minmax)[ind-1,*])
   
endif else begin

   ;ind was not set
   if keyword_set(all) then begin
   
      ;all was set -> get all parameter limits
      temp = *self.params_minmax
      
   endif else begin
   
      ;all was not set -> oops
      print,'Warning in prob_obj::get_minmax - arguments missing'
      temp = -1D
   endelse
endelse

;return the data
return, temp


end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_model
; returns the ordinate values of the model

   return, *self.model
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_n_par
; returns the number of parameter values

   return, self.n_par
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_param, ind=ind, all=all, best=best
; returns the current values of specific or all parameters
; ind : if set, gives the index + 1 of the parameter for which to return the current value
; all : if set, all parameter values are returned
; best : if set, the best parameter values yet found are returned

;check if best was set
if keyword_set(best) then begin
  
   ; best was set -> return the best parameter values yet
   temp = *self.params_best
   
endif else begin
  
   if keyword_set(ind) then begin
   
      ; ind was set -> return the parameter value for the parmameter with index !!! ind - 1 !!!
      temp = (*self.params)[ind-1]
      
   endif else begin   
      if keyword_set(all) then begin
      
         ; all was set -> return all current parameter values
         temp = *self.params
	 
      endif else begin
      
         ; oops ...	 
         print,'Warning in prob_obj::get_param - arguments missing'
         temp = -1D
      endelse
   endelse
endelse

;return the data
return, temp

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_par_descr, ind=ind, all=all
; returns the descriptions of specific or all parameters in string format
; ind : if set, gives the index + 1 of the parameter for which to return the description
; all : if set, all parameter descriptions are returned

;check if ind was set
if keyword_set(ind) then begin

   ; ind was set -> return the name of the parameter with the index !!! ind - 1 !!!
   temp = (*self.params_descr)[ind-1]
   
endif else begin   
   if keyword_set(all) then begin
   
      ; all was set -> return all parameter descriptions
      temp = *self.params_descr
      
   endif else begin
   
      ; oops
      print,'Warning in prob_obj::get_par_descr - arguments missinge'      
      temp = 'oops'
   endelse
endelse
; return the data
return, temp

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_probability
;returns the probability of the most recently evaluated parameter values

   return, self.prob
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_seed
;returns the objects seed variable for random number generation

   return, *self.seed
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_steps, ind=ind, all=all
; returns the current values of specific or all .9ol
; ind : if set, gives the index + 1 of the parameter for which to return the parameter stepwidth
; all : if set, all parameter stepwidths are returned

;check if ind was set
if keyword_set(ind) then begin

   ; ind was set -> return stepwidths for the parameter with index !!! ind - 1 !!!
   temp = (*self.params_step)[ind-1]   
endif else begin
   if keyword_set(all) then begin
   
      ;all was set -> return all parameter stepwidths 
      temp = *self.params_step
      
   endif else begin
   
      ;oops
      print,'Warning in prob_obj::get_steps - arguments missing'
      temp = -1D
   endelse
endelse

; return data
return, temp

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

function prob_obj::get_x
; returns the abscissa values of the observations (e.g., date, wavelength)

   return, *self.x_dat
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_ar, new_ar, ind=ind, all=all, clear=clear
; sets the acceptance for specific parameters or the complete model 
; new_ar : the state of acceptance to be stored (0 : rejected, 1 : accepted)
; ind : if set, gives the index + 1 of the parameter for which to set the acceptance 
; all : if set, the acceptance for the whole model is set 
; clear : if set, sets all acceptance records to zero before setting the new value


;check if ind was set
if keyword_set(ind) then begin
   
   ;if clear was set, set the acceptance/rejection parameters to 0
   if keyword_set(clear) then (*self.params_ar)[*,*] = 0L   
   
   ;check if (new_ar = 0), which means reject, or (new_ar = 1), which means accept, 
   ;and increment the corresponding variables
   if new_ar eq 0 then (*self.params_ar)[ind-1,1]++ $
   else (*self.params_ar)[ind-1,0]++
   
endif else begin
   if keyword_set(all) then begin
   
      ; all was set -> change the "global" fields accept and reject
      if keyword_set(clear) then begin         
      
         ; clear was set -> set the global fields to 0
         self.accept = 0L
	 self.reject = 0L	  
      endif
      
      if new_ar ne -1 then begin
         ; change the global fields
         if new_ar eq 0 then self.reject++ else self.accept++   
      endif
   endif else begin
   
      ; all was not set, ind was not set -> change acceptance variables for every parameter
      if keyword_set(clear) then begin 
      
         ; clear was set -> set all acceptance/rejection parameters to 0
         (*self.params_ar)[*,*] = 0L
      endif else begin
      
         ; oops ...
         print,'Warning in prob_obj::set_ar - arguments missing'
      endelse
   endelse
endelse

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_best, new_best
; sets the best probability achieved by the model
; new_best : the value to be set

   self.prob_best = new_best
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_data, new_data, init=init
; sets the (observed) data used for fitting
; new_data : the data to be set 
; init : if set, (re-)initializes the heap variable 

; check if init was set
if keyword_set(init) then begin

   ; init was set -> deallocate memory and delete any previous contents y_dat is pointing at
   ptr_free, self.y_dat
   
   ; allocate memory for new_data and let y_dat point at new_data
   self.y_dat = ptr_new(new_data)
   
endif else begin
  
   ; init was not set -> just change the contents of *y_dat
   *self.y_dat = new_data
   
endelse

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_minmax, new_minmax, ind=ind, all=all
; sets the upper and lower limits of specific or all parameters in the model
; new_minmax : the limits to be set
; ind : if set, gives the specific index + 1 of the parameter for which to set the limits (requires new_minmax to be a 1D array)
; all : if set, limits are set for all parameters (requires new_minmax to be a 2D array)

;check if ind was set
if keyword_set(ind) then begin

   ; ind was set -> return the lower and upper limits as a 1D-array for the parameter with the index !!! ind - 1 !!!
   (*self.params_minmax)[ind-1,*] = new_minmax
   
endif else begin       
   if keyword_set(all) then begin
   
     ; all was set -> return all lower and upper limits as a 2D-array
     ( *self.params_minmax)[*,*] = new_minmax
     
   endif else begin
   
      ; oops ...
      print,'Warning in prob_obj::set_minmax - arguments missing'
      
   endelse
endelse

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_model, new_model, init=init
; sets the model used for fitting
; new_model : the model to be set 
; init : if set, (re-)initializes the heap variable 

;check if init was set
if keyword_set(init) then begin
   
   ;init was set -> deallocate memory and delete any previous contents model is pointing at
   ptr_free, self.model
   
   ; allocate memory for new_data and let y_dat point at new_model
   self.model = ptr_new(new_model)
   
endif else begin

   ; init was not set -> just change the contents of *new_model
   *self.model = new_model
   
endelse

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_n_par, new_n_par
; sets the number of parameters 
; new_n_par : new number of parameters 

   self.n_par = n_par
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_param, new_param, ind=ind, all=all, init=init, best=best
; sets the value of specific or all parameters
; new_param : the new value(s) for the parameter(s)
; ind : if set, gives the index + 1 of the parameter for which to set the value
; all : if set, all parameter values are set
; init : if set, (re-)initializes the heap variable (only works if 'all' is also set) 
; best : if set, only stores the submitted parameters as "best parameters yet"

; check if best was set
if keyword_set(best) then begin

   ;best was set -> set new_param as new content of *params_best
   *self.params_best = new_param
   
endif else begin  

   ;best was not set
    if keyword_set(ind) then begin
      
      ;ind was set -> change the parameter value for the parameter with the index !!! ind - 1 !!! 
      (*self.params)[ind-1] = new_param
      
   endif else begin
   
      ; ind was not set
      if keyword_set(all) then begin
      
         ;all was set -> change all parameter values
         if keyword_set(init) then begin
	   
	    ;init was set -> re-initialize params
            ptr_free, self.params
            self.params = ptr_new(new_param)
	    
         endif else begin
	 
	    ; init was not set -> change the content of *params to new_param 
            *self.params = new_param
	    
         endelse
      endif else begin
      
         ; oops ...
         print, 'Warning in prob_obj::set_param - arguments missing'
	 
      endelse
   endelse
endelse

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_par_descr, new_par_descr, ind=ind, all=all, init=init
; sets the description of specific or all parameters
; new_param : the new description(s) for the parameter(s)
; ind : if set, gives the index + 1 of the parameter for which to set the description
; all : if set, all parameter descriptions are set
; init : if set, (re-)initializes the heap variable (only works if 'all' is also set) 
  
  
;check if ind was set  
if keyword_set(ind) then begin

   ; ind not set -> change the description of the parameter with the index !!! ind - 1 !!!
   (*self.params_descr)[ind-1] = new_par_descr
   
endif else begin

   ; ind was not set 
   if keyword_set(all) then begin
      
      ;all was set -> change all parameter descriptions
      if keyword_set(init) then begin
       
         ;init was set -> re-initialize all parameter descriptions
         ptr_free, self.params_descr	 
         self.params_descr = ptr_new(new_par_descr)
	 
      endif else begin
      
         ; change the content of *params_descr to new_par_descr
         *self.params_descr = new_par_descr
	 
      endelse
   endif else begin
   
      ; oops ...
      print, 'Warning in prob_obj::set_par_descr - arguments missing'
      
   endelse
endelse

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_probability, new_prob
; sets the probability of the current model
; new_prob : the value to be set

   self.prob = new_prob
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_seed, new_seed
; sets the seed used for random number generation
; new_seed : the value of the seed to set

   *self.seed = new_seed
end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_steps, new_steps, ind=ind, all=all, init=init
; sets the step widths of specific or all parameters
; new_steps : the new step width(s) for the parameter(s)
; ind : if set, gives the index + 1 of the parameter for which to set the step width
; all : if set, all parameter step widths are set
; init : if set, (re-)initializes the heap variable (only works if 'all' is also set) 
  
  
; check if ind was set  
if keyword_set(ind) then begin

   ; ind was set -> change the stepwidth for the parameter with the index !!! ind - 1 !!!
   (*self.params_step)[ind-1] = new_steps
   
endif else begin

   ; ind was not set 
   if keyword_set(all) then begin
   
      ; all was set -> change all parameter descriptions 
      if keyword_set(init) then begin
      
         ;init was set -> re-initialize all parameter descriptions
         ptr_free, self.params_step
         self.params_step = ptr_new(new_steps)
	 
      endif else begin
       
         ; init was not set -> change the content of *params_step to new_steps
         *self.params_step = new_steps
	 
      endelse
   endif else begin
   
      ; oops ...
      print, 'Warning in prob_obj::set_steps - arguments missing'
      
   endelse
endelse

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj::set_x, new_x, init=init
; sets the (observed) abscissa data used for fitting
; new_x : the abscissa data to be set 
; init : if set, (re-)initializes the heap variable 

; check if init was set
if keyword_set(init) then begin

   ; init was set -> re-initialize x_dat
   ptr_free, self.x_dat
   self.x_dat = ptr_new(new_x)
   
endif else begin

   ; init was not set -> change the content of *x_dat to x_dat
   *self.x_dat = new_x
   
endelse

end

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

;-------------------------------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-------------------------------------------------------------------------------------------

pro prob_obj__DEFINE
;provides the definition of the standard object "prob_obj":
;
;n_par		    : number of parameters
;accept     	    : number of accepted steps for MCMC (after calibration)
;reject     	    : number of rejected steps for MCMC (after calibration)
;prob	    	    : probability of the most recently evaluated parameter values
;prob_best  	    : probability of best parameter values yets
;seed	    	    : pointer to 1D-array as seed for random number generators 
;params     	    : pointer to 1D-array containing the parameter values
;params_best 	    : pointer to 1D-array containing the best parameters yet
;params_distr       : pointer to 2D-array containing the resulting values for each parameter
;params_distr_buf   ; pointer to 2D-array which acts as a buffer for storing parameter values for quickly
;params_descr	    : pointer to 1D-array containing description of parameters in string format
;params_priors	    : pointer to 1D-array containing the current value of the prior probability for each parameter
;params_ar  	    : pointer to 2D-array containing the number of accepted/rejected steps for individual parameters
;params_step 	    : pointer to 1D-array containing the current step width for individual parameters
;params_minmax 	    : pointer to 2D-array containing the lower and upper limits for each parameter
;x_dat      	    : pointer to 1D-array containing the abscissa of the observations (e.g., date, frequency, wavelength)
;y_dat	    	    : pointer to 1D-array containing the ordinate of the observations (e.g., intensity, power, radial vel.)
;model	    	    : pointer to 1D-array containing the model values corresponding to the observed abscissa values

;define the class 
struct = { prob_obj, n_par:0L, accept:0L, reject:0L, prob:0D, prob_best:0D, seed:ptr_new(/allocate_heap), params:ptr_new(/allocate_heap), $
           params_best  :ptr_new(/allocate_heap),params_distr:ptr_new(/allocate_heap), params_distr_buf:ptr_new(/allocate_heap), $
	   params_priors:ptr_new(/allocate_heap), params_descr:ptr_new(/allocate_heap), params_ar:ptr_new(/allocate_heap), params_step:ptr_new(/allocate_heap), $
	   params_minmax:ptr_new(/allocate_heap), x_dat:ptr_new(/allocate_heap), y_dat:ptr_new(/allocate_heap), model:ptr_new(/allocate_heap)} 

end


