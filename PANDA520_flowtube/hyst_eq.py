'''solution of deliquescence and efflorescence RH, generated by eqn_pars.py in fully functioning mode, or by ui_check.py in testing mode'''
# module to estimate deliquescence and efflorescence relative humidity as a function of temperature
# File Created at 2022-01-12 16:40:44.230081

# function for deliquescence
def drh(TEMP):
	
	# inputs: -----------------
	# TEMP - temperature in chamber (K)
	# ---------------------------
	
	# deliquescence relative humidity (fraction 0-1)
	DRH = 0.*TEMP
	return(DRH)

# function for efflorescence
def erh(TEMP):
	
	# inputs: -----------------
	# TEMP - temperature in chamber (K)
	# ---------------------------
	
	# efflorescence relative humidity (fraction 0-1)
	ERH = 0.*TEMP
	return(ERH)
