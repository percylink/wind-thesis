'''
PAL 7/15/13 mod 8/1/2013
modify geo_em file to change CA forests to savanna
'''

import numpy as np
import os


class ModifierGeo:

	def __init__(self, file_in=None, varsmod=[], mod_vals={}, mod_conditions={}):
		
		# open input file
		if file_in is None:
			raise Exception('You must specify an input file')
		self.filename_in = file_in

		dum = file_in.split('.')
		self.filename_out = '.'.join(dum[:-1])+'_MOD.'+dum[-1]
		self.file_out = open(self.filename_out, 'w')
		
		# save variables to modify, and values to change to, and conditions for changing
		self.varsmod = varsmod
		print('variables to modify: '+str(self.varsmod))
		self.mod_vals = mod_vals
		self.mod_conditions = mod_conditions

		if 'LU_INDEX' not in self.mod_vals:
			self.land_use_cat = 7  # default land use category if none specified
			print('land use category: 7 (default)')
		else:
			print('land use category: '+str(self.mod_vals['LU_INDEX']))
			self.land_use_cat = self.mod_vals['LU_INDEX']
		if 'SCT_DOM' not in self.mod_vals:
			self.soil_cat = 8  # default soil category if none specified
			print('soil category: 8 (default)')
		else:
			print('soil category: '+str(self.mod_vals['SCT_DOM']))
			self.soil_cat = self.mod_vals['SCT_DOM']
		self.make_default_mod_vals()
		
		# open output file
		with open(self.filename_in, 'r') as f:
			self.alltxt = f.read()
		ixhead = self.alltxt.find("data:")
		self.header = self.alltxt[:ixhead]
		self.body = self.alltxt[ixhead:]
		self.get_dimensions()

	def make_default_mod_vals(self):
		self.default_mod_vals = {'HGT_M': 0, 'LANDMASK': 1, 'SOILTEMP': 289., \
								 'ALBEDO12M': 15., 'GREENFRAC': 0.8, 'LAI12M': 2., \
								 'SLOPECAT': 2, 'CON': 0.1, 'VAR': 10., 'OA1': 0., \
								 'OA2': 0., 'OA3':0., 'OA4': 0., 'VAR_SSO': 0., \
								 'LU_INDEX': self.land_use_cat, \
								 'SCT_DOM': self.soil_cat, 'SCB_DOM': self.soil_cat}

	def get_dimensions(self):

		h = self.header.split('\n')
		self.south_north = None
		self.west_east = None
		for hh in h:
			htmp = hh.strip().split(' = ')
			if htmp[0] == 'south_north':
				self.south_north = int(htmp[1].replace(';', '').strip())
			elif htmp[0] == 'west_east':
				self.west_east = int(htmp[1].replace(';', '').strip())
		if self.south_north is None:
			raise Exception('Could not find south_north in the header')
		if self.west_east is None:
			raise Exception('Could not find west_east in the header')

	def write_header(self):

		self.file_out.write(self.header)
		self.file_out.write('data:\n\n Times =\n  "0000-00-00_00:00:00" ;\n\n')

	def get_varnames(self):

		# get list of variable names
		headtmp = self.header
		variables = []
		while headtmp.find("float ")>0:
			ixtmp = headtmp.find("float ")
			headtmp = headtmp[ixtmp+6:]
			ixtmp2 = headtmp.find("(")
			variables.append(headtmp[:ixtmp2].strip())
		self.vars = variables

	def get_var_vals(self, variable, body):

		ixvarS = body.find(variable)
		ixvarE = body.find(';')+1
		#print vv, ixvarS, ixvarE
		vartmp = body[ixvarS:ixvarE-1]
		varlines = vartmp.splitlines()
		#print len(varlines)
		varlist = []
		for jj,var in enumerate(varlines[1:]):  # skip first line, which is just "<var> ="
			var = var.strip()
			if var[-1] == ',':  #jj<len(varlines)-2:
				varlist = varlist+var[:var.rfind(',')].split(',')
			else:
				varlist = varlist+var.split(',')
		currvar = np.array(varlist).astype(float)

		vartxt = body[ixvarS:ixvarE]
		body = body[ixvarE:]

		return body, vartxt, currvar

	def modify_variables(self, variable, currvar):

		if variable in ['HGT_M', 'LANDMASK', 'SOILTEMP', 'ALBEDO12M', 'GREENFRAC', \
						'LAI12M', 'SLOPECAT', 'CON', 'VAR', 'OA1', 'OA2', 'OA3', \
						'OA4', 'VAR_SSO', 'LU_INDEX', 'SCT_DOM', 'SCB_DOM']:

			if variable in self.mod_vals:
				mod_val = self.mod_vals[variable]
			else:
				mod_val = self.default_mod_vals[variable]

			if variable in self.mod_conditions:
				mask = exec(self.mod_conditions[variable])
				currvar[mask] = mod_val
				print("modifying "+variable+" to "+str(mod_val)+" where "+self.mod_conditions[variable])
			else:
				currvar = currvar * 0. + mod_val
				print("modifying "+variable+" to "+str(mod_val))

			if variable in ['LANDMASK', 'SLOPECAT', 'LU_INDEX', 'SCT_DOM', 'SCB_DOM']:
				currvar = currvar.astype(int)

		elif variable == 'LANDUSEF':
			print("modifying LANDUSEF to be 1 in category "+str(self.land_use_cat)+" and 0 for other categories")
			for i in xrange(24):
				ixstart = self.west_east*self.south_north*i
				ixend = self.west_east*self.south_north*(i+1)
				if i == self.land_use_cat-1:
					currvar[ixstart:ixend] = currvar[ixstart:ixend]*0. + 1.
				else:
					currvar[ixstart:ixend] = currvar[ixstart:ixend]*0.

		elif (variable == 'SOILCTOP') or (variable == 'SOILCBOT'):
			print("modifying "+variable+" to be 1 in category "+str(self.soil_cat)+" and 0 for other categories")
			for i in xrange(16):
				ixstart = self.west_east*self.south_north*i
				ixend = self.west_east*self.south_north*(i+1)
				if i == self.soil_cat-1:
					currvar[ixstart:ixend] = currvar[ixstart:ixend]*0. + 1.
				else:
					currvar[ixstart:ixend] = currvar[ixstart:ixend]*0.

		else:
			raise Exception('Variable '+variable+' has no code for modifying')

		return self.make_vartxt(variable, currvar)
		
	def make_vartxt(self, variable, var_array):
		# using user-modified array, print lines with 6 numbers, separated by commas,
		# until the end of the array, then print the remainder, and end with semicolon.
		vartxt = ' '+variable+' =\n'
		k = 6
		while len(var_array) > kk:
			linetmp = var_array[:kk]
			# convert to comma-separated string and write to vartxt
			strtmp = ''
			for ss in linetmp: 
				strtmp = strtmp+' '+str(ss)+','
			vartxt += strtmp+'\n'
			var_array = var_array[kk:]
		# output last line
		strtmp = ''
		for ix,ss in enumerate(var_array):
			if ix<len(var_array)-1:
				strtmp = strtmp+' '+str(ss)+','
			else:
				strtmp = strtmp+' '+str(ss)+';'
		vartxt += strtmp+'\n\n'
		return vartxt

	def get_variables(self):
		body = self.body
		body = body[body.find(';')+1:]
		self.orig_var_text = {}
		self.var_arrays = {}
		for ii,vv in enumerate(self.vars):
			# find variable
			body, vartxt, currvar = self.get_var_vals(vv, body)
			self.orig_var_text[variable] = vartxt
			self.var_arrays[variable] = currvar

	def save_key_variables(self):
		self.landmask = self.var_arrays['LANDMASK']
		self.xlat = self.var_arrays['XLAT_M']
		self.xlong = self.var_arrays['XLONG_M']

	def write_variables(self):
		# for each variable, either write original lines
		# or write new values, defined by user
		for ii,vv in enumerate(self.vars):
			# print variable to new file
			print vv,
			if vv in self.varsmod:
				print(' - modifying')
				vartxt = self.modify_variables(vv, self.var_arrays[vv])
			else:
				print(' - leaving as is')
				vartxt = self.orig_var_text[variable]
			self.file_out.write(' '+vartxt+'\n\n')
	
	def finish_file(self):
		self.file_out.write("}")
		self.file_out.close()

	def run(self):

		self.write_header()
		self.get_varnames()
		self.get_variables()
		self.save_key_variables()
		self.write_variables()
		self.finish_file()


if __name__ == "__main__":

	land_use_cat = 7  # grassland
	soil_cat = 8  # silty clay loam

	file_list = ['geo_em.d01', 'geo_em.d02', 'geo_em.d03']
	varsmod = ['HGT_M', 'LANDMASK', 'LANDUSEF', 'LU_INDEX', 'SOILTEMP', 'SOILCTOP', \
			   'SCT_DOM', 'SOILCBOT', 'SCB_DOM', 'ALBEDO12M', 'GREENFRAC', 'LAI12M', \
			   'SLOPECAT', 'CON', 'VAR', 'OA1', 'OA2', 'OA3', 'OA4', 'VAR_SSO']

	for f in file_list:
		
		print f
		os.system("ncdump "+f+".nc > "+f+".cdl")

		m = ModifierGeo(file_in=f+'.cdl', varsmod=varsmod, land_use_cat=land_use_cat, soil_cat=soil_cat)
		m.run()

		print "ncgen -o "+f+"_MOD.nc "+f+"_MOD.cdl"

		os.system("rm "+f+"_MOD.nc")
		os.system("ncgen -o "+f+"_MOD.nc "+f+"_MOD.cdl")