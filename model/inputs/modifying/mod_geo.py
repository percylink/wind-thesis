'''
PAL 7/15/13 mod 8/1/2013
modify geo_em file to change CA forests to savanna
'''

import numpy as np
import os


class ModifierGeo:

	def __init__(self, file_in=None, varsmod=[], land_use_cat=7, soil_cat=8):
		
		if file_in is None:
			raise Exception('You must specify an input file')
		self.filename_in = file_in

		dum = file_in.split('.')
		self.filename_out = '.'.join(dum[:-1])+'_MOD.'+dum[-1]
		self.file_out = open(self.filename_out, 'w')
		
		self.varsmod = varsmod
		print('variables to modify: '+str(self.varsmod))
		self.land_use_cat = land_use_cat
		print('land use category: '+str(land_use_cat))
		self.soil_cat = soil_cat
		print('soil category: '+str(soil_cat))
		
		with open(self.filename_in, 'r') as f:
			self.alltxt = f.read()
		ixhead = self.alltxt.find("data:")
		self.header = self.alltxt[:ixhead]
		self.body = self.alltxt[ixhead:]
		self.get_dimensions()

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

		if variable == 'HGT_M':
			currvar = currvar * 0.
		elif variable == 'LANDMASK':
			currvar = currvar * 1.
			currvar = currvar.astype(int)
		elif variable == 'LU_INDEX':
			currvar = currvar*0 + self.land_use_cat
			currvar = currvar.astype(int)
		elif variable == 'LANDUSEF':
			for i in xrange(24):
				ixstart = self.west_east*self.south_north*i
				ixend = self.west_east*self.south_north*(i+1)
				if i == self.land_use_cat-1:
					currvar[ixstart:ixend] = currvar[ixstart:ixend]*0. + 1.
				else:
					currvar[ixstart:ixend] = currvar[ixstart:ixend]*0.
		elif variable == 'SOILTEMP':
			currvar = currvar*0. + 289.
		elif (variable == 'SOILCTOP') or (variable == 'SOILCBOT'):
			for i in xrange(16):
				ixstart = self.west_east*self.south_north*i
				ixend = self.west_east*self.south_north*(i+1)
				if i == self.soil_cat-1:
					currvar[ixstart:ixend] = currvar[ixstart:ixend]*0. + 1.
				else:
					currvar[ixstart:ixend] = currvar[ixstart:ixend]*0.
		elif (variable == 'SCT_DOM') or (variable == 'SCB_DOM'):
			currvar = currvar*0 + self.soil_cat
			currvar = currvar.astype(int)
		elif variable == 'ALBEDO12M':
			currvar = currvar*0. + 15.
		elif variable == 'GREENFRAC':
			currvar = currvar*0. + 0.8
		elif variable == 'LAI12M':
			currvar = currvar*0. + 2.
		elif variable == 'SLOPECAT':
			currvar = currvar*0. + 2
			currvar = currvar.astype(int)
		elif variable == 'CON':
			currvar = currvar*0. + 0.1
		elif variable == 'VAR':
			currvar = currvar*0. + 10.
		elif (variable == 'OA1') or (variable == 'OA2') or (variable == 'OA3') or (variable == 'OA4'):
			currvar = currvar*0.
		elif variable == 'VAR_SSO':
			currvar = currvar*0.
		else:
			raise Exception('Variable '+variable+' has no code for modifying')
			
		self.file_out.write(' '+variable+' =\n')
		# using user-modified array, print lines with 6 numbers, separated by commas,
		# until the end of the array, then print the remainder, and end with semicolon.
		kk = 6
		while len(currvar) > kk:
			linetmp = currvar[:kk]
			# convert to comma-separated string and write to file
			strtmp = ''
			for ss in linetmp: strtmp = strtmp+' '+str(ss)+','
			self.file_out.write(strtmp+'\n')
			currvar = currvar[kk:]
		# output last line
		strtmp = ''
		for ix,ss in enumerate(currvar):
			if ix<len(currvar)-1:
				strtmp = strtmp+' '+str(ss)+','
			else:
				strtmp = strtmp+' '+str(ss)+';'
		self.file_out.write(strtmp+'\n\n')

	def get_data(self):
		# for each variable, either write original lines
		# or write new values, defined by user
		body = self.body
		body = body[body.find(';')+1:]
		for ii,vv in enumerate(self.vars):

			print vv,
			
			# find variable
			body, vartxt, currvar = self.get_var_vals(vv, body)

			# print variable to new file
			if vv in self.varsmod:
				print(' - modifying')
				self.modify_variables(vv, currvar)
			else:
				print(' - leaving as is')
				self.file_out.write(' '+vartxt+'\n\n')
	
	def finish_file(self):
		self.file_out.write("}")
		self.file_out.close()

	def run(self):

		self.write_header()
		self.get_varnames()
		self.get_data()
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