# @File: IOframe.py
# @Author: Saeed Ahmad (email:saeediitb@gmail.com)
# @Date:   2023-03-20 16:57:48
# @Last Modified time: 2023-07-03 12:50:39
import numpy as np
import pandas as pd 
from ast import literal_eval
import csv

def frame_out(n,time,sys,name_file='shell.csv'):
	shell = sys
	#genome=sys.sys.particles[0]
	with open(name_file,'a') as f:
		f.write("#nthItr:%d\n"%n)
		f.write("#FPT:%f\n"%time)
		#f.write("#genome:{},{},{}\n".format(genome.position[0],genome.position[1],genome.position[2]))
		f.write("#vertex:%d\n"%len(shell.vertex))
		shell.vertex.to_csv(f)
		f.write("#line:%d\n"%len(shell.line))
		shell.line.to_csv(f)
		f.write("#triangle:%d\n"%len(shell.triangle))
		shell.triangle.to_csv(f)
def frame_in(param,shell,genome):
	with open('../{}.csv'.format(param['infile']),'r') as f:
		frame_row=0
		for row_num, line in enumerate(f, 1):
			if 'frame' in line:
				frame_row = row_num
				print('frame_row:',frame_row)
			if row_num>frame_row and line[0]=='#':
				data_row=row_num
				data_name,data_nrow=line.split(':')
				if data_name[1:]=='genome':
					genome.particles=np.array([[float(x) for x in data_nrow.split(',')]])
				else:
					print(data_name[1:],int(data_nrow))
					df=pd.read_csv('../{}.csv'.format(param['infile']), skiprows = data_row, index_col=0,header=0, nrows=int(data_nrow))
					if data_name[1:]=='vertex':
						for vi,row in df.iterrows():
							[x,y,z,l,t,edge]=row
							if l=='inf':
								shell.vertex.loc[vi]=[x,y,z,np.inf,np.inf,np.inf]
							else:
								shell.vertex.loc[vi]=[x,y,z,literal_eval(l),literal_eval(t),edge]
					elif data_name[1:]=='line':shell.line=df
					elif data_name[1:]=='triangle':shell.triangle=df
					else: raise NameError('Unknown input data')
def energy_frame(frame,shell_Energy,rates = [0],rate=0, back_rate=[0], mc_type = 'None',  col_name=True):
	with open ('enegy_w_time.csv', 'a') as f:
		shell_Energy=shell_Energy.rename(index={0:frame})
		# Define data for multiple columns
		new_columns = {
    		'rates': [rates],
			'rate': [rate],
			'back_rate':[back_rate],
    		'type': [mc_type]
		}
		shell_Energy = shell_Energy.assign(**new_columns)
		shell_Energy.to_csv(f,header=col_name)
def First_passage_Time(n, time, name = False):
	data=[n,time]
	with open ('FPT.csv', 'a') as f:
		csv_writer = csv.writer(f)
		if name:
			csv_writer.writerow(('attempt','time'))
		else:
			csv_writer.writerow((data))

