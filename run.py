# @File: run.py
# @Author: Saeed AHmad (email:saeediitb@gmail.com)
# @Date:   2023-06-20 18:08:51
# @Last Modified (Saeed Ahmad) time: 2023-08-13 13:49:54
# This code is for versions python=3.10.8, panda=1.5.1, numpy=1.23.5, hoomd=3.5.0 
import hoomd
from hoomd import md
import numpy as np
import pandas as pd
import json
import glob
import virus
from virus.IOframe import *
from virus.system import System
from virus.system import Update_Snap
import sys,os,glob
from termcolor import colored,cprint
import pickle
import os
ONEDGE=32767

def test_merge(sys):
	#cprint("test merging...",'yellow')
	mergestatus = sys.shell.shell_merge()
	while mergestatus == 1:
		cprint('Force Merge Happen','yellow')
		sys.try_relax()
		mergestatus = sys.shell.shell_merge()
	if mergestatus == -1:
		#cprint("Assemble done from merge")
		sys.try_relax()
		#sys.setup_operation("final1.gsd")
	#else: cprint("no more merge!",'yellow')
	return mergestatus

def run():
	virus_sys = System(param)
	virus_sys.initial()
	snap = hoomd.Snapshot()
	snap.configuration.box = [1.5*param['L'], 1.5*param['L'], 1.5*param['L'], 0, 0, 0]
	snap.particles.types = ['A','B']
	snap.bonds.types = ['out_ver_to_out_ver']
	snap.dihedrals.types = ['capsomer']
	snap.angles.types = ['ds']
	virus_sys.sim.create_state_from_snapshot(snap)
	virus_sys.setup_potential() #
	virus_sys.setup_operation() #-----------------------------------
	virus_sys.setup_integrator() #-------------------------
	max_run = virus_sys.try_relax() #---------------------------
	#frame_out(virus_sys.time, virus_sys.shell)
	n=0
	shell=pickle.loads(pickle.dumps(virus_sys.shell))
	shellE = pickle.loads(pickle.dumps(virus_sys.shellE))
	for n in range(100000):
		virus_sys.time = 0.0
		virus_sys.shell = shell
		virus_sys.shellE = shellE
		energy_frame(virus_sys.time,virus_sys.shellE)#,rate=0.0,mc_type='type')
		while (not virus_sys.shell.test_shell_close()):
			virus_sys.kmc_moves()
			mergestatus = test_merge(virus_sys)
			if len(virus_sys.shell.triangle)>=25:
				break
		frame_out(n,virus_sys.time, virus_sys.shell)
		#frame_out(virus_sys.time, virus_sys.shell,name_file='final_shell.csv')
		if len(virus_sys.shell.triangle)==20:
			print(n, virus_sys.time)
			energy_frame(virus_sys.time,virus_sys.shellE,col_name=False)
			First_passage_Time(n, virus_sys.time)



	cprint(" Final Assemble Done!\nFinal Shell:",'white', attrs = ['reverse'])
	print('Congratulations: Run Completed')
class Logger(object):
	def __init__(self):
		self.terminal = sys.stdout
		self.log = open(param['outfile'], "a")
	def write (self, massage):
		self.terminal.write(massage)
		self.log.write(massage)
	def flush(self):
		# this flush method is needed for python3 compatibility.
		# this handles the flush command by doing nothing.
		# you might want to
		pass
param = sys.argv[1]
#cprint(param,'red')
with open(param) as f: param = json.load(f)
if param['outfile'] != 'NA':
	sys.stdout = Logger()
run()
