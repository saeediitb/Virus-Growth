import hoomd
import virus.capsid as cp
import numpy as np
from termcolor import colored, cprint
import copy 
import pickle
import math
from virus.IOframe import *
ONEDGE = 32767

def kmc_moves(self):
	shell = pickle.loads(pickle.dumps(self.shell))
	shellE = pickle.loads(pickle.dumps(self.shellE))
	#-------------------Diffusion------rates----------------#
	diff_multiplier =1
	growvlist= self.mc_info(self.shell)
	#-------------------Vertex---------close------rates-----#
	mc_vrates = []
	if len(growvlist) > 0:
		kmc_vmove_delE, kmc_vmove_Index = self.kmc_vmove(shell, shellE, event=None)
		mc_vrates = map_List(kmc_vmove_delE, multiplier = diff_multiplier)
		self.av_kmc_vmove_delE = kmc_vmove_delE
		self.av_mc_vrates = mc_vrates
		self.av_mc_vrates_invr = Inverse_rate(mc_vrates,kmc_vmove_delE)

	#-------------------Trimer-----------Attach--------------#
	grow_multiplier = 1
	mc_grow_rates = []
	growlist = shell.get_all_grows(grow=1)
	if len(growlist) > 0:
		kmc_grow_delE, kmc_grow_Index = self.kmc_grow(shell, shellE, event=None)
		mc_grow_rates = map_List(kmc_grow_delE, multiplier = grow_multiplier)
		self.av_mc_grow_rates = mc_grow_rates
		self.av_mc_grow_rates_invr = Inverse_rate(mc_grow_rates,kmc_grow_delE)
	#-------------------Trimer-----------Detach--------------#
	detach_multiplier = 1
	mc_detach_rates = []
	detach_tlist = shell.detach_tlist()
	if len(detach_tlist) and len(shell.triangle) > 1:
		kmc_detach_delE, kmc_detach_Index = self.kmc_detach(shell,shellE,event=None)
		mc_detach_rates = map_List(kmc_detach_delE, multiplier = detach_multiplier)
		self.av_mc_detach_rates = mc_detach_rates
		self.av_mc_detach_rates_invr = Inverse_rate(mc_detach_rates,kmc_detach_delE)
	#-------All possible rates-------------------------------#
	rates = mc_vrates + mc_grow_rates + mc_detach_rates
	partial_sum = [(sum(rates[:j+1]))/(sum(rates)) for j in range(len(rates))] # making a list for selecting jth event
	
	num = np.random.rand()
	jth_event = track_jth_event(num, partial_sum)# selecting jth event
	self.event_h_rate = rates[jth_event]
	lists_dict, list_name, position = find_which_event(mc_vrates, mc_grow_rates, mc_detach_rates,target_number=rates[jth_event])
	dt = -(math.log(1-np.random.rand()))/(np.sum(rates))
	self.time += dt
	if list_name == 'mc_vrates':
		event_index = kmc_vmove_Index[position]
		self.kmc_vmove(shell, shellE, event = event_index)
	elif list_name == 'mc_grow_rates':
		event_index = kmc_grow_Index[position]
		self.kmc_grow(shell, shellE, event = event_index)
	else:
		event_index = kmc_detach_Index[position]
		self.kmc_detach(shell, shellE, event = event_index)

def kmc_vmove(self, shell, shellE, event):
	growvlist = self.mc_info(shell)
	growlist = growvlist
	if len(growlist) > 0:
		if event == None:
			kmc_vmove_delE = []
			kmc_vmove_Index = []
			for i in growlist.index:
				growi = i
				self.shell = shell.shell_grow(growlist, growi)
				max_run = self.try_relax()
				if max_run < self.param['maxrun']:
					delE = self.shellE.loc[0,'Eshell']-shellE.loc[0,'Eshell']
					kmc_vmove_delE.append(delE)
					kmc_vmove_Index.append(growi)
			return kmc_vmove_delE, kmc_vmove_Index
		else:
			for i in growlist.index:
				if i == event:
					growi = event
					self.shell = shell.shell_grow(growlist, growi)
					max_run = self.try_relax()
					#cprint("Vertex move accepted",'green')
					#frame_out(self.time, self.shell)
					energy_frame(self.time, self.shellE,rates=self.av_mc_vrates,rate=self.event_h_rate,back_rate=self.av_mc_vrates_invr,mc_type='vc',col_name=False)
					return 1		
def kmc_grow(self, shell, shellE, event):
	growlist = shell.get_all_grows(grow=1)
	if len(growlist) > 0:
		if event == None:
			kmc_grow_delE = []
			kmc_grow_Index = []
			for i in growlist.index:
				growi = i
				self.shell = shell.shell_grow(growlist, growi)
				max_run = self.try_relax()
				if max_run < self.param['maxrun']:
					delE = self.shellE.loc[0,'Eshell']-shellE.loc[0,'Eshell']-self.param['mu']
					kmc_grow_delE.append(delE)
					kmc_grow_Index.append(growi)
			return kmc_grow_delE, kmc_grow_Index
		else:
			for i in growlist.index:
				if i == event:
					growi = event
					self.shell = shell.shell_grow(growlist, growi)
					max_run = self.try_relax()
					#cprint("Grow accepted",'green')
					#frame_out(self.time, self.shell)
					energy_frame(self.time,self.shellE,rates=self.av_mc_grow_rates,rate=self.event_h_rate,back_rate=self.av_mc_grow_rates_invr,mc_type='g',col_name=False)
					return 3
def kmc_detach(self, shell, shellE, event):
	detach_tlist = shell.detach_tlist()
	if len(detach_tlist) and len(shell.triangle) > 1:
		if event == None:
			kmc_detach_delE = []
			kmc_detach_Index = []
			for i in detach_tlist:
				ti = i
				self.shell = shell.shell_detach(ti)
				max_run = self.try_relax()
				if max_run < self.param['maxrun']:
					delE = self.shellE.loc[0,'Eshell']-shellE.loc[0,'Eshell']+self.param['mu']
					kmc_detach_delE.append(delE)
					kmc_detach_Index.append(ti)
			return kmc_detach_delE, kmc_detach_Index
		else:
			for i in detach_tlist:
				if i == event:
					ti = event
					self.shell = shell.shell_detach(ti)
					max_run = self.try_relax()
					#cprint("Detach accepted",'green')
					#frame_out(self.time, self.shell)
					energy_frame(self.time,self.shellE,rates=self.av_mc_detach_rates,rate=self.event_h_rate,back_rate=self.av_mc_detach_rates_invr,mc_type='d',col_name=False)
					return -1
def map_List(kmc_list, multiplier = 1):
	rates = list(map(lambda x: multiplier*min(1,np.exp(-x)), kmc_list))
	return rates
def Inverse_rate(kmc_list,kmc_delE):
	multiplier = list(map(lambda x: np.exp(-x), kmc_delE))
	rates = [a * b for a, b in zip(kmc_list, multiplier)]
	return rates
def track_jth_event(num, partial_sum):
	num_array = np.array(partial_sum)
	sorted_indices = np.argsort(num_array)
	index = np.searchsorted(num_array[sorted_indices],num)
	return sorted_indices[index]
def	find_which_event(mc_vrates, mc_grow_rates, mc_detach_rates,target_number=None):
	lists_dict = {
		'mc_vrates': mc_vrates,
		'mc_grow_rates' : mc_grow_rates,
		'mc_detach_rates': mc_detach_rates,
	}
	for list_name, lst in lists_dict.items():
		if target_number in lst:
			position = lst.index(target_number)
			return lists_dict, list_name, position
