# @File: capsid.py
# @Author: Siyu Li (email:lisiyu211@gmail.com)
# @Date:   2018-03-20 18:05:02
# @Last Modified Saeed time: 2023-04-03 14:57:45 
import virus
import numpy as np
import pandas as pd
import json
import copy
from termcolor import colored,cprint
import pickle 
import csv

ONEDGE=32767

class shell():
	"""docstring for Shell"""
	def __init__(self,vertex,line,triangle,param):
		self.R0=param['R0']
		self.R1=param['R1']
		self.R2=param['R2']
		self.R3=param['R3']
		self.R=param['R']
		self.l0=param['l0']
		self.theta0=param['theta0']
		self.thetagrow=param['thetagrow']
		self.theta_lower=param['theta_lower']
		self.gamma=param['gamma']
		self.vertex=vertex
		self.line=line
		self.triangle=triangle
		self.n=len(self.vertex)
		self.particles=self.vertex[list('xyz')].values
		self.typeid=0
		self.bodyid=-1
		self.bondid=1
		self.dihedralid=0
		self.bonds=np.int_(self.line[['v0','v1']])
		self.dihedrals=[]
		self.bondE=param['bondE']
		self.pentE=param['pentE']
		self.l_tension=param['l_tension']
		self.hp_eps=param['hp_eps']
		self.thres_Ehp=-param['thres_T']*param['hp_eps']
		self.frame_start = 0
		self.frame_size = 0
		self.r1 = param['r1']
		self.r2 = param['r2']
	#****************For complete Capsid T1 or T3 ************#
	def capsid_shape(self):
		self.vertex=pd.DataFrame(columns=['x','y','z','l','t','edge'])
		self.line=pd.DataFrame(columns=['v0','v1','t','edge']).astype('int32')
		self.triangle=pd.DataFrame(columns=['v0','v1','v2','l0','l1','l2','edge']).astype('int32')
		self.vertex_value,self.line_value,self.triangle_value, self.frame_start,self.frame_size = get_footer('Initial_shell_T_30.csv',frame_end=self.frame_start)
		# For vertex
		for vi, element in enumerate(self.vertex_value, start = 1):
			self.vertex.loc[element[0]] = [element[1], element[2], element[3], element[4], element[5], element[6]]
		# For line
		for vi, element in enumerate(self.line_value, start = 1):
			self.line.loc[element[0]] = [element[1], element[2], element[3], element[4]]
		# For triangle
		for vi, element in enumerate(self.triangle_value, start = 1):
			self.triangle.loc[element[0]] = [element[1], element[2], element[3], element[4], element[5], element[6], element[7]]
	def initial(self):
		initialcor=[[-0.5*self.l0,-(np.sqrt(3)/6.)*self.l0, 0.,[],[]],\
					[0.,(np.sqrt(3)/3.)*self.l0, 0.,[],[]],\
					[0.5*self.l0, -(np.sqrt(3)/6.)*self.l0, 0.,[],[]]]
		for vi in [1,2,3]: self.add_vertex(vi,initialcor[vi-1],ONEDGE)
		for li,[v0,v1] in enumerate([[1,2],[2,3],[3,1]],start=1):self.add_line(li,v0,v1,1,ONEDGE)
		self.add_triangle(1,1,2,3,1,2,3,ONEDGE)
	def create_new_cor(self,l,theta=None):
		#create a new point for shell growth
		#suppose the new triangle attach to line l
		if theta is None:
			theta=self.thetagrow
		[v0,v1,t,edge]=self.line.loc[l]
		v2=self.other_vertex(t,v0,v1)
		[cor0,cor1,cor2]=np.array(self.vertex.loc[[v0,v1,v2],list('xyz')])
		a=cor1-cor0
		c=(cor0+cor1)*0.5
		n=self.tri_norm_vector(t)
		b=np.cross(a,n)
		#print(self.l0**2-(0.5*np.linalg.norm(a))**2)
		#print('value of tria_h')
		tria_h=0.2*np.sqrt(abs(self.l0**2-(0.5*np.linalg.norm(a))**2))
		b_hat=b/np.linalg.norm(b)
		n_hat=n/np.linalg.norm(n)
		new_p=c+(b_hat*tria_h*np.cos(theta))-(n_hat*tria_h*np.sin(theta));
		# print(colored('new_p:','yellow'),new_p)
		return new_p
	# *********************Add functions***************************
	def add_vertex(self,vi,cor,edge):
		#add new vertex coordinate to vertex table
		self.vertex.loc[vi]=[cor[0],cor[1],cor[2],[],[],edge]
	def add_line(self,li,v0,v1,t,edge):
		#add new line to line table
		self.line.loc[li]=[v0,v1,t,edge]
		for vi in [v0,v1]:self.vertex.loc[vi,'l'].append(li)
	def add_triangle(self,ti,v0,v1,v2,l0,l1,l2,edge):
		#add new triangle to triangle table
		self.triangle.loc[ti]=[v0,v1,v2,l0,l1,l2,edge]
		for vi in [v0,v1,v2]:self.vertex.loc[vi,'t'].append(ti)	
	# *********************Auxilary functions***************************
	def other_vertex(self,t,v0,v1):
		# function to find the third vertex in triangle t
		[tv0,tv1,tv2]=self.triangle.loc[t,['v0','v1','v2']]
		return tv0+tv1+tv2-v0-v1
	def vertex_near_vlt(self,v):
		left=self.line[(self.line.v1==v) & (self.line.edge==ONEDGE)].iloc[0]
		right=self.line[(self.line.v0==v) & (self.line.edge==ONEDGE)].iloc[0]
		return left.v0,left.name,left.t,right.v1,right.name,right.t
	def vertex_center(self,v0,v1):
		[cor0,cor1]=np.array(self.vertex.loc[[v0,v1],list('xyz')])
		return 0.5*(cor0+cor1)
	def update_vertex_edge(self,v):
		self.vertex.loc[v,'edge']=0
		l=self.vertex.loc[v,'l']
		if sum(self.line.loc[l,'edge']==ONEDGE)>0: self.vertex.loc[v,'edge']=ONEDGE
	def vertex_line(self,v,t):
		# function to find the line connected to vertex v in triangle t
		# line=[]
		# for li in self.triangle.loc[t,['l0','l1','l2']]:
		# 	if self.line.loc[li,'v0']==v or self.line.loc[li,'v1']==v:
		# 		line=line+[li]
		return np.intersect1d(self.triangle.loc[t,['l0','l1','l2']],self.vertex.loc[v,'l'])
	def line_side_tl(self,v,l):
		lp=l      #line pointer
		lsub,tsub=[],[]
		while True:
			lsub=lsub+[lp]
			t=self.line.loc[lp,'t']
			tsub=tsub+[t]
			lp=sum(self.vertex_line(v,t))-lp
			lsub=lsub+[lp]
			ledge=self.line.loc[lp,'edge']
			if ledge!=ONEDGE: 
				lp=np.abs(ledge)
			else:
				break
		return lsub,tsub
	def set_line_edge(self,l,lindex):
		#set line edge infomation:
		#32767: line is on edge
		#i(not 32767): line is paired with line |i|
		self.line.loc[l,'edge']=lindex
	def link_line(self,l0,l1):
		self.set_line_edge(l0,l1)
		self.set_line_edge(l1,-l0)
	def unlink_line(self,l0,l1):
		self.set_line_edge(l0,ONEDGE)
		self.set_line_edge(l1,ONEDGE)
	def update_tria_edge(self,t):
		#triangle is on edge as long as there is one line on edge
		#32767:on edge
		#0: not on edge
		self.triangle.loc[t,'edge']=0
		l=self.triangle.loc[t,['l0','l1','l2']].values
		if sum(self.line.loc[l,'edge']==ONEDGE)>0: self.triangle.loc[t,'edge']=ONEDGE
	def tri_norm_vector(self,t):
		#get triangle's norm vector
		[v0,v1,v2]=self.triangle.loc[t,['v0','v1','v2']]
		[cor0,cor1,cor2]=np.array(self.vertex.loc[[v0,v1,v2],list('xyz')])
		return np.cross(cor0,cor1)+np.cross(cor1,cor2)+np.cross(cor2,cor0)
	# *********************Open Angle functions***************************
	def vertex_open_angle(self,v):
		vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
		nl=self.tri_norm_vector(tl)
		nr=self.tri_norm_vector(tr)
		[vcor,vlcor,vrcor]=np.array(self.vertex.loc[[v,vl,vr],list('xyz')])
		al=vlcor-vcor
		ar=vrcor-vcor
		n_ave=nl+nr
		n_new=np.cross(al,ar)
		al_hat=al/np.linalg.norm(al)
		ar_hat=ar/np.linalg.norm(ar)
		cos_theta=np.dot(al_hat,ar_hat)
		if (cos_theta>1.) and (cos_theta<1.+1.e-8): cos_theta=1.
		if (cos_theta<-1.) and (cos_theta>-1.-1.e-8): cos_theta=-1.
		angle=np.arccos(cos_theta)
		if np.dot(n_ave,n_new)>=-1.0E-12: return angle
		elif len(self.vertex.loc[v,'t'])>4:
			#print('two triangle intersect,angle rule changed')
			return -angle
		#three neighbors are easily perpenticuler with norm vector, define always as smaller angle(<180degree)
		elif len(self.vertex.loc[v,'t'])==3:
			return angle
		else: 
			#angle is large that left*right point to inside(opposie as n)
			return 2*np.pi-angle
	# *********************Test functions***************************
	def test_line_attach(self,l):
		[v0,v1]=self.line.loc[l,['v0','v1']]
		if (len(self.vertex.loc[v0,'t'])<=4) and (len(self.vertex.loc[v1,'t'])<=4):
			return True
		else: return False
	def test_shell_close(self):
		return sum(self.triangle['edge'])==0
	def test_vertex_close(self,v):
		vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
	    # usually vl=vr when merge at last step
		if (vl!=vr and len(self.vertex.loc[vl,'t'])+len(self.vertex.loc[vr,'t']) > 6):
			return False
		else: 
			return True
	def test_doublev_close(self,v0,v1):
		v0l,l0l,t0l,v0r,l0r,t0r=self.vertex_near_vlt(v0)
		v2=v0l+v0r-v1
		v1l,l1l,t1l,v1r,l1r,t1r=self.vertex_near_vlt(v1)
		v3=v1l+v1r-v0
		if v2 in self.vertex.index.values and v3 in self.vertex.index.values:
			if len(self.vertex.loc[v2,'t'])+len(self.vertex.loc[v3,'t'])<6:
				return True
		else:
			return False
	def test_vertex_merge(self,v):
		vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
	    # usually vl=vr when merge at last step
		if (len(self.vertex.loc[vl,'t'])>=2 and len(self.vertex.loc[vr,'t']) >=2):
			return True
		else: 
			return False
	def test_vertex_insert(self,v):
		vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
		if (vl==vr or len(self.vertex.loc[v,'t'])>=6 or len(self.vertex.loc[vl,'t'])>=6 or len(self.vertex.loc[vr,'t'])>=6):
			return False
		else: 
			return True
	def test_line_tear(self,l):
		v0,v1=self.line.loc[l,['v0','v1']]
		v0edge,v1edge=np.int_(self.vertex.loc[[v0,v1],'edge']==ONEDGE)
		if v0edge+v1edge==1:
			v=v0 if v1edge else v1
			# lpair=np.abs(self.line.loc[l,'edge'])
			# t0,t1=self.line.loc[[l,lpair],'t']
			# if self.triangle.loc[t0,'edge']!=ONEDGE and self.triangle.loc[t1,'edge']!=ONEDGE:
			if len(self.vertex.loc[v,'t'])<5:
				return True
		return False
	# *********************Action functions***************************
	def attach(self,l):
		#shell grow algorithm:
		#attach new triangle to line l
		[v0,v1,t]=self.line.loc[l,['v0','v1','t']]
		# print('attach to line: {} v({},{})'.format(l,self.pdic[v0]+1,self.pdic[v1]+1))
		v_new_index=self.vertex.index[-1]+1
		v_new_cor=self.create_new_cor(l)
		[l_new0,l_new1,l_new2]=self.line.index[-1]+np.int_([1,2,3])
		t_new=self.triangle.index[-1]+1
		self.add_vertex(v_new_index,v_new_cor,ONEDGE)
		self.add_line(l_new0,v0,v_new_index,t_new,ONEDGE)
		self.add_line(l_new1,v_new_index,v1,t_new,ONEDGE)
		self.add_line(l_new2,v1,v0,t_new,ONEDGE)
		self.add_triangle(t_new,v0,v_new_index,v1,l_new0,l_new1,l_new2,ONEDGE)	
		self.link_line(l,l_new2)                 #update l edge
		self.update_tria_edge(t)
	def detach(self,t):
		#shell remove algorithm:
		#remove a triangle
		v0,v1,v2,l0,l1,l2,edge=self.triangle.loc[t]
		#cprint("deleting triangle t{}, vinfo: {}, {}, {} ".format(t,self.pdic[v0]+1,self.pdic[v1]+1,self.pdic[v2]+1),"yellow")
		self.triangle=self.triangle.drop(t)
		for lpair in np.abs(self.line.loc[[l0,l1,l2],'edge']):
			if not(lpair==ONEDGE):
				self.set_line_edge(lpair,ONEDGE)
				self.update_tria_edge(self.line.loc[lpair,'t'])
		self.line=self.line.drop([l0,l1,l2])
		for vi in [v0,v1,v2]:
			self.vertex.loc[vi,'t'].remove(t)
			for li in [l0,l1,l2]:
				try:
					self.vertex.loc[vi,'l'].remove(li)
				except ValueError:
					pass
			if not self.vertex.loc[vi,'t']:
				self.vertex=self.vertex.drop(vi)
				# self.vertex.loc[vi]=np.inf
			else: self.update_vertex_edge(vi)
	def vertex_merge(self,v0,v1):
		#is part of vertex close
		if v0!=v1:
			vkeep=min(v0,v1)
			vdrop=max(v0,v1)
			self.vertex.loc[vkeep,'l'].extend(self.vertex.loc[vdrop,'l'])
			self.vertex.loc[vkeep,'t'].extend(self.vertex.loc[vdrop,'t'])
			self.line[self.line[['v0','v1']]==vdrop]=vkeep
			self.triangle[self.triangle[['v0','v1','v2']]==vdrop]=vkeep
			self.vertex.loc[vkeep,list('xyz')]=self.vertex_center(vkeep,vdrop)
			self.vertex=self.vertex.drop(vdrop)
		else:
			self.update_vertex_edge(v0)
	def vertex_close(self,v):
		#print('vertex {} ({}) closing'.format(v,self.pdic[v]+1))
		vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
		#print('merge (v{}/l{}) and (v{}/l{})'.format(vl,ll,vr,lr))
		self.link_line(ll,lr)
		self.vertex_merge(vl,vr)
		self.update_vertex_edge(v)
		self.update_tria_edge(tl)
		self.update_tria_edge(tr)
		#print('done;')
	def vertex_insert(self,v):
		#print('vertex {} ({}) inserting'.format(v,self.pdic[v]+1))
		vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
		[l_new0,l_new1,l_new2]=self.line.index[-1]+np.int_([1,2,3])
		t_new=self.triangle.index[-1]+1
		#print('vl:{},ll:{},tl:{},vr:{},lr:{},tr:{},l_new0:{},lnew1:{},lnew2:{},tnew:{}'.format(vl,ll,tl,vr,lr,tr,l_new0,l_new1,l_new2,t_new))
		self.add_line(l_new0,v,vl,t_new,ONEDGE)
		self.add_line(l_new1,vl,vr,t_new,ONEDGE)
		self.add_line(l_new2,vr,v,t_new,ONEDGE)
		self.add_triangle(t_new,v,vl,vr,l_new0,l_new1,l_new2,ONEDGE)
		self.link_line(ll,l_new0)
		self.link_line(lr,l_new2)
		self.update_vertex_edge(v)
		self.update_tria_edge(tl)
		self.update_tria_edge(tr)
	def vertex_divorce(self,v,l):
		lsub,tsub=self.line_side_tl(v,l)
		# print("tsub:"," ",tsub,"lsub:",lsub)

		vnew=self.vertex.index[-1]+1
		self.add_vertex(vnew,self.vertex.loc[v,['x','y','z']],ONEDGE)
		
		self.triangle[self.triangle.loc[tsub,['v0','v1','v2']]==v]=vnew
		self.line[self.line.loc[lsub,['v0','v1']]==v]=vnew
		for ti in tsub:
			self.vertex.loc[v,'t'].remove(ti)
			self.vertex.loc[vnew,'t'].append(ti)
		for li in lsub:
			self.vertex.loc[v,'l'].remove(li)
			self.vertex.loc[vnew,'l'].append(li)
		#print("after divorce: vt:",self.vertex.loc[v,'t'],' ',self.vertex.loc[vnew,'t'])
		#print("after divorce: vl:",self.vertex.loc[v,'l'],' ',self.vertex.loc[vnew,'l'])
	def line_tear(self,l):
		#print("tear line:", l)
		v0,v1,lpair=np.abs(self.line.loc[l,['v0','v1','edge']])
		t,tpair=self.line.loc[[l,lpair],'t']
		v0edge,v1edge=(self.vertex.loc[[v0,v1],'edge']==ONEDGE)
		v=v0 if v0edge else v1
		# vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
		#print("tear vertex {} ({})".format(v, self.pdic[v]+1))
		self.unlink_line(l,lpair)
		self.vertex_divorce(v,l)
		self.update_tria_edge(t)
		self.update_tria_edge(tpair)
		self.update_vertex_edge(v0+v1-v)
		# self.shellprint()
	# *********************Hoomd related functions*************************
	def update_particles(self):
		self.n=len(self.vertex)
		# self.pids=self.typeid*self.n
		self.pdic=dict(zip(self.vertex.index, range(self.n)))
		self.particles=np.array(self.vertex[list('xyz')])
		# for trashi in self.vertex[self.vertex.edge==np.inf].index:
		# 	self.pids[trashi-1]=2
		# 	self.vertex.loc[trashi,['x','y','z']]=np.random.rand(3)
	def correct_index(self,obj):
		for i,row in enumerate(obj):
			obj[i]=list(map(lambda x: self.pdic[x], row))
		return np.array(obj)
	def update_bonds(self):
		bonds=np.int_(self.line[['v0','v1']])
		self.bonds=self.correct_index(bonds)
	def update_dihedrals(self):
		dihedralgroup=[]
		for l in self.line[self.line.edge!=ONEDGE].index:
			l_pair=self.line.loc[l,'edge']
			if l_pair>0:
				[v0,v1,t0]=self.line.loc[l,['v0','v1','t']]
				[v10,v11,t1]=self.line.loc[l_pair,['v0','v1','t']]
				v2=self.other_vertex(t0,v0,v1)
				v12=self.other_vertex(t1,v10,v11)
				if len(set([v2,v1,v11,v12]))==len([v2,v1,v11,v12]):
					dihedralgroup=dihedralgroup+[[v2,v1,v11,v12]]
		self.dihedrals=self.correct_index(dihedralgroup)
	def update_shell_info(self):
		# self.update_shell_radius()
		self.update_particles()
		self.update_bonds()
		self.update_dihedrals()
	def copy_hoomd(self,snap):
		list_particle = []
		for i, cor in enumerate(snap.particles.position, start = 0):
			list_particle.append(list(cor))
			if i >= len(self.vertex)-1:
				break
		self.vertex[['x','y','z']]=np.array(list_particle)
	# *********************Shell Info functions*************************
	def pentnum(self):
		pentnum=0
		for vi in self.vertex[self.vertex.edge==0.].index:
			if len(self.vertex.loc[vi,'t'])==5: pentnum+=1
		return pentnum
	def pentenergy(self):
		return self.pentE * self.pentnum()
	def line_tension(self):
		return self.l_tension * sum(self.line.edge==ONEDGE)
	def EperS(self,E):
		return E / len(self.triangle)
	def vertex_Ehp(self,v,merge=False):
		vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
		t0,t1,t2=self.vertex.loc[[v,vl,vr],'t']
		if merge:
			return -(len(t1)*len(t2))*2*self.hp_eps
		else:
			return -(len(t0)+len(t1)+len(t2))*2*self.hp_eps
	def doublev_Ehp(self,v0,v1):
		v0l,l0l,t0l,v0r,l0r,t0r=self.vertex_near_vlt(v0)
		v2=v0l+v0r-v1
		v1l,l1l,t1l,v1r,l1r,t1r=self.vertex_near_vlt(v1)
		v3=v1l+v1r-v0
		t0,t1,t2,t3=self.vertex.loc[[v0,v1,v2,v3],'t']
		return -(len(t0)+len(t1)+len(t2)+len(t3)+len(t2)*len(t3))*2*self.hp_eps
	def line_Ehp(self,l,merge=False):
		v0,v1=self.line.loc[l,['v0','v1']]
		t0,t1=self.vertex.loc[[v0,v1],'t']
		if merge:
			v0edge,v1edge=(self.vertex.loc[[v0,v1],'edge']==ONEDGE)
			v=v0 if v0edge else v1
			t=self.vertex.loc[v,'t']
			lsub,tsub=self.line_side_tl(v,l)
			return (len(t)-len(tsub))*len(tsub)*2*self.hp_eps
		else:
			return -(len(t0)+len(t1))*2*self.hp_eps
	def tri_Ehp(self,ti=-1):
		if ti==-1: ti=self.triangle.index[-1]
		v0,v1,v2=self.triangle.loc[ti,['v0','v1','v2']]
		t0,t1,t2=self.vertex.loc[[v0,v1,v2],'t']
		return -(len(t0)+len(t1)+len(t2)-3)*self.hp_eps
	def shell_Ehp(self):
		neighNlist=self.vertex['t'].apply(lambda x: len(x))
		return -np.dot(neighNlist,neighNlist-1)*self.hp_eps
		# np.array(self.vertex['t'].tolist()).size
	def shellprint(self):
		cprint('Shell Structure:\n','yellow')
		#print(self.pdic)
		print(self.vertex)
		print(self.line)
		print(self.triangle)
		#print(self.R)
	def shell_copy(self):
		return pickle.loads(pickle.dumps(self))#copy.deepcopy(self)
	# *********************Grow functions*************************
	def get_all_grows(self,Ekey=[],iniE=[],grow=3):
		growlist=pd.DataFrame(columns=['choice','growtype','growindex','Ehp']+Ekey)
		for vi in self.vertex[self.vertex.edge==ONEDGE].index:
			vt=len(self.vertex.loc[vi,'t'])
			if vt==5:
				if self.test_vertex_insert(vi):
					growlist.loc[len(growlist)]=['G','v',vi,self.vertex_Ehp(vi)]+iniE
				# if self.test_vertex_merge(vi):
				if self.test_vertex_close(vi):
					growlist.loc[len(growlist)]=['S','vc',vi,self.vertex_Ehp(vi,True)]+iniE
			#vertex move
			if vt==6:
				if self.test_vertex_close(vi):
					growlist.loc[len(growlist)]=['S','vc',vi,self.vertex_Ehp(vi,True)]+iniE
		vindex=list(growlist[growlist.growtype=='v']['growindex'])
		#print("vindex:",vindex)
		for li in self.line[self.line.edge==ONEDGE].index:
			if self.test_line_attach(li):
				growlist.loc[len(growlist)]=['G','l',li,self.line_Ehp(li)]+iniE
			[v0,v1]=self.line.loc[li,['v0','v1']]
			# print("line info:",li,v0,v1,(v0 in vindex) and (v1 in vindex))
			if (v0 in vindex) and (v1 in vindex):
				if self.test_doublev_close(v0,v1):
					growlist.loc[len(growlist)]=['G','vv',li,self.doublev_Ehp(v0,v1)]+iniE
		for li in self.line[self.line.edge!=ONEDGE].index:
			if self.line.loc[li,'edge']>0:
				if self.test_line_tear(li):
					growlist.loc[len(growlist)]=['S','vt',li,self.line_Ehp(li,True)]+iniE
		
		if grow==0: 
			return growlist[growlist['choice']=='S']
		elif grow==1: 
			# growlist=self.grow_filter(growlist)
			return growlist[growlist['choice']=='G']
		elif grow==2:
			# growlist=self.grow_filter(growlist)
			return growlist
		else: 
			return growlist
	def shell_grow(self,growlist,growi):
		stwin=self.shell_copy()
		choice,growtype,growindex=growlist.loc[growi,['choice','growtype','growindex']]
		#cprint('Grow {}, choice:{}, growtype:{}, growindex:{}'.format(growi,choice,growtype,growindex), 'yellow')
		if growtype=='v':
			stwin.vertex_insert(growindex)
		elif growtype=='vv':
			v0,v1=stwin.line.loc[growindex,['v0','v1']]
			stwin.attach(growindex)
			stwin.vertex_close(v0)
			stwin.vertex_close(v1)
		elif growtype=='vc':
			stwin.vertex_close(growindex)
		elif growtype=='vt':
			stwin.line_tear(growindex)
		elif growtype=='l':
			stwin.attach(growindex)
		else:
			cprint('unknown growtype!','red')
		return stwin
	def detach_tlist(self):
		tlist=[]
		for ti,row in self.triangle[self.triangle.edge==ONEDGE].iterrows():
			[v0,v1,v2,l0,l1,l2,tedge]=row
			vedgenum=sum(np.int_(self.vertex.loc[[v0,v1,v2],'edge'])==ONEDGE)
			ledgenum=sum(np.int_(self.line.loc[[l0,l1,l2],'edge'])==ONEDGE)
			# if (self.tri_open_angle(ti)>self.theta_tri):
			if (vedgenum==2 and ledgenum==1) or (vedgenum==3 and ledgenum==2):
					tlist+=[ti]
		return tlist
	def diffuse_tlist(self):
		tlist=[]
		for ti,row in self.triangle[self.triangle.edge==ONEDGE].iterrows():
			[v0,v1,v2,l0,l1,l2,tedge]=row
			vedgenum=sum(np.int_(self.vertex.loc[[v0,v1,v2],'edge'])==ONEDGE)
			ledgenum=sum(np.int_(self.line.loc[[l0,l1,l2],'edge'])==ONEDGE)
			# if (self.tri_open_angle(ti)>self.theta_tri):
			if (vedgenum==2 and ledgenum==1) or (vedgenum==3 and ledgenum==2):
				if self.tri_Ehp(ti)>=self.thres_Ehp:
					tlist+=[ti]
		return tlist
	def shell_detach(self,ti):
		stwin=self.shell_copy()
		stwin.detach(ti)
		return stwin
	def shell_merge(self):
		status=0
		for vi in self.vertex[self.vertex.edge==ONEDGE].index:
			if (np.abs(self.vertex_open_angle(vi))<self.theta_lower) and self.test_vertex_close(vi):
				self.vertex_close(vi)
				status=1
				break
		if sum(self.triangle['edge'])==0.0:status=-1
		return status
	######################################################
	#####################################################
	######################################################
	def Saeed_detach_tlist(self):
		detachlist=pd.DataFrame(columns=['choice','detachtype','detachindex','weights'])
		for ti,row in self.triangle[self.triangle.edge==ONEDGE].iterrows():
			[v0,v1,v2,l0,l1,l2,tedge]=row
			vedgenum=sum(np.int_(self.vertex.loc[[v0,v1,v2],'edge'])==ONEDGE)
			ledgenum=sum(np.int_(self.line.loc[[l0,l1,l2],'edge'])==ONEDGE)
			if (vedgenum==2 and ledgenum==1):
				detachlist.loc[len(detachlist)]=['D','v',ti,1]
			elif (vedgenum==3 and ledgenum==2):
				detachlist.loc[len(detachlist)]=['DD','vv',ti,100]
		return detachlist
	def pent_vertex_merge(self,v0,v1):
		#is part of vertex close
		if v0!=v1:
			vkeep=min(v0,v1)
			vdrop=max(v0,v1)
			self.vertex.loc[vkeep,'l'].extend(self.vertex.loc[vdrop,'l'])
			self.vertex.loc[vkeep,'t'].extend(self.vertex.loc[vdrop,'t'])
			self.line[self.line[['v0','v1']]==vdrop]=vkeep
			
			self.triangle[self.triangle[['v0','v1','v2']]==vdrop]=vkeep
			
			self.vertex.loc[vkeep,list('xyz')]=self.vertex_center(vkeep,vdrop)
			if vdrop in self.vertex.index:
				self.vertex=self.vertex.drop(vdrop)
		else:
			self.update_vertex_edge(v0)
	def pentamer_vertex_close(self,v):
		vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
		self.link_line(ll,lr)	
		self.pent_vertex_merge(vl,vr)	
		self.update_vertex_edge(v)
		self.update_tria_edge(tl)
		self.update_tria_edge(tr)
		#print('done;')
	def impose_pentamer_merge(self):
		for vi in self.vertex[self.vertex.edge==ONEDGE].index:
			if vi in self.vertex.index:
				if len(self.vertex.loc[vi,'t'])==5 and (np.abs(self.vertex_open_angle(vi))<self.theta_lower):
					print('vetex closed')
					self.pentamer_vertex_close(vi)
# *********************Outside Class functions*************************
def create_shell(param):
	vertex=pd.DataFrame(columns=['x','y','z','l','t','edge'])
	line=pd.DataFrame(columns=['v0','v1','t','edge']).astype('int32')
	triangle=pd.DataFrame(columns=['v0','v1','v2','l0','l1','l2','edge']).astype('int32')
	return shell(vertex,line,triangle,param)
# ***************************** Initiate vertex line triangle From CSV *************** #

def get_footer(file_, all_row = 0, vertex_start = 1, line_start = 6, tri_start = 11,frame_end=14):
		with open(file_) as f:
			csv_reader = csv.reader(f,delimiter =',')
			# Findind all row and starting point of ver, line, tri,
			num_of_hat = 0
			for index, row in enumerate(csv_reader):
				#print(index,row)
				if frame_end < index:
					if num_of_hat < 4:
						if row[0][0:1] == '#':
							num_of_hat += 1 
							all_row = index
							if row[0][0:4] == '#ver':
								vertex_start = index
							elif row[0][0:5] == '#line':
								line_start = index
							elif row[0][0:5] == '#tria':
								tri_start = index
							elif tri_start < index:
								frame_end = index
			frame_size = index
			f.seek(0)
			vertex_value = []
			line_value = []
			triangle_value = []
			num_v = 0
			num_l = 0
			num_t = 0
			for index, row in enumerate(csv_reader):
				# For Vertex#
				dummy_value = []		
				if vertex_start +1 < index < line_start:
					res_line = json.loads(row[4]) # Converting string into list
					res_tri = json.loads(row[5])
					dummy_value = [int(row[0]),float(row[1]),float(row[2]),float(row[3]),res_line,res_tri,int(row[6])]
					vertex_value.insert(num_v,dummy_value)
					num_v += 1
				# For Line#
				dummy_value = []
				if line_start+1 < index < tri_start:
					dummy_value = [int(row[0]),int(row[1]),int(row[2]),int(row[3]),int(row[4])]
					line_value.insert(num_l,dummy_value)
					num_l += 1
				# For triangle#
				dummy_value = []
				if tri_start+1 < index < frame_end:
					num_t += 1
					dummy_value = [int(row[0]),int(row[1]),int(row[2]),int(row[3]),int(row[4]),int(row[5]),int(row[6]),int(row[7])]
					triangle_value.insert(num_t,dummy_value)
				
			return vertex_value, line_value, triangle_value,frame_end,frame_size

'''
	def tri_centers(self):
		cors=[]
		for ti,row in self.triangle.iterrows():
			[v0,v1,v2]=row[:3]
			cors=cors+[np.average(self.vertex.loc[[v0,v1,v2],list('xyz')],axis=0)]
		return cors
'''
'''
	def grow_filter(self,growlist,filterl=True,filterv=True):
		tricenters=self.tri_centers()
		if filterl:
			for gi,row in growlist[growlist['growtype']=='l'].iterrows():
				if not self.test_line_attach_update(row['growindex'],tricenters):
					growlist.loc[gi,'choice']='F'
		if filterv:
			for gi,row in growlist[growlist['growtype']=='v'].iterrows():
				if not self.test_vertex_insert_update(row['growindex'],tricenters):
					growlist.loc[gi,'choice']='F'
		# print("refined growlist:",growlist)
		return growlist[growlist.choice!='F']
'''
'''
	def line_open_angle(self,l):
		v0,v1=self.line.loc[l,['v0','v1']]
		return min(self.vertex_open_angle(v0),self.vertex_open_angle(v1))
'''
'''
	def tri_open_angle(self,t):
		angle=0.
		for v in self.triangle.loc[t,['v0','v1','v2']]:
			if self.vertex.loc[v,'edge']==ONEDGE:
				angle+=self.vertex_open_angle(v)
		return angle
'''
'''
	def test_line_attach_update(self,l,tricenter):
		newp=self.create_new_cor(l=l,theta=self.theta0)
		[v0,v1]=self.line.loc[l,['v0','v1']]
		[cor0,cor1]=np.array(self.vertex.loc[[v0,v1],list('xyz')])
		newt_c=(newp+cor0+cor1)/3.
		# print("edge triangle centers:",cors)
		centdist=np.array(list(map(lambda x: np.linalg.norm(newt_c-x),tricenter)))
		pdist1=np.array(list(map(lambda x: np.linalg.norm(newt_c-x),np.array(self.vertex[list('xyz')]))))
		pdist2=np.array(list(map(lambda x: np.linalg.norm(newp-x),tricenter)))
		# print("centdist:",centdist)
		# print("pdis:",pdist)
		if np.any(centdist<2*self.rc) or np.any(pdist1<self.rp+self.rc) or np.any(pdist2<self.rp+self.rc):
			return False
		else: return True
	
	def test_vertex_insert_update(self,v,tricenter):
		vl,ll,tl,vr,lr,tr=self.vertex_near_vlt(v)
		[cor0,cor1,cor2]=np.array(self.vertex.loc[[v,vl,vr],list('xyz')])
		newt_c=(cor0+cor1+cor2)/3.
		centdist=np.array(list(map(lambda x: np.linalg.norm(newt_c-x),tricenter)))
		pdist=np.array(list(map(lambda x: np.linalg.norm(newt_c-x),np.array(self.vertex[list('xyz')]))))
		# print("test_insert centdist{}\npdist:{}".format(centdist,pdist))
		if np.any(centdist<2*self.rc) or np.any(pdist<self.rp+self.rc):
			return False
		else: return True
'''
'''
	def dihedral_angles(self):
		angles=[]
		for l in self.line[self.line.edge!=ONEDGE].index:
			l_pair=self.line.loc[l,'edge']
			if l_pair>0:
				[v0,v1,t0]=self.line.loc[l,['v0','v1','t']]
				[v10,v11,t1]=self.line.loc[l_pair,['v0','v1','t']]
				n0=self.tri_norm_vector(t0)
				n1=self.tri_norm_vector(t1)
				[cor0,cor1]=np.array(self.vertex.loc[[v0,v1],list('xyz')])
				vecl=cor1-cor0
				dotn=np.dot(n0,n1)
				crossn=np.cross(n0,n1)
				theta=np.arccos(dotn/np.linalg.norm(n0)/np.linalg.norm(n1))
				if np.dot(crossn,vecl)<0:
					theta*=-1
				angles=angles+[theta]
		#cprint("dihedral angles:{}".format(angles),'yellow')
	
	def shell_radius(self):
		vcor=self.vertex[self.vertex.edge!=np.inf][list('xyz')]
		center=vcor.sum()/len(vcor)
		return np.mean(np.linalg.norm(vcor-center,axis=1))
	
	def bondenergy(self):
		return self.bondE * len(self.dihedrals)
'''
'''
	def get_all_detach(self):
		detachlist=pd.DataFrame(columns=['choice','detachtype','detachindex','angle']+Ekey)
		if len(self.triangle)>2:
			for ti,row in self.triangle[self.triangle.edge==ONEDGE].iterrows():
				[v0,v1,v2,l0,l1,l2,tedge]=row
				vedgenum=sum(np.int_(self.vertex.loc[[v0,v1,v2],'edge'])==ONEDGE)
				ledgenum=sum(np.int_(self.line.loc[[l0,l1,l2],'edge'])==ONEDGE)
				tEhp=self.tri_Ehp(ti)
				if tEhp>self.thres_Ehp:
					if (vedgenum==2 and ledgenum==1):
						growlist.loc[len(growlist)]=['D','tv',ti,tEhp]+iniE
					elif (vedgenum==3 and ledgenum==2):
						growlist.loc[len(growlist)]=['D','tl',ti,tEhp]+iniE
		return detachlist
'''