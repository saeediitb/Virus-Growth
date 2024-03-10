import hoomd
from hoomd import md
import numpy as np
import pandas as pd
import random
import json
import copy,sys
from virus.capsid import create_shell
from . import kmc_system
from ast import literal_eval
from termcolor import colored, cprint
import pickle 
from virus.IOframe import *

ONEDGE=32767
NORLX = 123456
class System(object): #why except object word it is not accepting   
    def __init__(self, param):
        self.shell = create_shell(param)
        self.param = param
        self.theta0 = param['theta0']
        self.time = 0.0
        self.nl = hoomd.md.nlist.Tree(buffer = 1.0, exclusions = ('angle','bond','dihedral','special_pair'))
    def initial(self):
        self.shell.initial()
        self.sim = hoomd.Simulation(device = hoomd.device.CPU(), seed = 1)
        #self.shell.capsid_shape()
        self.shell.update_shell_info()
    def Capsid_tria_remove(self):
        ti_row = self.shell.triangle.sample()
        ti = ti_row.index.values.astype(int)[0]
        self.shell = self.shell.shell_detach(ti) 
        self.try_relax()
        frame_out(self.time, self.shell)
        energy_frame(self.time, self.shellE, col_name=False)
    def setup_potential(self):
        if self.param['MIE']:
            nn,mm=6,1
            self.MIE = hoomd.md.pair.Mie(nlist = self.nl,mode='shift')
            self.MIE.params[('A','A')] = dict(sigma = 2.0*self.param['R0'], epsilon = self.param['lj_eps'],n=nn,m=mm)
            self.MIE.r_cut[('A','A')] = 2.0*(pow(nn/mm,1/(nn-mm)))*self.param['R0']
            self.MIE.params[('A','B')] = dict(sigma = 2.0*self.param['R0'], epsilon = self.param['lj_eps'],n=nn,m=mm)
            self.MIE.r_cut[('A','B')] = 2.0*(pow(nn/mm,1/(nn-mm)))*self.param['R0']
            self.MIE.params[('B','B')] = dict(sigma = 2.0*self.param['R0'], epsilon = self.param['lj_eps'],n=nn,m=mm)
            self.MIE.r_cut[('B','B')] = 2.0*(pow(nn/mm,1/(nn-mm)))*self.param['R0']
        if self.param['BONDHARM']: #True
            self.harmonic_force = hoomd.md.bond.Harmonic()
            self.harmonic_force.params['out_ver_to_out_ver'] = dict(k = self.param['sh_bondk'], r0 = self.param['l0'])
        if self.param['DIHEDHARM']:
            self.dihedral_force = hoomd.md.dihedral.Harmonic()
            self.dihedral_force.params['capsomer'] = dict(k = 2*self.param['sh_dihedk'], d = -1, n = 1, phi0 = np.pi-self.theta0)
    def setup_operation(self, file = None):
         # GSD
        if file is None: file = self.param['gsdfile']
        gsd_writer = hoomd.write.GSD(filename = file, filter=hoomd.filter.All(), trigger=hoomd.trigger.Periodic(10), mode='ab', dynamic=['property', 'attribute', 'topology'])
        self.sim.operations.writers.append(gsd_writer)
        thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(
            filter=hoomd.filter.All())
        self.sim.operations.computes.append(thermodynamic_properties)
        log_rec = np.array(['num_particles','potential_energy','kinetic_energy']) # quantities to be extracted
        self.logger = hoomd.logging.Logger(categories=['scalar', 'string'])
        self.logger.add(thermodynamic_properties, quantities=log_rec)
        self.logger.add(self.sim, quantities=['timestep', 'walltime'])
        # self.logger.add(self.lj, quantities=['energy'])
        self.logger.add(self.harmonic_force, quantities=['energy'])
        self.logger.add(self.dihedral_force, quantities=['energy'])
        self.logger.add(self.MIE, quantities=['energy'])
    def setup_integrator(self):
        # Integrator
        self.fire = md.minimize.FIRE(dt = self.param['dt'],
                        force_tol = self.param['ftol'],
                        angmom_tol = 1e-2,
                        energy_tol = self.param['Etol'],
                        forces = [self.harmonic_force, self.dihedral_force, self.MIE])
        self.fire.methods.append(md.methods.NVE(hoomd.filter.All()))
        self.sim.operations.integrator = self.fire
    def try_relax(self):
        num = 0
        time_list = [0.005,0.001,0.0005,0.0001,0.00005,0.00001]
        while num < len(time_list):
            self.Dt = time_list[num]
            try:
                max_run = self.relax()
                if max_run == self.param['maxrun']:
                    return self.param['maxrun']
                else:
                    return max_run
            except RuntimeError:
                print('Particles are going very far:',self.Dt)
                num += 1
        if num == len(time_list):
            print('System does not even relaxed at:',self.Dt)
            return NORLX # For not relaxing at dt
    def relax(self):
        snap = self.sim.state.get_snapshot()
        self.shell.update_shell_info()
        Update_Snap(snap, self.shell, self.param['R'], dihedral=True)
        self.sim.state.set_snapshot(snap)
        self.fire.dt = self.Dt
        count = 0
        while not(self.fire.converged) and count<self.param['maxrun']:
            self.fire.dt = self.Dt
            self.sim.run(10)#1*self.param['total'])
            count+=1
        self.fire.reset()
        snap = self.sim.state.get_snapshot()
        self.shell.copy_hoomd(snap)
        self.shellE = self.energy()
        if count>=self.param['maxrun']:
            return self.param['maxrun']
        else:
            return count
    def energy(self):
        logger = self.logger
        Ekey = self.param['growrec']
        shell = pickle.loads(pickle.dumps(self.shell))
        Elist = pd.DataFrame(np.zeros((1,len(Ekey))), columns = Ekey)
        Energy_tol = logger.log()['md']['compute']['ThermodynamicQuantities']['potential_energy'][0]
        Stretch = logger.log()['md']['bond']['Harmonic']['energy'][0]
        Bending = logger.log()['md']['dihedral']['Harmonic']['energy'][0]
        #Lj_energy = logger.log()['md']['pair']['LJ0804']['energy'][0]
        Lj_energy = logger.log()['md']['pair']['Mie']['energy'][0]
        Eltension = shell.line_tension()
        Epent = shell.pentenergy()
        Ehp = shell.shell_Ehp()
        for i, key in enumerate(Ekey):
            if key == 'tot': Elist[key] = Energy_tol + Eltension + Epent + Ehp
            elif key == 'Eshell': Elist[key] = Energy_tol + Eltension + Epent + Ehp
            elif key == 'Stretching': Elist[key] = Stretch
            elif key == "Bending": Elist[key] = Bending
            elif key == "lj": Elist[key] = Lj_energy 
            elif key == "Sub_units": Elist[key] = len(self.shell.triangle)
            elif key == "hp": Elist[key] = Ehp
        return Elist
    def printenergy(self,status):
        cprint('shellE:\n{}\nSubunits:{}\nEshell:{}\nEshell/Subunits:{}\nElasticity:{}\nElasticity/Subunits:{}'.format(self.shellE,len(self.shell.triangle),self.shellE.loc[0,'Eshell'],self.shellE.loc[0,'Eshell']/len(self.shell.triangle),self.shellE.loc[0,'Stretching'],self.shellE.loc[0,'Stretching']/len(self.shell.triangle)),'green',attrs=['reverse'] if status else [])
        cprint('Current loop is running for:\nks:   {};  ep_hp:    {}' .format(self.param["sh_bondk"],self.shell.hp_eps))
        if status: self.shell.shellprint()
    def mc_info(self,shell):
        growlist = shell.get_all_grows(grow = 2)
        #cprint("\n\n\nMC system growlist:\n{}\n\diffusetlist:{}".format(growlist,diffusetlist),'green')
        growvlist = growlist[growlist['choice'] == 'S']
        return growvlist
    #------------------KMC---step----system----########
    def kmc_moves(self):
        return kmc_system.kmc_moves(self)
    def kmc_vmove(self, shell, shellE, event):
        return kmc_system.kmc_vmove(self, shell, shellE, event)
    def kmc_tmove(self,shell, shellE, event,a_i):
        return kmc_system.kmc_tmove(self, shell, shellE, event,a_i)
    def kmc_grow(self, shell, shellE, event):
        return kmc_system.kmc_grow(self, shell, shellE, event)
    def kmc_detach(self, shell, shellE, event):
        return kmc_system.kmc_detach(self, shell, shellE, event)
def assign_snap_particles(snap, obj, R):
    pend = obj.n
    snap.particles.N = pend+1
    l0 = list(np.array(obj.vertex[list('xyz')])) # outer triangle vertex
    total_particle = np.array(l0)
    core_pos = Core_position(total_particle[:3, :3], R)
    for i, cor in enumerate(total_particle):
        snap.particles.diameter[i] = 2*obj.R0
        snap.particles.position[i] = cor
        snap.particles.typeid[i] = 0
    # Core
    snap.particles.diameter[i+1] = 2*R
    snap.particles.position[i+1] = core_pos
    snap.particles.typeid[i+1] = 1
def assign_snap_bonds(snap, obj):
    l0 = list(obj.bonds)
    snap.bonds.N = len(obj.bonds) 
    totl_bonds=np.array(l0)
    for i, bond in enumerate(totl_bonds):
        snap.bonds.group[i] = bond
        snap.bonds.typeid[i] = 0
def assign_snap_dihedrals(snap, obj):
    snap.dihedrals.N = len(obj.dihedrals)
    for i, dihedral in enumerate(obj.dihedrals):
        snap.dihedrals.group[i] = dihedral
        snap.dihedrals.typeid[i] = obj.dihedralid
def Update_Snap(snap, obj, R, particle = True, bond = True, dihedral = False):
    #!!assign bond and dihedral before particle
    #particle number needed for re-index bond and dihedral
    if particle:
        assign_snap_particles(snap, obj, R)
    if bond:
        assign_snap_bonds(snap, obj)
    if dihedral:
        assign_snap_dihedrals(snap, obj)

def Core_position(triangle,R):
    centroid = np.mean(triangle, axis=0)
    AB = triangle[1] - triangle[0]
    AC = triangle[2] - triangle[0]
    normal_vector = np.cross(AB, AC)
    normalized_normal = normal_vector / np.linalg.norm(normal_vector)
    distance_from_centroid = R
    core_centre = centroid + distance_from_centroid * normalized_normal
    return core_centre