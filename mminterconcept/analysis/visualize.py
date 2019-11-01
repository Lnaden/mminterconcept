from ..construction.tools import ImportPDB
import nglview
import mdtraj
import numpy

import matplotlib.pylab as plt
plt.rc('text', usetex=True)

def plot_com(time_mdt, com_mdt, time_mda, com_mda):

	ncomp = numpy.size(com_mdt,1)
	fig, ax = plt.subplots(ncomp,1)

	fig.set_figheight(11)
	fig.set_figwidth(5)
	fig.set_dpi(120)

	comp = ('x', 'y', 'z')

	for i in range(ncomp):
		ax[i].plot(time_mdt, com_mdt[:,i], '-o')
		ax[i].plot(time_mda, com_mda[:,i], '-*')
		ax[i].set_xlabel('Time (ps)')
		ax[i].set_ylabel(f'COM-{comp[i]} (nm)')
		ax[i].grid(linestyle=':')
		ax[i].legend(['MDTraj', 'MDAnalysis'])

	plt.show()

def plot_rmsd(time_mdt, rmsd_mdt, time_mda, rmsd_mda):

	fig = plt.figure(figsize=(5,4), dpi=120)
	plt.plot(time_mdt, rmsd_mdt, '-o')
	plt.plot(time_mda, rmsd_mda, '-*')
	plt.xlabel('Time (ps)')
	plt.ylabel('RMSD (nm)')
	plt.grid(linestyle=':')
	plt.legend(['MDTraj', 'MDAnalysis'])
	plt.show()

def plot_rdf(dist_mdt, rdf_mdt, dist_mda, rdf_mda):
        fig = plt.figure(figsize=(5,4), dpi=120)
        plt.plot(dist_mdt, rdf_mdt, '--')
        plt.plot(dist_mda[1:], rdf_mda[1:], '-.')
        plt.xlabel('Distance (nm)')
        plt.ylabel('RDF')
        plt.grid(linestyle=':')
        plt.legend(['MDTraj', 'MDAnalysis'])
        plt.show()

def plot_den(time_mdt, den_mdt, time_mda, den_mda):
        fig = plt.figure(figsize=(5,4), dpi=120)
        plt.plot(time_mdt, den_mdt, '-o')
        plt.plot(time_mda, den_mda, '-*')
        plt.xlabel('Time (ps)')
        plt.ylabel(r'Density ($Kg/m^3$)')
        plt.grid(linestyle=':')
        plt.legend(['MDTraj', 'MDAnalysis'])
        plt.show()

def plot_rg(time_mdt, rg_mdt, time_mda, rg_mda):
        fig = plt.figure(figsize=(5,4), dpi=120)
        plt.plot(time_mdt, rg_mdt, '-o')
        plt.plot(time_mda, rg_mda, '-*')
        plt.xlabel('Time (ps)')
        plt.ylabel('RG (nm)')
        plt.grid(linestyle=':')
        plt.legend(['MDTraj', 'MDAnalysis'])
        plt.show()

def ngl_view(**args):

	if 'pdbFile' in args:
		traj = mdtraj.load(args['pdbFile'])
		view = nglview.show_mdtraj(traj)

	elif 'pdbID' in args:
		with ImportPDB(args['pdbID']) as traj:
    			view = nglview.show_mdtraj(traj)

	return view
