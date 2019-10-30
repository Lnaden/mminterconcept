import matplotlib.pylab as plt
plt.rc('text', usetex=True)
plt.rcParams['figure.figsize'] = (16, 18)

def plot_com(time_mdt, com_mdt, time_mda, com_mda):
        plt.plot(time_mdt, com_mdt, '-o')
        plt.plot(time_mda, com_mda, '-*')
        plt.xlabel('Time (ps)')
        plt.ylabel('COM (nm)')
        plt.grid(linestyle=':')
        plt.legend(['MDTraj', 'MDAnalysis'])
        plt.show()

def plot_rmsd(time_mdt, rmsd_mdt, time_mda, rmsd_mda):
        plt.plot(time_mdt, rmsd_mdt, '-o')
        plt.plot(time_mda, rmsd_mda, '-*')
        plt.xlabel('Time (ps)')
        plt.ylabel('RMSD (nm)')
        plt.grid(linestyle=':')
        plt.legend(['MDTraj', 'MDAnalysis'])
        plt.show()

def plot_rdf(dist_mdt, rdf_mdt, dist_mda, rdf_mda):
        plt.plot(dist_mdt, rdf_mdt, '--')
        plt.plot(dist_mda[1:], rdf_mda[1:], '-.')
        plt.xlabel('Distance (nm)')
        plt.ylabel('RDF')
        plt.grid(linestyle=':')
        plt.legend(['MDTraj', 'MDAnalysis'])
        plt.show()

def plot_den(time_mdt, den_mdt, time_mda, den_mda):
        plt.plot(time_mdt, den_mdt, '-o')
        plt.plot(time_mda, den_mda, '-*')
        plt.xlabel('Time (ps)', fontsize=16)
        plt.ylabel(r'Density ($Kg/m^3$)', fontsize=16)
        plt.grid(linestyle=':')
        plt.legend(['MDTraj', 'MDAnalysis'])
        plt.show()
