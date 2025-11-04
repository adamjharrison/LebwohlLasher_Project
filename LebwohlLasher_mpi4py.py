"""
Basic Python Lebwohl-Lasher code.  Based on the paper 
P.A. Lebwohl and G. Lasher, Phys. Rev. A, 6, 426-429 (1972).
This version in 2D.

Run at the command line by typing:

python LebwohlLasher.py <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG>

where:
  ITERATIONS = number of Monte Carlo steps, where 1MCS is when each cell
      has attempted a change once on average (i.e. SIZE*SIZE attempts)
  SIZE = side length of square lattice
  TEMPERATURE = reduced temperature in range 0.0 - 2.0.
  PLOTFLAG = 0 for no plot, 1 for energy plot and 2 for angle plot.
  
The initial configuration is set at random. The boundaries
are periodic throughout the simulation.  During the
time-stepping, an array containing two domains is used; these
domains alternate between old data and new data.

SH 16-Oct-23
"""

import sys
import time
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpi4py import MPI

#=======================================================================
def initdat(nmax):
    """
    Arguments:
      nmax (int) = size of lattice to create (nmax,nmax).
    Description:
      Function to create and initialise the main data array that holds
      the lattice.  Will return a square lattice (size nmax x nmax)
	  initialised with random orientations in the range [0,2pi].
	Returns:
	  arr (float(nmax,nmax)) = array to hold lattice.
    """
    arr = np.random.random_sample((nmax,nmax))*2.0*np.pi
    return arr
#=======================================================================
def plotdat(arr,pflag,nmax):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  pflag (int) = parameter to control plotting;
      nmax (int) = side length of square lattice.
    Description:
      Function to make a pretty plot of the data array.  Makes use of the
      quiver plot style in matplotlib.  Use pflag to control style:
        pflag = 0 for no plot (for scripted operation);
        pflag = 1 for energy plot;
        pflag = 2 for angles plot;
        pflag = 3 for black plot.
	  The angles plot uses a cyclic color map representing the range from
	  0 to pi.  The energy plot is normalised to the energy range of the
	  current frame.
	Returns:
      NULL
    """
    if pflag==0:
        return
    u = np.cos(arr)
    v = np.sin(arr)
    x = np.arange(nmax)
    y = np.arange(nmax)
    cols = np.zeros((nmax,nmax))
    if pflag==1: # colour the arrows according to energy
        mpl.rc('image', cmap='rainbow')
        for i in range(nmax):
            for j in range(nmax):
                cols[i,j] = one_energy(arr,i,j,nmax,rows=0)
        norm = plt.Normalize(cols.min(), cols.max())
    elif pflag==2: # colour the arrows according to angle
        mpl.rc('image', cmap='hsv')
        cols = arr%np.pi
        norm = plt.Normalize(vmin=0, vmax=np.pi)
    else:
        mpl.rc('image', cmap='gist_gray')
        cols = np.zeros_like(arr)
        norm = plt.Normalize(vmin=0, vmax=1)

    quiveropts = dict(headlength=0,pivot='middle',headwidth=1,scale=1.1*nmax)
    fig, ax = plt.subplots()
    q = ax.quiver(x, y, u, v, cols,norm=norm, **quiveropts)
    ax.set_aspect('equal')
    plt.show()
#=======================================================================
def savedat(arr,nsteps,Ts,runtime,ratio,energy,order,nmax):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  nsteps (int) = number of Monte Carlo steps (MCS) performed;
	  Ts (float) = reduced temperature (range 0 to 2);
	  ratio (float(nsteps)) = array of acceptance ratios per MCS;
	  energy (float(nsteps)) = array of reduced energies per MCS;
	  order (float(nsteps)) = array of order parameters per MCS;
      nmax (int) = side length of square lattice to simulated.
    Description:
      Function to save the energy, order and acceptance ratio
      per Monte Carlo step to text file.  Also saves run data in the
      header.  Filenames are generated automatically based on
      date and time at beginning of execution.
	Returns:
	  NULL
    """
    # Create filename based on current date and time.
    current_datetime = datetime.datetime.now().strftime("%a-%d-%b-%Y-at-%I-%M-%S%p")
    filename = "LL-Output-{:s}.txt".format(current_datetime)
    FileOut = open(filename,"w")
    # Write a header with run parameters
    print("#=====================================================",file=FileOut)
    print("# File created:        {:s}".format(current_datetime),file=FileOut)
    print("# Size of lattice:     {:d}x{:d}".format(nmax,nmax),file=FileOut)
    print("# Number of MC steps:  {:d}".format(nsteps),file=FileOut)
    print("# Reduced temperature: {:5.3f}".format(Ts),file=FileOut)
    print("# Run time (s):        {:8.6f}".format(runtime),file=FileOut)
    print("#=====================================================",file=FileOut)
    print("# MC step:  Ratio:     Energy:   Order:",file=FileOut)
    print("#=====================================================",file=FileOut)
    # Write the columns of data
    for i in range(nsteps+1):
        print("   {:05d}    {:6.4f} {:12.4f}  {:6.4f} ".format(i,ratio[i],energy[i],order[i]),file=FileOut)
    FileOut.close()
#=======================================================================
def one_energy(arr,ix,iy,nmax,neighbours=None,rows=None):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  ix (int) = x lattice coordinate of cell;
	  iy (int) = y lattice coordinate of cell;
      nmax (int) = side length of square lattice.
    Description:
      Function that computes the energy of a single cell of the
      lattice taking into account periodic boundaries.  Working with
      reduced energy (U/epsilon), equivalent to setting epsilon=1 in
      equation (1) in the project notes.
	Returns:
	  en (float) = reduced energy of cell.
    """
    en = 0.0
    iyp = (iy+1)%nmax # neighbour with wraparound
    iym = (iy-1)%nmax #
#
# Add together the 4 neighbour contributions
# to the energy
#
    ang = arr[ix,iy]-arr[ix,iyp]
    en += 0.5*(1.0 - 3.0*np.cos(ang)**2)
    ang = arr[ix,iy]-arr[ix,iym]
    en += 0.5*(1.0 - 3.0*np.cos(ang)**2)

    if (rows > 0) and (rows<nmax):
      if ix==0:
        ang = arr[ix,iy]-neighbours[0,iy]
      else:
        ang = arr[ix,iy]-arr[ix-1,iy]
      en += 0.5*(1.0 - 3.0*np.cos(ang)**2)
      
      if ix==rows-1:
        ang = arr[ix,iy]-neighbours[1,iy]
      else:
        ang = arr[ix,iy]-arr[ix+1,iy]
      en += 0.5*(1.0 - 3.0*np.cos(ang)**2)
    else:
      ixm = (ix-1)%nmax
      ixp = (ix+1)%nmax
      ang = arr[ix,iy]-arr[ixm,iy]
      en += 0.5*(1.0 - 3.0*np.cos(ang)**2)
      ang = arr[ix,iy]-arr[ixp,iy]
      en += 0.5*(1.0 - 3.0*np.cos(ang)**2)
    return en
#=======================================================================
def all_energy(arr,nmax,rows,neighbours):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
      nmax (int) = side length of square lattice.
    Description:
      Function to compute the energy of the entire lattice. Output
      is in reduced units (U/epsilon).
	Returns:
	  enall (float) = reduced energy of lattice.
    """
    enall = 0.0
    for j in range(nmax):
        for i in range(rows):
            enall += one_energy(arr,i,j,nmax,neighbours,rows)
    return enall
#=======================================================================
def get_order(arr,nmax,rows):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
      nmax (int) = side length of square lattice.
    Description:
      Function to calculate the order parameter of a lattice
      using the Q tensor approach, as in equation (3) of the
      project notes.  Function returns S_lattice = max(eigenvalues(Q_ab)).
	Returns:
	  max(eigenvalues(Qab)) (float) = order parameter for lattice.
    """
    Qab = np.zeros((3,3))
    delta = np.eye(3,3)
    #
    # Generate a 3D unit vector for each cell (i,j) and
    # put it in a (3,i,j) array.
    #
    lab = np.vstack((np.cos(arr),np.sin(arr),np.zeros_like(arr))).reshape(3,rows,nmax)
    for a in range(3):
        for b in range(3):
            for i in range(rows):
                for j in range(nmax):
                    Qab[a,b] += 3*lab[a,i,j]*lab[b,i,j] - delta[a,b]
    return Qab
#=======================================================================
def MC_half_step(arr,Ts,nmax,rows,rng,neighbours):
    """
    Arguments:
	  arr (float(nmax,nmax)) = array that contains lattice data;
	  Ts (float) = reduced temperature (range 0 to 2);
      nmax (int) = side length of square lattice.
    Description:
      Function to perform one MC step, which consists of an average
      of 1 attempted change per lattice site.  Working with reduced
      temperature Ts = kT/epsilon.  Function returns the acceptance
      ratio for information.  This is the fraction of attempted changes
      that are successful.  Generally aim to keep this around 0.5 for
      efficient simulation.
	Returns:
	  accept/(nmax**2) (float) = acceptance ratio for current MCS.
    """
    #
    # Pre-compute some random numbers.  This is faster than
    # using lots of individual calls.  "scale" sets the width
    # of the distribution for the angle changes - increases
    # with temperature.
    scale=0.1+Ts
    accept = 0
    xran = rng.integers(0,high=rows, size=(rows,nmax))
    yran = rng.integers(0,high=nmax, size=(rows,nmax))
    aran = rng.normal(scale=scale, size=(rows,nmax))
    uran = rng.uniform(size=(rows,nmax))
    for i in range(rows):
        for j in range(nmax):
            ix = xran[i,j]
            iy = yran[i,j]
            ang = aran[i,j]
            en0 = one_energy(arr,ix,iy,nmax,neighbours,rows)
            arr[ix,iy] += ang
            en1 = one_energy(arr,ix,iy,nmax,neighbours,rows)
            if en1<=en0:
                accept += 1
            else:
            # Now apply the Monte Carlo test - compare
            # exp( -(E_new - E_old) / T* ) >= rand(0,1)
                boltz = np.exp( -(en1 - en0) / Ts )

                if boltz >= uran[i,j]:
                    accept += 1
                else:
                    arr[ix,iy] -= ang
    return accept

#=======================================================================
def main(program, nsteps, nmax, temp, pflag,save_file):
    """
    Arguments:
	  program (string) = the name of the program;
	  nsteps (int) = number of Monte Carlo steps (MCS) to perform;
      nmax (int) = side length of square lattice to simulate;
	  temp (float) = reduced temperature (range 0 to 2);
	  pflag (int) = a flag to control plotting.
    Description:
      This is the main function running the Lebwohl-Lasher simulation.
    Returns:
      NULL
    """
    comm = MPI.COMM_WORLD
    id = comm.Get_rank()
    ntasks = comm.Get_size()
    rng = np.random.default_rng(seed=id*int(time.time()))
    # Create and initialise lattice
    temp_rows = nmax // ntasks
    rem = nmax%ntasks
    if(id<rem)and(ntasks>1):
      rows=temp_rows + 1
    else:
      rows=temp_rows
    sub_lattice = np.zeros((rows,nmax),dtype=np.dtype('d'))
    up= (id-1)%ntasks
    down = (id+1)%ntasks
    if (id==0):
      lattice = np.zeros((nmax,nmax),dtype=np.dtype('d'))
      lattice = initdat(nmax)
      plotdat(lattice,pflag,nmax)
      if ntasks>1:
        counts = [0]*ntasks
        disp = [0]*ntasks
        for i in range(ntasks):
          rows_i = temp_rows
          if(i<rem):
            rows_i=temp_rows + 1
          counts[i] = rows_i*nmax
          disp[i] = sum(counts[:i])
        sendbuf = [lattice,counts,disp,MPI.DOUBLE]
    else:
      sendbuf = None
      lattice = None
    neighbours = np.zeros((2,nmax),dtype=np.dtype('d'))
    if ntasks>1:
      comm.Scatterv(sendbuf,sub_lattice[0:rows,:],root=0)
      comm.Sendrecv(sendbuf=sub_lattice[0,:],dest=up,sendtag=00,recvbuf=neighbours[1],source=down,recvtag=00)
      comm.Sendrecv(sendbuf=sub_lattice[rows-1,:],dest=down,sendtag=11,recvbuf=neighbours[0],source=up,recvtag=11)
    # Plot initial frame of lattice
    # Create arrays to store energy, acceptance ratio and order parameter
    energy = np.zeros(nsteps+1,dtype=np.dtype('d'))
    ratio = np.zeros(nsteps+1,dtype=np.dtype('d'))
    Qab = np.zeros((3,3))
    Qab_local = get_order(sub_lattice,nmax,rows)
    if ntasks>1:
      energy_local = np.zeros(nsteps+1,dtype=np.dtype('d'))
      ratio_local = np.zeros(nsteps+1,dtype=np.dtype('d'))
    
    # Set initial values in arrays
      energy_local[0] = all_energy(sub_lattice,nmax,rows,neighbours)
      comm.Reduce(Qab_local,Qab,op=MPI.SUM,root=0)
    else:
      energy[0] = all_energy(sub_lattice,nmax,rows,neighbours)
    if (id==0):
      order = np.zeros(nsteps+1,dtype=np.dtype('d'))
      Qab=Qab_local/(2*nmax*nmax)
      eigenvalues,eigenvectors = np.linalg.eig(Qab)
      order[0] = eigenvalues.max()
    comm.Barrier()
    # Begin doing and timing some MC steps.
    if (id==0):
      initial = MPI.Wtime()
    for it in range(1,nsteps+1):
        if ntasks>1:
          if(id%2==0):
            ratio_local[it] = MC_half_step(sub_lattice,temp,nmax,rows,rng,neighbours)
            comm.Send(sub_lattice[0,:],dest=up)
            comm.Send(sub_lattice[rows-1,:],dest=down)
          if(id%2==1):
            comm.Recv(neighbours[0],source=up)
            comm.Recv(neighbours[1],source=down)
            
          comm.Barrier()
          
          if(id%2==1):
            ratio_local[it] = MC_half_step(sub_lattice,temp,nmax,rows,rng,neighbours)
            comm.Send(sub_lattice[0,:],dest=up)
            comm.Send(sub_lattice[rows-1,:],dest=down)
          if(id%2==0):
            comm.Recv(neighbours[0],source=up)
            comm.Recv(neighbours[1],source=down)
            
          comm.Barrier()
          
          energy_local[it] = all_energy(sub_lattice,nmax,rows,neighbours)
          Qab_local = get_order(sub_lattice,nmax,rows)
          comm.Reduce(Qab_local,Qab,op=MPI.SUM,root=0)
          if (id==0):
            Qab=Qab/(2*nmax*nmax)
            eigenvalues,eigenvectors = np.linalg.eig(Qab)
            order[it] = eigenvalues.max()
        else:
          ratio[it] = MC_half_step(sub_lattice,temp,nmax,rows,rng,neighbours)
          energy[it] = all_energy(sub_lattice,nmax,rows,neighbours)
          Qab_local = get_order(sub_lattice,nmax,rows)
          Qab=Qab_local/(2*nmax*nmax)
          eigenvalues,eigenvectors = np.linalg.eig(Qab)
          order[it] = eigenvalues.max()
    comm.Barrier()
    if ntasks>1:
      comm.Reduce(energy_local,energy,op=MPI.SUM,root=0)
      comm.Reduce(ratio_local,ratio,op=MPI.SUM,root=0)
      recvbuf = [lattice,counts,disp,MPI.DOUBLE] if id==0 else 0
      comm.Gatherv(sub_lattice[0:rows,:],recvbuf,root=0)
    else:
      lattice = sub_lattice
    if (id==0):
      ratio = ratio/(nmax*nmax)
      ratio[0] = 0.5 # ideal value
      final = MPI.Wtime()
      runtime = final-initial
    
    # Final outputs
      print("{}: Size: {:d}, Steps: {:d}, T*: {:5.3f}: Order: {:5.3f}, Time: {:8.6f} s".format(program, nmax,nsteps,temp,order[nsteps-1],runtime))
    # Plot final frame of lattice and generate output file
      if save_file==0:
        savedat(lattice,nsteps,temp,runtime,ratio,energy,order,nmax)
      plotdat(lattice,pflag,nmax)
#=======================================================================
# Main part of program, getting command line arguments and calling
# main simulation function.
#
if __name__ == '__main__':
    if int(len(sys.argv)) == 6:
        PROGNAME = sys.argv[0]
        ITERATIONS = int(sys.argv[1])
        SIZE = int(sys.argv[2])
        TEMPERATURE = float(sys.argv[3])
        PLOTFLAG = int(sys.argv[4])
        SAVEFILE = int(sys.argv[5])
        main(PROGNAME, ITERATIONS, SIZE, TEMPERATURE, PLOTFLAG, SAVEFILE)
    else:
        print("Usage: python {} <ITERATIONS> <SIZE> <TEMPERATURE> <PLOTFLAG> <SAVEFILE>".format(sys.argv[0]))
#=======================================================================