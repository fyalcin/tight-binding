# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 12:35:01 2019

@author: fyalcin
"""

#
#tot = []
#for i in range(len(eigvecs)):
#    for j in range(i):
#        tot.append(np.dot(eigvecs[i].conj(),eigvecs[j]))
#print(np.round(tot,2))
#
#tot = [np.dot(eigvecs[i],eigvecs[i]) for i in range(len(eigvecs))]
#print(tot)



#plt.figure(figsize=(10,10))
#plt.plot(E, dos1)
#plt.savefig(figname, bbox_inches='tight')

#H,dosbef,dosaft,pdosbef,pdosaft,pdosdiff = plotdosandpdos(cell, gamma, delta, sigma, center, slices)
#plotcell(cell)

#graphene = geom.graphene().tile(6,0).tile(6,1)
##for i in range(len(graphene.xyz)):
##    graphene.xyz[i] = cell.atoms[i].coord
#xyz_remove = graphene.xyz[29]
#system = graphene.remove(29)
#H = Hamiltonian(system)
#r = (0.1,  1.44)
#t = (0. , -2.7)
#H.construct([r, t])
#N = len(system.xyz)
#h = np.zeros((N,N))
#for i in range(N):
#    for j in range(N):
#        h[i][j] = H[i,j]
##plt.scatter(system.xyz[:, 0], system.xyz[:, 1]);
##for ia in system: # loop over atomic indices
##    plt.annotate(str(ia), (system.xyz[ia, 0]+0.1, system.xyz[ia, 1]+0.1))
#es = H.eigenstate()
##es_fermi = es.sub(range(len(H) // 2 - 1, len(H) // 2 + 2))
##_, ax = plt.subplots(1, 3, figsize=(14, 2));
##for i, (n, c) in enumerate(zip(es_fermi.norm2(sum=False), 'rbg')):
##    ax[i].scatter(system.xyz[:, 0], system.xyz[:, 1], 600 * n, facecolor=c, alpha=0.5);
##    ax[i].scatter(xyz_remove[0], xyz_remove[1], c='k', marker='*'); # mark the removed atom
#    
##E = np.linspace(-4, 4, 400)
##plt.plot(E, es.DOS(E));
##plt.xlabel(r'$E - E_F$ [eV]'); plt.ylabel(r'DOS at $\Gamma$ [1/eV]');
#
#E = np.linspace(-1, -.5, 100)
#dE = E[1] - E[0]
#PDOS = es.PDOS(E).sum(1) * dE # perform integration
#plt.scatter(system.xyz[:, 0], system.xyz[:, 1], 500 * PDOS);
#plt.scatter(xyz_remove[0], xyz_remove[1], c='k', marker='*'); # mark the removed atom

#plotcell(cell)
#cell = SimulationCell(1,1,Lattice="PG")
#H, S, ads, eigs = solveSecularSO(cell, gamma, alpha = 0)
#for i in range(len(ads[0])):
#    plt.plot(range(len(ads[:,i])),ads[:,i])
#plotcell(cell)
#S = np.array([[1, 7.31],[7.31, 1]])
#S = np.identity(2)
#pdos1 = PDOS(band[0], -1, -0.5, eigs, np.identity(len(band[0])))
#a = np.random.rand(10000, 10000)
#b = np.random.rand(10000, 10000)
#np.dot(a, b)
#H, S, ads, eigs = solveSecularSO(cell, gamma, alpha = 5)

#pdos2 = PDOS(ads[0], -100, 500, eigs)
#
#
#plotPDOS(cell, pdos1)
#
#kpoints = kpath(4000)
#
#def helper(kpoint):
#    return solveSecularAO(cell, [kpoint])
#
#def pool_factorizer_map(nums, nprocs):
#    bands = []
#    with ProcessPoolExecutor(max_workers=nprocs) as executor:
#        return executor.map(helper, nums)
#        
#def mp_factorizer_map(nums, nprocs):
#    with mp.Pool(nprocs) as pool:
#        return pool.map(helper, nums)
#
#start = time.time()
##a = pool_factorizer_map(kpoints, 4)
##a = mp_factorizer_map(kpoints, 4)
#H, S, ads, eigs = solveSecularAO(cell, kpoints, alpha=0)
#end = time.time()
#
#print(end-start)

#H, ads, eigs = solveSecularSO(cell1, kpoints, alpha=0)
#pdos = PDOS(ads[0], 3, 5, eigs)
#plotPDOS(cell1, pdos)
#cell1.populateWlatticepts(BuildJunction(7,18))
#plotcell(cell1)
#end = time.time()
#print(end-start)
#for _ in range(10):
#    start = time.time()
#    H,ads,eigs = solveSecularAO(cell, kpoints, alpha = 0)
#    end = time.time()
#    print(end-start)
#    tottime += end-start
#print(tottime)
#H,ads,eigs = solveSecularAO(cell, kpath(3), alpha = 0)
#end = time.time()
#cell = SimulationCell(1,1, Lattice = "PG")
#H, ads, eigs = solveSecularAO(cell, kpath(50))

               
#plt.figure(figsize=(10,10))
#for i in range(len(ads[0])):
#    plt.plot(range(len(ads[:,i])),ads[:,i])
#plt.savefig(figname, bbox_inches='tight')
    
#figname = "dos.png"
#xarray, yarray = gaussianize(ads[0], np.linspace(-3,3,1000))
#plt.figure(figsize=(10,10))
#plt.plot(xarray, yarray)
#plt.savefig(figname, bbox_inches='tight')

#end = time.time()
#
#print(end-start)

#cell = MakeJunction(2)
#plotcell(cell, neighbors = True)
#Htot,ads = solveSecular(cell, kpoint, alpha=0.1)
#cell = SimulationCell(6,6,Lattice="flake",q=1)

#plt.figure(figsize=(15,15))
#plt.hist(ads[0],len(ads[0]) )

#for alpha in np.linspace(0,0.1,11):
#    plt.close()
#    Htot, ads = solveSecular(cell, kpoint, alpha)
#    figname = "dos" + str(alpha) + ".png"
#    xarray, yarray = gaussianize(ads[0], np.linspace(-3,3,1000))
#    plt.figure(figsize=(10,10))
#    plt.plot(xarray, yarray)
#    plt.savefig(figname, bbox_inches='tight')

#xcoord = [atom.coord[0] for atom in cell.atoms]
#ycoord = [atom.coord[1] for atom in cell.atoms]
#plt.figure(figsize=(30,30))
#ax = plt.axes()
#plt.scatter(xcoord, ycoord, s=30, c="r")
#
#for item in cell.gridxy:
#    plt.annotate("["+str(item[1][0])+","+str(item[1][1])+"]", (item[0][0]-0.1,item[0][1]-0.1), size = 6)
#    Hex = RegularPolygon((item[0][0],item[0][1]), numVertices=6, radius=latConst,orientation=np.radians(30), facecolor=(1,0,0,0.1*item[2]), edgecolor = (0,0,0,0.1))
#    ax.add_patch(Hex)

#nxcoord = [atom.coord[0] for atom in cell.neighborCell]
#nycoord = [atom.coord[1] for atom in cell.neighborCell]
#plt.scatter(nxcoord, nycoord, s=2, c="b")

#for atom in cell.atoms:
##    plt.annotate(str(atom.UCindex[0]), (atom.coord[0]-0.2, atom.coord[1]+0.2), size=6)
##    plt.annotate(str(atom.UCindex[1])+")", (atom.coord[0], atom.coord[1]+0.2), size=6)
##    plt.annotate(str(atom.MUCindex[0]), (atom.coord[0]-0.4, atom.coord[1]-0.2), size=6)
#    plt.annotate(str(len(atom.neighbors)), (atom.coord[0]-0.4, atom.coord[1]+0.2), size=6)
#    plt.annotate("["+str(atom.Spaceindex[0])+","+str(atom.Spaceindex[1])+"]", (atom.coord[0]-0.3, atom.coord[1]+0.2), size=6)
#    plt.annotate(str(atom.MUCindex[1]), (atom.coord[0]-0.4, atom.coord[1]-0.2), size=6)
#plt.savefig('junction2.pdf', format='pdf')
#figname = "indices"+cell.Lattice+str(cell.N1)+str(cell.N2)+"-"+str(cell.q)+".png"
#plt.savefig(figname,bbox_inches='tight')
#for atom in cell.neighborCell:
#    plt.annotate(str(atom.UCindex[0]), (atom.coord[0]-0.2, atom.coord[1]+0.2), size=6)
#    plt.annotate(str(atom.Spaceindex), (atom.coord[0]-0.2, atom.coord[1]+0.2), size=6)
#    plt.annotate(str(len(atom.neighbors)), (atom.coord[0]-0.4, atom.coord[1]+0.2), size=6)
#    plt.annotate(str(atom.MUCindex[1]), (atom.coord[0]-0.4, atom.coord[1]-0.2), size=6)
#plt.title('Hexagonal Lattice')
#plt.grid(True)
#plt.axes().set_aspect(aspect='auto')
#plt.axes().set_aspect('equal', 'datalim')
#
#plt.ylim(0,50)
#plt.xlim(0,50)
#plt.show()

#plt.subplot(1,2,2)
#Htot, ads, cell, eigfuncs = butterfly('PG',1,1,251,orbitals="s")
#center = 0
#delta = 0.1
#sigma = 0.04
#slices = int(2*delta/sigma)
#cell = SimulationCell(6, 6, Lattice="PG", periodic = True)
#E = np.linspace(-10,10,1000)
#dos1 = DOS(E,eigvals[0],E[2]-E[0])
#plt.plot(E,dos1)
#plt.title("before")
#dosbef = dos1
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, orthogonal = True)
#cell.populateWlatticepts(BuildJunction(5,12))
#print(DataFrame(np.round(H,2)))
#fig2 = plt.figure(2)
#plotcell(cell)

#cell.remove(31)
#fig1 = plt.figure(1)
#plt.figure(figsize=(12,9))
#plt.subplot(1,2,1)