from __future__ import division
# Overlapping Branch and Bound, data parallel algorithm with one-at-a-time hashing
def runpar(f, g, H, Lg, Lh, l, u, bound, circle, A=None, b=None, E=None, d=None, Tol=1e-2, Heur=0, TolType='r', Vis=0, nsize=12, SD=1, qpsolver='cvxopt'):

    # Optional Inputs
    # Tolerance
    # Heuristic lattice (0 - off, 1 - on)
    # Tolerance type (r - relative, a - absolute)
    # Visualisation (0 - off, 1 - on)
    # Step Debugging (0 - off, 1 - on)
    # Node Size
    # QP Solver (cvxopt, quadprog, nag)

    # MPI
    from mpi4py import MPI

    # MPI comm
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Number of processors
    numprocs = comm.Get_size()

    # Print node names for each processors
    print('Processor %i is on node %s') % (rank,str(MPI.Get_processor_name()))

    # Create local groups for load balacing on each node
    from numpy import ceil
    nodes = [[0]]
    bins = int(ceil(numprocs/nsize))
    for i in range(bins):
        nodes.append([])
    for p in range(1,numprocs):
        nodes[(p % bins)+1].append(p)
    print nodes

    # Get D
    D = len(l)
    if(D > 25):
        raise RuntimeError('Hashing for greater than 25D not yet implemented.')

    # Balance Test
    def balance_test(lst):

        for p1 in reversed(range(0,len(lst))):

            for p2 in range(0,p1):

                if( abs(lst[p1]-lst[p2])/(max(min(lst[p1],lst[p2]),1))  > 0.1):

                    return True

        return False

    if(Heur == 0):

        ksp = 3**D
        kspi = ksp

    else:
        # Load relevant normalised lattice
        from numpy import loadtxt, array, sqrt
        from pkg_resources import resource_stream
        if(D == 2):
            lat = array([[0., 0.], [1., 0.], [-1., 0.],
            [0.5, sqrt(3.)/2.], [-0.5, sqrt(3.)/2.],
            [0.5, -sqrt(3.)/2.], [-0.5, -sqrt(3.)/2.]])
        elif(D == 3):
            lat = loadtxt(resource_stream('obb','lattices/d3'))
        elif(D == 4):
            lat = loadtxt(resource_stream('obb','lattices/d4'))
        elif(D == 5):
            lat = loadtxt(resource_stream('obb','lattices/d5'))
        elif(D == 6):
            lat = loadtxt(resource_stream('obb','lattices/e6'))
        elif(D == 7):
            lat = loadtxt(resource_stream('obb','lattices/e7'))
        elif(D == 8):
            lat = loadtxt(resource_stream('obb','lattices/e8'))
        else:
            raise RuntimeError('A lattice for '+str(D)+' Dimensions has yet to be provided.')

        # Get kissing number + 1
        ksp = lat.shape[0]
        kspi = 2**D

    # Master Process
    if(rank == 0):

        # Necessary functions
        from itertools import product
        from numpy import array, sqrt, empty, zeros, inf, hstack, isinf, argmin, argmax
        from numpy.linalg import norm
        from math import floor
        from time import time

        # Timer
        timer = time()

        # Set up initial circle
        r = norm(u-l)/2 # Radius
        xc = (u+l)/2 # Centre

        if(Heur == 0):

            # Split circle into more circles
            inc = r/sqrt(D) # Increment
            rn = r/2 # New radius

            # Create square spoke configuration
            spk = array([p for p in product([-inc,0,inc], repeat=D)])

        else:

            # Split circle into more circles
            inc = r # Increment
            rn = r/2 # New radius

            # Scale configuration
            spk = inc*lat

        # Array to distribute amongst processes
        darray = empty((kspi,D+1))
        reqidsend = []

        # Create surrounding circles
        for i in range(0,kspi):

            # Create centre
            xcn = xc + spk[i,:]

            # Add to distribution list
            darray[i,:] = hstack([array([rn]),xcn])

        # Distribute data evenly amongst processes
        ihi = 0
        trem = kspi
        prem = numprocs-1

        for p in range(0,numprocs-1):

            tproc = int(round(trem/prem))

            ilo = ihi + 1
            ihi = ihi + tproc

            print('Master sending data elements [%i:%i] to processor %i') % (ilo-1,ihi,p+1)
            comm.Ssend(array([ihi-ilo+1],dtype='i'), dest=p+1, tag=9) # send array size
            reqidsend.append(comm.Issend(darray[ilo-1:ihi,:], dest=p+1, tag=10)) # send data array

            prem = prem - 1
            trem = trem - tproc


        # Initial communication setup
        lclists = [1 for i in range(1,numprocs)]
        zerolist = [0 for i in range(1,numprocs)]
        confd = empty(1,dtype='i')
        reqbrecv = []
        reqbsend = []
        rarr = empty((numprocs-1,2))
        reqdreq = [MPI.REQUEST_NULL for i in range(1,numprocs)]
        reqconf = [MPI.REQUEST_NULL for i in range(1,numprocs)]
        reqdreqn = [MPI.REQUEST_NULL for i in range(1,numprocs)]
        reqconfn = [MPI.REQUEST_NULL for i in range(1,numprocs)]

        # Initial local upper bound
        U = inf

        # Initial rlist, xclist and comm setup
        rlist = []
        xclist = []
        rxcndat = empty(3,dtype='i')
        dwait = 0
        reqxcnsend = MPI.REQUEST_NULL

        # Loop counter
        #itr = 0

        # Asynchronously receive U and len(clist)
        for p in range(0,numprocs-1):

            # New Receive
            reqbrecv.append(comm.Irecv(rarr[p], source=p+1, tag=2))

            # Send U
            reqbsend.append(comm.Issend(array([U]), dest=p+1, tag=3))

        # Make sure data has been sent
        print('Master waiting for initial data sends to complete...')
        for p in range(0,numprocs-1):
            reqidsend[p].Wait()
        print('done.')

        # Ansynchronously receive xcn
        reqxcndat = comm.Irecv(rxcndat, source=MPI.ANY_SOURCE, tag=12)

        #and(rad > 1e-15)
        #and((time()-timer) < 3000)
        while((lclists != zerolist)):

            # Update iteration count
            #itr = itr + 1

                 # Asynchronously receive U and len(clist)
            for p in range(0,numprocs-1):

                if(reqbrecv[p].Test()):

                    # Update U and len(clist)
                    if(rarr[p,0] < U):
                        U = rarr[p,0]
                    lclists[p] = rarr[p,1]

                    # New Receive
                    reqbrecv[p] = comm.Irecv(rarr[p], source=p+1, tag=2)

            # Asynchronously send U
            for p in range(0,numprocs-1):

                # If previous communication has finished start new
                if(reqbsend[p].Test()):

                    # Send U
                    reqbsend[p] = comm.Issend(array([U]), dest=p+1, tag=3)


            # If data check requested, process accordingly
            # If data size received and previous send has completed
            if((reqxcndat.Test())and(reqxcnsend.Test())):

                # If not waiting to receive data
                if(dwait == 0):

                    #print('Data check from processor %i received.') % rxcndat[0]

                    # Receive data from processor
                    rxcn = empty(rxcndat[1],dtype='i')
                    reqxcnrecv = comm.Irecv(rxcn, source=rxcndat[0], tag=14)

                    # Waiting for data
                    dwait = 1

                # If data received
                if(reqxcnrecv.Test()):

                    # Array to send
                    nc = zeros(rxcndat[1],dtype='i')

                    # If radius already in list
                    if rxcndat[2] in rlist:

                        # Check if xcns already exist
                        ridx = rlist.index(rxcndat[2])
                        for i in range(0,rxcndat[1]):
                            for arr in xclist[ridx]:
                                if(rxcn[i] == arr):
                                    nc[i] = 1
                            if(nc[i] == 0):
                                xclist[ridx].append(rxcn[i])

                    # Else add radius and xcns to respective lists
                    else:
                        rlist.append(rxcndat[2])
                        xclist.append(list(rxcn))

                    # Send result of check back
                    reqxcnsend = comm.Issend(nc, dest=rxcndat[0], tag=14)
                    #print('Data check returned to processor %i.') % rxcndat[0]

                    # Ansynchronously receive xcn
                    reqxcndat = comm.Irecv(rxcndat, source=MPI.ANY_SOURCE, tag=12)

                    # Not waiting for data
                    dwait = 0

            # Load Balancing

            # Load Balancing spread across nodes

            # Take snapshot of load
            slist = lclists[:]

            # Get load on each node
            lnode = []
            for node in nodes:

                # Master proc, skip
                if(len(node)==1):

                    pass

                # Actual Node with workers
                else:

                    # Create list
                    blist = []
                    for n in node:
                        blist.append(slist[n-1])

                    # Get load per node
                    lnode.append(blist)

            # Get total load per node
            tload = [sum(lnode[i]) for i in range(len(lnode))]

            # Iteration counter
            it = 0

            while((balance_test(tload))and(it < numprocs)):

                # Subtract smallest load
                ml = min(tload)
                tload = map(lambda x: x-ml, tload)

                # Find amount to redistribute
                load = int(floor(sum(tload)/len(lnode)))

                # Find nodes with smallest & largest load
                minn = argmin(tload)
                maxn = argmax(tload)

                # Calculate node load to send
                nload = min(tload[maxn],max(load-tload[minn],0))

                # Find procs with smallest & largest load on max/min nodes
                lminp = argmin(lnode[minn])
                lmaxp = argmax(lnode[maxn])

                # Get proc number
                gminp = nodes[minn+1][lminp]
                gmaxp = nodes[maxn+1][lmaxp]

                # Calculate actual load to send
                lfrac = 3
                aload = min(int(lnode[maxn][lmaxp]/lfrac),max(nload-lnode[minn][lminp],0))

                # Fraction of processors on node
                dfrac = len(nodes[minn+1])

                # Send load from largest to smallest if previous send has completed and aload/dfrac > 0
                if((reqdreqn[gmaxp-1].Test())and(reqconfn[gmaxp-1].Test())and(int(aload/dfrac) > 0)):

                    # Send load from largest to smallest
                    print('Intercomm: Master allocating %i subproblems from processor %i to node %i processors ') % (aload,gmaxp,minn+1),
                    print nodes[minn+1][0:]

                    # Send data allocation request
                    reqdreqn[gmaxp-1] = comm.Issend(array([minn+1,aload],dtype='i'), dest=gmaxp, tag=24)
                    print('Intercomm: Data allocation request sent asynchronously.')

                    # Wait for confirmation of data send
                    reqconfn[gmaxp-1] = comm.Irecv(confd, source=gmaxp, tag=26)

                # Update load per node
                for i in range(len(nodes[minn+1])):
                    lnode[minn][i]+= min(int(aload/dfrac),lnode[maxn][lmaxp])
                lnode[maxn][lmaxp] -= min(aload,lnode[maxn][lmaxp])

                # Update tload
                tload = [sum(lnode[i]) for i in range(len(lnode))]

                # Iteration counter
                it += 1

            # Load Balancing on each node

            # For each node (excl master)
            for node in nodes:

                # Master proc, skip
                if(len(node)==1):

                    pass

                # Actual Node with workers
                else:

                    # Copy list
                    blist = []
                    for n in node:
                        blist.append(lclists[n-1])

                    # Subtract smallest load
                    ml = min(blist)
                    blist = map(lambda x: x-ml, blist)

                    # Find load amount to redistribute
                    load = int(floor(sum(blist)/len(blist)))

                    # Iteration counter
                    biter = 0

                    while((balance_test(blist))and(biter < len(blist))):

                        # Find procs with smallest & largest load
                        minp = argmin(blist)
                        maxp = argmax(blist)

                        # Convert to global proc numbers
                        gminp = node[minp]
                        gmaxp = node[maxp]

                        # Calculate actual load to send
                        aload = max(load-int(blist[minp]),0)

                        # If previous send has completed and aload > 0
                        if((reqdreq[gmaxp-1].Test())and(reqconf[gmaxp-1].Test())and(aload > 0)):

                            # Send load from largest to smallest
                            print('Master allocating %i subproblems from processor %i to processor %i') % (aload,gmaxp,gminp)

                            # Send data allocation request
                            reqdreq[gmaxp-1] = comm.Issend(array([gminp,aload],dtype='i'), dest=gmaxp, tag=0)
                            print('Data allocation request sent asynchronously.')

                            # Wait for confirmation of data send
                            reqconf[gmaxp-1] = comm.Irecv(confd, source=gmaxp, tag=22)

                        # Update blist
                        blist[minp] += min(aload,blist[maxp])
                        blist[maxp] -= min(aload,blist[maxp])

                        # Update iteration counter
                        biter += 1


        # If all finished, send kill signals
        print('Master sending kill signals...')
        for p in range(0,numprocs-1):
            comm.Ssend(array([1],dtype='i'), dest=p+1, tag=4)

        # Clean up redundant messages
        for p in range(0,numprocs-1):
            #reqbsend[p].Cancel()
            reqbrecv[p].Cancel()
            #reqdreq[p].Cancel()
            #reqbsend[p].Free()
            reqbrecv[p].Free()
            #reqdreq[p].Free()

        # Receive xopt from all processors
        reqe = []
        earr = empty((numprocs-1,D))
        U = inf
        for p in range(0,numprocs-1):

            # New Receive
            #print('Master posting asynchronous receive of final data from processor %i') % (p+1)
            reqe.append(comm.Irecv(earr[p], source=p+1, tag=5))
            #print('posted.')

        for p in range(0,numprocs-1):

            print('Master waiting to receive final data from processor %i') % (p+1)
            reqe[p].Wait()
            print('Final data from processor %i received') % (p+1)

            # Find xopt with smallest U
            if(not(any(isinf(earr[p])))):
                if(f(earr[p]) < U):
                    xopt = earr[p]
                    U = f(xopt)
        #xopt = empty(D)
        # How do you count iterations?
        itr = None

        # Output end result
        print('Minimum value of %f at (') % U, #f(xopt),
        for i in range(0,D-1):
            print('%f,') % xopt[i],
        print('%f)') % xopt[D-1]
        print('found with'),
        if(TolType == 'a'):
            print('absolute'),
        else:
            print('relative'),
        print('tolerance %f.') % Tol
        tol = Tol
        xs = xopt
        fxs = U #f(xopt)
        print('Elapsed time %f seconds.') % (time()-timer)

        return xs, fxs, tol, itr

    # Worker processes
    else:

        # import necessary functions
        from heapq import heappush, heappop
        from numpy import array, pi, sqrt, dot, identity, hstack, vstack, inf, empty, ones, zeros, tile
        from numpy.linalg import norm
        from itertools import product
        from sys import float_info
        from math import floor
        from time import time

        # Vector Hash function: http://stackoverflow.com/questions/5928725/hashing-2d-3d-and-nd-vectors
        #htsize = 499 # Hash Table size
        res = 1e-5 # Resolution
        # Twenty five arbitrarily chosen 8-digit primes
        primes = [73856093, 19349663, 83492791, 15485863, 86028121, 32452843, 67867967, 49979687,
                  86028157, 13769629, 95189161, 29223841, 83972821, 39271703, 79186817, 47286739,
                  98767549, 14139199, 64244023, 18551173, 55621541, 10920859, 53615137, 24405817,
                  43176253 ]
        # Hash function
        def hashnd(x):

            result = 0
            for i in range(0,D):
                result = result ^ (int(floor(x[i]/res))*primes[i])
            return result

        # Hash radius
        def hashrn(r):

            return int(r*1e8)

        # Initialise empty arrays
        if(A is None):
            A = empty((0,D))
            b = empty(0)

        if(E is None):
            E = empty((0,D))
            d = empty(0)

        # QuadProg++ solver (fast!)
        if(qpsolver == 'quadprog'):

            # Import QuadProg++ solver
            from PyQuadProg import PyQuadProg

            # Check if circle has feasible point
            def mfeasible(c):

                # Solve QP to check feasibility
                sol = PyQuadProg(2*identity(D),-2*c.xc,E.transpose(),-1*d,vstack([identity(D),-1*identity(D),-1*A]).transpose(),hstack([-l,u,b]))
                mxopt = sol.x.getArray()
                mr = dot(mxopt,mxopt) + dot(mxopt,-2*c.xc) + dot(c.xc,c.xc)

                # Check if point lies inside domain
                if(mr < (c.r**2)):
                    mf = 1
                else:
                    mf = 0
                    mxopt = None

                return mf, mxopt

        # CVXOPT QP solver (slow)
        elif(qpsolver == 'cvxopt'):

            # Import cvxopt solver
            from cvxopt import matrix
            from cvxopt.solvers import qp, options

            # Set tolerance options
            options['show_progress'] = False
            options['abstol'] = 1e-9
            options['reltol'] = 1e-8
            options['feastol'] = 1e-9

            # Check if circle has feasible point
            def mfeasible(c):

                # Solve QP to check feasibility
                sol = qp(matrix(2*identity(D)),matrix(-2*c.xc),matrix(vstack([-1*identity(D),identity(D),A])),matrix(hstack([-l,u,b])),matrix(E),matrix(d))
                mxopt = array(sol['x']).flatten()
                mr = sol['primal objective']
                mr = mr + dot(c.xc,c.xc)

                # Check if point lies inside domain
                if(mr < (c.r**2)):
                    mf = 1
                else:
                    mf = 0
                    mxopt = None

                return mf, mxopt

        # NAG QP solver (fast!)
        elif(qpsolver == 'nag'):

            # Import QP Solver
            from qpsolver_lincon import qpsolver_lincon

            # Check if circle has feasible point
            def mfeasible(c):

                # Solve QP to check feasibility
                mxopt = c.xc.copy()
                mr = qpsolver_lincon(2*identity(D),-2*c.xc,hstack([l,-inf,d]),hstack([u,b,d]),mxopt,vstack([A,E]),D,A.shape[0]+E.shape[0])
                mr = mr + dot(c.xc,c.xc)

                # Check if point lies inside domain
                if(mr < (c.r**2)):
                    mf = 1
                else:
                    mf = 0
                    mxopt = None

                return mf, mxopt

        # Visualisation
        if(Vis == 1):

            from matplotlib.pyplot import figure, Rectangle, Circle, gca, show, title, axis, draw
            from colorsys import hsv_to_rgb
            from random import random

            # Draw bound constraints [l,u]
            fig = figure('Processor '+str(rank))
            gca().add_patch(Rectangle((l[0],l[1]), u[0]-l[0], u[1]-l[1], fill=False))
            axis([l[0]-1,u[0]+1,l[1]-1,u[1]+1])
            title('Processor '+str(rank))
            show(block=False)

            # Circle drawing procedure
            def drawc(c):

                # Draw circle colured according to partition
                gca().add_patch(Circle((c.xc[0],c.xc[1]), radius=c.r , color=hsv_to_rgb(random(),1,1), fill=False))
                axis('equal')
                draw()

        # Priority queue
        clist = []                      # list of entries arranged in a heap

        # Add a new task
        def add_task(task):
            heappush(clist, (task.lbound, task))

        # Remove and return the lowest priority task. Raise KeyError if empty.
        def pop_task():
            while clist:
                return heappop(clist)
            raise KeyError('Pop from an empty priority queue!')

        # Remove and return (in an mpi transfer array) the n lowest priority tasks
        def npop_task(n):
            carray = empty((n,D+2))
            for i in range(0,n):
                celm = pop_task()
                carray[i,:] = hstack([array([celm[0],celm[1].r]),celm[1].xc])
            return carray

        # Pick up scattered list
        rsize = array([0],dtype='i')
        comm.Recv(rsize, source=0, tag=9)
        rarray = empty((rsize[0],D+1))
        comm.Recv(rarray, source=0, tag=10)
        print('Processor %i received initial data') % rank

        # Initial local upper bound
        U = inf

        # Bound communication setup
        rarr = empty(1)
        karr = empty(1,dtype='i')
        kill = 0
        dwait = 0
        dwaits = 0
        reqdsend = MPI.REQUEST_NULL
        reqdsizesend = MPI.REQUEST_NULL
        reqconf = MPI.REQUEST_NULL
        dwaitn = 0
        dwaitsn = 0
        reqdsendn = [MPI.REQUEST_NULL for i in range(1,numprocs)]
        reqdsizesendn = [MPI.REQUEST_NULL for i in range(1,numprocs)]
        reqconfn = MPI.REQUEST_NULL

        # Bound intial list and convert to priority queue
        for k in range(0,rsize[0]):

            # Create circle
            cn = circle(rarray[k,1:D+1],rarray[k,0])

            # Check if circle is feasible
            mfeas, cxopt = mfeasible(cn)

            # If circle has a feasible point (i.e. it's feasible)
            if(mfeas != 0):

                # Bound circle
                cn.lbound = bound(cn,Lg,Lh,f,g,H,D)

                # Upper bound
                ubound = f(cxopt)

                # Update U and xopt
                if(ubound < U):
                    U = ubound
                    xopt = cxopt

                # Add to priority queue
                add_task(cn)

                # Visualise
                if(Vis ==1):
                    drawc(cn)

        # Asynchronously send U and len(clist)
        reqbsend = comm.Issend(array([U,len(clist)]), dest=0, tag=2)

        # Asynchronously receive U
        reqbrecv = comm.Irecv(rarr, source=0, tag=3)

        # Debug Output
        if (SD == 1):
            #if(rank == 0):
            print('Time: %f') % (time()),
            print('Processor %i Number of elements: %i') % (rank,len(clist))
            print('U: %f') % U
            #print('L: %f') % L
            #print('Circle radius: %e') % rad
            print('----------------------------')

        # Intial data requests
        rq_proc_load = empty(2,dtype='i')
        rq_proc_loadn = empty(2,dtype='i')
        reqd = comm.Irecv(rq_proc_load, source=0, tag=0)
        reqdn = comm.Irecv(rq_proc_loadn, source=0, tag=24)

        # Initial kill signal request
        reqk = comm.Irecv(karr, source=0, tag=4)

        while(kill == 0):

            # Set tolerance
            if((TolType == 'r')and(abs(U) > float_info.epsilon)):
                cutoff = U - Tol*abs(U)
            else:
                cutoff = U - Tol

            # Prune list
            i = 0
            while(i < len(clist)):
                if(clist[i][0] > cutoff):
                    del clist[i]
                    i=i-1
                i = i+1

            # If there is work left to do
            if(len(clist) != 0):

                # Get circle with smallest lbound
                L, cslb = pop_task()

                if(Heur == 0):

                    # Split circle into more circles
                    inc = (cslb.r)/sqrt(D) # Increment
                    rn = (cslb.r)/2 # New radius
                    xc = cslb.xc # Centre

                    # Create square spoke configuration
                    spk = array([p for p in product([-inc,0,inc], repeat=D)])

                else:

                    # Split circle into more circles
                    inc = cslb.r # Increment
                    rn = (cslb.r)/2 # New radius
                    xc = cslb.xc # Centre

                    # Scale configuration
                    spk = inc*lat

                # Create surrounding circles
                # Create centre
                nc = zeros(ksp)
                xcn = tile(xc,(ksp,1)) + spk

                # Check if circles exists locally
                for k in range(0,len(clist)):
                    for i in range(0,ksp):
                        if(all(xcn[i,:] == clist[k][1].xc)):
                            nc[i] = 1

                # Check if circles exists globally
                # Remove local duplicates
                dxcn = []
                hdxcn = []
                for i in range(0,ksp):
                    if(nc[i] == 0):
                        dxcn.append(xcn[i,:])
                        hdxcn.append(hashnd(xcn[i,:]))
                xcn = array(dxcn)
                hxcn = array(hdxcn,dtype='i')
                del dxcn[:]
                del hdxcn[:]

                # Send data asynchronously
                #print('Processor %i sending data check to master processor.') % rank
                comm.Send(array([rank,len(hxcn),hashrn(rn)],dtype='i'), dest=0, tag=12) # send array size
                _ = comm.Issend(hxcn, dest=0, tag=14) # send data array
                #print('Data check sent. Awaiting response.')

                # Receive data asynchronously
                nc = empty(len(hxcn),dtype='i')
                reqdcheck = comm.Irecv(nc, source=0, tag=14)
                #print('Data check completed successfully.')

                # Wait for data check
                print('Processor %i started waiting for data check.') % rank
                wtimer = time() # Time wait time
                while(not((reqdcheck.Test())or(reqk.Test()))):
                    pass

                # Break loop if kill signal received
                if(reqk.Test()):
                    break

                print('Processor %i spent %f seconds waiting for data check from master.') % (rank,time()-wtimer)

                # Create circles which don't already exist
                for i in range(0,len(xcn)):

                    # If circle doesn't already exist somewhere
                    if(nc[i] == 0):

                        # Create circle
                        cn = circle(xcn[i,:],rn)
                        mfeas, cxopt = mfeasible(cn)

                        # If circle has a feasible point (i.e. it's feasible)
                        if(mfeas != 0):

                            # Lower bound
                            cn.lbound = bound(cn,Lg,Lh,f,g,H,D)

                            # Upper bound
                            ubound = f(cxopt)

                            # Update U and xopt
                            if(ubound < U):
                                U = ubound
                                xopt = cxopt

                            # Add to priority queue
                            add_task(cn)

                            # Visualise
                            if(Vis ==1):
                                drawc(cn)

                # Intercomm: Send data if requested and previous data send has completed
                if((reqdn.Test())and(MPI.Request.Testall(reqdsendn))and(MPI.Request.Testall(reqdsizesendn))and(reqconfn.Test())):

                    print('Intercomm: Node %i needs data from processor %i') % (rq_proc_loadn[0],rank)
                    # Workload to send
                    nit = min(int(rq_proc_loadn[1]/len(nodes[rq_proc_loadn[0]])),len(clist)/len(nodes[rq_proc_loadn[0]]))

                    if(nit > 0):

                        print('Intercomm: Sending %i subproblems from processor %i to node %i each of processors') % (nit,rank,rq_proc_loadn[0]),
                        print nodes[rq_proc_loadn[0]][0:]

                        # Send data asynchronously
                        for p in nodes[rq_proc_loadn[0]]:
                            reqdsizesendn[p-1] = comm.Issend(array([nit,rank],dtype='i'), dest=p, tag=30) # send array size
                            reqdsendn[p-1] = comm.Issend(npop_task(nit), dest=p, tag=28) # send data array

                        print('Intercomm: Data sent asynchronously.')
                        # Send confirmation to master
                        reqconfn = comm.Issend(array([1],dtype='i'), dest=0, tag=26)

                    else:
                        print('Intercomm: Not enough data to send.')
                        # Send confirmation to master
                        reqconfn = comm.Issend(array([1],dtype='i'), dest=0, tag=26)

                    # Wait for other data requests
                    reqdn = comm.Irecv(rq_proc_loadn, source=0, tag=24)

                # Send data if requested and previous data send has completed
                if((reqd.Test())and(reqdsend.Test())and(reqdsizesend.Test())and(reqconf.Test())):

                    print('Processor %i needs data from processor %i') % (rq_proc_load[0],rank)
                    # Workload to send
                    nit = min(rq_proc_load[1],len(clist))

                    print('Sending %i subproblems from processor %i to processor %i') % (nit,rank,rq_proc_load[0])
                    # Send data asynchronously
                    reqdsizesend = comm.Issend(array([nit,rank],dtype='i'), dest=rq_proc_load[0], tag=11) # send array size
                    reqdsend = comm.Issend(npop_task(nit), dest=rq_proc_load[0], tag=1) # send data array

                    print('Data sent asynchronously.')
                    # Send confirmation to master
                    reqconf = comm.Issend(array([1],dtype='i'), dest=0, tag=22)

                    # Wait for other data requests
                    reqd = comm.Irecv(rq_proc_load, source=0, tag=0)

                # If not waiting to receive data size
                if(dwaits == 0):

                    # Receive data
                    print('Processor %i asynchronously waiting to receive data.') % rank
                    darrsize = array([0,0],dtype='i')
                    reqdsize = comm.Irecv(darrsize, source=MPI.ANY_SOURCE, tag=11)

                    # Waiting for data size
                    dwaits = 1

                # If data size received
                if(reqdsize.Test()):

                    # If not waiting to receive data
                    if(dwait == 0):

                        darr = empty((darrsize[0],D+2))
                        reqds = comm.Irecv(darr,source=darrsize[1], tag=1)

                        # Waiting for data
                        dwait = 1

                    # If data received
                    if(reqds.Test()):

                        for i in range(0,darrsize[0]):
                            ctemp = circle(darr[i,2:2+D],darr[i,1])
                            ctemp.lbound = darr[i,0]
                            add_task(ctemp)

                        print('Data received succesfully.')

                        # Not waiting for data size or data
                        dwaits = 0
                        dwait = 0

                # Intercomm: If not waiting to receive data size
                if(dwaitsn == 0):

                    # Receive data
                    print('Intercomm: Processor %i asynchronously waiting to receive data.') % rank
                    darrsizen = array([0,0],dtype='i')
                    reqdsizen = comm.Irecv(darrsizen, source=MPI.ANY_SOURCE, tag=30)

                    # Waiting for data size
                    dwaitsn = 1

                # Intercomm: If data size received
                if(reqdsizen.Test()):

                    # If not waiting to receive data
                    if(dwaitn == 0):

                        darrn = empty((darrsizen[0],D+2))
                        reqdsn = comm.Irecv(darrn,source=darrsizen[1], tag=28)

                        # Waiting for data
                        dwaitn = 1

                    # If data received
                    if(reqdsn.Test()):

                        for i in range(0,darrsizen[0]):
                            ctemp = circle(darrn[i,2:2+D],darrn[i,1])
                            ctemp.lbound = darrn[i,0]
                            add_task(ctemp)

                        print('Intercomm: Data received succesfully.')

                        # Not waiting for data size or data
                        dwaitsn = 0
                        dwaitn = 0

                # Debug Output
                if (SD == 1):
                    #if(rank == 0):
                    print('Time: %f') % (time()),
                    print('Processor %i Number of elements: %i') % (rank,len(clist))
                    print('U: %f') % U
                    print('LL: %f') % L
                    #print('Circle radius: %e') % rad
                    print('----------------------------')

            # Run out of work
            else:

                    # If not waiting to receive data size
                if(dwaits == 0):

                    print('Time: %f') % (time()),
                    print('Processor %i ran out of work!') % rank

                    # Receive data
                    print('Processor %i asynchronously waiting to receive data.') % rank
                    darrsize = array([0,0],dtype='i')
                    reqdsize = comm.Irecv(darrsize, source=MPI.ANY_SOURCE, tag=11)

                    # Waiting for data size
                    dwaits = 1

                # If data size received
                if(reqdsize.Test()):

                    # If not waiting to receive data
                    if(dwait == 0):

                        darr = empty((darrsize[0],D+2))
                        reqds = comm.Irecv(darr,source=darrsize[1], tag=1)

                        # Waiting for data
                        dwait = 1

                    # If data received
                    if(reqds.Test()):

                        for i in range(0,darrsize[0]):
                            ctemp = circle(darr[i,2:2+D],darr[i,1])
                            ctemp.lbound = darr[i,0]
                            add_task(ctemp)

                            if(Vis ==1):
                                drawc(ctemp)

                        print('Data received succesfully.')

                        # Not waiting for data size or data
                        dwaits = 0
                        dwait = 0

                # Intercomm: If not waiting to receive data size
                if(dwaitsn == 0):

                    # Receive data
                    print('Intercomm: Processor %i asynchronously waiting to receive data.') % rank
                    darrsizen = array([0,0],dtype='i')
                    reqdsizen = comm.Irecv(darrsizen, source=MPI.ANY_SOURCE, tag=30)

                    # Waiting for data size
                    dwaitsn = 1

                # Intercomm: If data size received
                if(reqdsizen.Test()):

                    # If not waiting to receive data
                    if(dwaitn == 0):

                        darrn = empty((darrsizen[0],D+2))
                        reqdsn = comm.Irecv(darrn,source=darrsizen[1], tag=28)

                        # Waiting for data
                        dwaitn = 1

                    # If data received
                    if(reqdsn.Test()):

                        for i in range(0,darrsizen[0]):
                            ctemp = circle(darrn[i,2:2+D],darrn[i,1])
                            ctemp.lbound = darrn[i,0]
                            add_task(ctemp)

                        print('Intercomm: Data received succesfully.')

                        # Not waiting for data size or data
                        dwaitsn = 0
                        dwaitn = 0

            # Asynchronously send UL and len(clist)
            # If previous communication has finished start new
            if(reqbsend.Test()):

                # Send UL and len(clist)
                reqbsend = comm.Issend(array([U,len(clist)]), dest=0, tag=2)

            if(reqbrecv.Test()):

                # Update U
                if(rarr[0] < U):
                    U = rarr[0]

                # New Receive
                reqbrecv = comm.Irecv(rarr, source=0, tag=3)

            # Asynchronously receive kill signal
            if(reqk.Test()):
                kill = 1

        # Say if finished!
        print('Processor %i has finished!') % rank

        # Clean up redundant messages
        #reqbsend.Cancel()
        reqbrecv.Cancel()
        #reqd.Cancel()
        #reqbsend.Free()
        reqbrecv.Free()
        #reqd.Free()

        # Send xopt to master
        print('Processor %i sending final data to master') % rank
        if('xopt' not in locals()):  # No xopt
            xopt = inf*ones(D)
        comm.Ssend(xopt, dest=0, tag=5)

        # Display figures and wait
        if(Vis == 1):
            show()

        return None, None, None, None
