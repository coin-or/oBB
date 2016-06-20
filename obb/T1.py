from __future__ import division
# Overlapping Branch and Bound, data parallel algorithm
def runpar(f, g, H, Lg, Lh, l, u, bound, circle, A=None, b=None, E=None, d=None, Tol=1e-2, Heur=0, TolType='r', Vis=0, SD=1, Rtol=1e-15, TimeQP=0, TimeWidle=0, qpsolver='cvxopt'):

    # Optional Inputs
    # Tolerance
    # Heuristic lattice (0 - off, 1 - on)
    # Tolerance type (r - relative, a - absolute)
    # Visualisation (0 - off, 1 - on)
    # Step Debugging (0 - off, 1 - on)
    # Radius Tolerance
    # Time Worker Idling (0 - off, 1 - on)
    # Time QP solves (0 - off, 1 - on)
    # QP Solver (quadprog, cvxopt, nag)

    # MPI
    from mpi4py import MPI

    # MPI comm
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Get D
    D = len(l)

    # Master Process
    if(rank == 0):

        # Number of processors
        numprocs = comm.Get_size()

        # Import necessary functions
        from numpy import array, pi, sqrt, dot, identity, hstack, vstack, empty, inf
        from numpy.linalg import norm
        from time import time
        from itertools import product
        from sys import float_info

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

            # Draw bound constraints [l,u]
            fig = figure('Processor '+str(rank))
            gca().add_patch(Rectangle((l[0],l[1]), u[0]-l[0], u[1]-l[1], fill=False))
            axis([l[0]-1,u[0]+1,l[1]-1,u[1]+1])
            title('Master Processor')
            show(block=False)

            # Circle drawing procedure
            def drawc(c,col):

                # Draw circle
                gca().add_patch(Circle((c.xc[0],c.xc[1]), radius=c.r , color=col, fill=False))
                axis('equal')
                draw()

        if(Heur == 0):

            ksp = 3**D

        else:
            # Load relevant normalised lattice
            from numpy import loadtxt
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

        # Set up initial circle
        c0 = circle((u+l)/2,norm(u-l)/2)
        c0.xopt = c0.xc

        # Bound circle
        c0.lbound = bound(c0,Lg,Lh,f,g,H,D)

        # Upper bound
        c0.ubound = f(c0.xopt)

        # Set up circle list
        clist = [c0]
        cslb = clist[0]

        # Update global bounds
        U = cslb.ubound
        L = cslb.lbound
        xopt = cslb.xopt
        rad = cslb.r

        # Loop counter
        itr = 0
        timer = time()
        if(TimeQP == 1):
            qptime = 0

        # Debug Output
        if (SD == 1):
            print('----------------------------')
            print('Number of elements: %i') % len(clist)
            print('U: %f') % U
            print('L: %f') % L
            print('Circle radius: %e') % rad
            print('----------------------------')

        # Set tolerance
        if((TolType == 'r')and(abs(U) > float_info.epsilon)):
            cutoff = (U - L)/abs(U)
        else:
            cutoff = U - L

        #and((time()-timer) < 3000)
        while((cutoff > Tol)and(rad > Rtol)):

            # Update iteration count
            itr = itr + 1

            # Prune list
            i = 0
            while(i < len(clist)):
                if(clist[i].lbound > U):

                         # Visualise
                    if(Vis == 1):
                        drawc(clist[i],'k')

                    del clist[i]
                    i=i-1
                i = i+1

            if(Heur == 0):

                # Split circle into more circles
                inc = (cslb.r)/sqrt(D) # Increment
                rn = (cslb.r)/2 # New radius
                xc = cslb.xc # Centre

                # Visualise
                if(Vis == 1):
                    drawc(cslb,'k')

                clist.remove(cslb) # Remove split circle from list

                # Create square spoke configuration
                spk = array([p for p in product([-inc,0,inc], repeat=D)])

            else:

                # Split circle into more circles
                inc = cslb.r # Increment
                rn = (cslb.r)/2 # New radius
                xc = cslb.xc # Centre

                # Visualise
                if(Vis == 1):
                    drawc(cslb,'k')

                clist.remove(cslb) # Remove split circle from list

                # Scale configuration
                spk = inc*lat

            # List to distribute amongst processes
            dlist = []

            # Create surrounding circles
            for i in range(0,ksp):

                # Create centre
                xcn = xc + spk[i,:]

                # Check if circle exists
                nc = 0;
                for k in range(0,len(clist)):
                    if(all(xcn == clist[k].xc)):
                        nc = 1

                # If circle doesn't exist
                if(nc == 0):

                    # Create circle
                    cn = circle(xcn,rn)

                    # Time QP Solve
                    if(TimeQP == 1):
                        eltime = MPI.Wtime()
                    mfeas, cn.xopt = mfeasible(cn)
                    if(TimeQP == 1):
                        qptime += MPI.Wtime() - eltime

                    # If circle has a feasible point (i.e. it's feasible)
                    if(mfeas != 0):

                        # Add to distribution list
                        dlist.append(cn)

                        # Visualise
                        if(Vis == 1):
                            drawc(cn,'b')

            # Distribute data evenly amongst processes
            ihi = 0
            trem = len(dlist)
            prem = numprocs
            req = []

            for p in range(0,numprocs-1):

                tproc = int(round(trem/prem))

                ilo = ihi + 1
                ihi = ihi + tproc

                req.append(comm.isend(dlist[ilo-1:ihi], dest=numprocs-1-p, tag=0))

                prem = prem - 1
                trem = trem - tproc

            # Distribute remaining data to self
            tproc = int(round(trem/prem))

            # Bound each item allocated to self
            for k in range(ihi,ihi+tproc):

                # Bound circle
                dlist[k].ubound = f(dlist[k].xopt)

                dlist[k].lbound = bound(dlist[k],Lg,Lh,f,g,H,D)

                # Add to circle list
                clist.append(dlist[k])

            # Gather data back up
            for p in range(0,numprocs-1):

                # Make sure data has been sent
                req[p].Wait()

                # Get list of bounded circles from other processes
                dlist = comm.recv(source=numprocs-1-p, tag=1)

                # Add to clist
                clist += dlist

            # Find circle c with smallest ubound
            cslb = min(clist,key=lambda x: x.ubound)

            # Update global feasible upper bound
            U = cslb.ubound
            xopt = cslb.xopt
            rad = cslb.r

            # Find circle with smallest lbound
            cslb = min(clist,key=lambda x: x.lbound)
            L = cslb.lbound

            # Set tolerance
            if((TolType == 'r')and(abs(U) > float_info.epsilon)):
                cutoff = (U - L)/abs(U)
            else:
                cutoff = U - L

            # Debug Output
            if (SD == 1):
                print('Number of elements: %i') % len(clist)
                print('U: %f') % U
                print('L: %f') % L
                print('Circle radius: %e') % rad
                print('----------------------------')

        # Output end result
        print('Minimum value of %f at (') % f(xopt),
        for i in range(0,D-1):
            print('%f,') % xopt[i],
        print('%f)') % xopt[D-1]
        print('found with'),
        if((TolType == 'r')and(abs(U) > float_info.epsilon)):
            print('relative'),
        else:
            print('absolute'),
        print('tolerance %f in %i iterations.') % (cutoff, itr)
        tol = cutoff
        xs = xopt
        fxs = f(xopt)
        print('Elapsed time %f seconds.') % (time()-timer)
        if(TimeQP == 1):
            print('Time taken for QP solves is %f seconds') % qptime

        # Kill worker processes
        for p in range(0,numprocs-1):

            comm.send(None, dest=numprocs-1-p, tag=0)

        # Display figures and wait
        if(Vis == 1):
            show()

        return xs, fxs, tol, itr

    # Worker processes
    else:

        # Idle time
        if(TimeWidle == 1):
            eltime = MPI.Wtime()
            itime = 0

        # Pick up scattered list
        rlist = comm.recv(source=0, tag=0)

        # Idle time
        if(TimeWidle == 1):
            itime += MPI.Wtime() - eltime

        while(rlist != None):

            # Bound each item in the list
            for k in range(0,len(rlist)):

                # Bound circle
                rlist[k].ubound = f(rlist[k].xopt)

                rlist[k].lbound = bound(rlist[k],Lg,Lh,f,g,H,D)

            # Send bounded list
            comm.send(rlist, dest=0, tag=1)

            # Idle time
            if(TimeWidle == 1):
                eltime = MPI.Wtime()

            # Pick up next scattered list
            rlist = comm.recv(source=0, tag=0)

            # Idle time
            if(TimeWidle == 1):
                itime += MPI.Wtime() - eltime

        # Output idle time
        if(TimeWidle == 1):
            print('Processor %i has been idle for %f seconds') % (rank,itime)
        return None, None, None, None
