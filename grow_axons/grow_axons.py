import numpy as np

def distance(A,P,Q):
    num = (Q[0]-P[0])*(P[1]-A[1]) - (P[0]-A[0])*(Q[1]-P[1])
    denom = np.sqrt((Q[0]-P[0])**2+(Q[1]-P[1])**2)
    return num/denom

def get_slope(A,B):
    return np.arctan2(B[1]-A[1],B[0]-A[0])

def get_angle(A,B,P,Q):
    ''' angle of AB w.resp. to PQ '''
    return np.angle(np.exp(1j*(get_slope(A,B)-get_slope(P,Q))))

def get_normal(P,Q):
    aPQ = get_slope(P,Q)
    return aPQ + 0.5*np.pi

def ipoint(A,B,P,Q):
    Ax,Ay = A; Bx,By = B
    Px,Py = P; Qx,Qy = Q

    D = (Px-Qx)*(Ay-By) - (Py-Qy)*(Ax-Bx)
    x = ((Px*Qy-Py*Qx)*(Ax-Bx)-(Px-Qx)*(Ax*By-Ay*Bx)) / D
    y = ((Px*Qy-Py*Qx)*(Ay-By)-(Py-Qy)*(Ax*By-Ay*Bx)) / D

    return np.array([x,y])
    
def get_rect(x,y, w_d,w_u):
    ''' d, u, d, u, ... '''
    w = w_d+w_u
    col = 2*int(x/w) + ((w*(x/float(w)-int(x/w))) >= w_d)
    row = 2*int(y/w) + ((w*(y/float(w)-int(y/w))) >= w_d)

    return row,col

def get_H(row,col,H):
    if len(H) == 0:
        return 0
    if row<0: row = 0
    if row >= np.shape(H)[0]: row = np.shape(H)[0]-1
    if col<0: col = 0
    if col >= np.shape(H)[1]: col = np.shape(H)[1]-1

    return H[row,col]
        
def grow_NC_grid(X0,Y0,H=[],h=0.0,Pup=1.0,Pdown=1.0,w_d=300e-3,w_u=200e-3,
                Pe=0.8, alphaE=0.2, alphaI=0.4, Dl=10e-3, phi_sd=0.1,
                r_d_mu_E=150e-3, r_d_sd_E=20e-3, r_d_mu_I=150e-3, r_d_sd_I=20e-3,
                L_mu_E=1.0/np.sqrt(np.pi/2.0), L_mu_I=1.0/np.sqrt(np.pi/2.0)):
    M = len(X0)
    Me= int(Pe*M)

    r_d = np.append(
        r_d_mu_E + r_d_sd_E*np.random.normal(0,1,Me),
        r_d_mu_I + r_d_sd_I*np.random.normal(0,1,M-Me))

    Ln = np.mat(np.append(
        L_mu_E*np.sqrt(-2*np.log(1-np.random.rand(Me))),
        L_mu_I*np.sqrt(-2*np.log(1-np.random.rand(M-Me))))).T

    Nl = int(np.max(Ln)/Dl)
    Xi = np.mat(np.zeros((M,Nl))); Xi[:,0] = np.mat(X0).T
    Yi = np.mat(np.zeros((M,Nl))); Yi[:,0] = np.mat(Y0).T

    phi = np.mat(2*np.pi*np.random.rand(M,1))

    W = np.mat(np.zeros((M,M)))
    Wp = np.mat(np.zeros((M,M)),dtype=int)

    for n in np.arange(1,Nl):
        print('%d / %d' % (n,Nl))

        Dx = np.multiply(Dl*np.cos(phi), (n*Dl)<Ln)
        Dy = np.multiply(Dl*np.sin(phi), (n*Dl)<Ln)


        xi = Xi[:,n-1]; yi = Yi[:,n-1]
        xj = xi + Dx;   yj = yi + Dy

        for i in np.arange(M)[np.array((n*Dl)<Ln).T[0]]:
            if len(H) > 0:
                re_check = True

                while re_check:
                    re_check = False
                    Arow,Acol = get_rect(xi[i,0],yi[i,0],w_d,w_u)
                    Brow,Bcol = get_rect(xj[i,0],yj[i,0],w_d,w_u)

                    Ha = get_H(Arow,Acol,H)
                    Hb = get_H(Brow,Bcol,H)
                    Pa = Pup if Ha == 1 else Pdown
                    Pb = Pup if Hb == 1 else Pdown

                    if Ha != Hb and Pb != 1:
                        A = np.array([xi[i,0],yi[i,0]])
                        B = np.array([xj[i,0],yj[i,0]])
                        if Acol < Bcol:
                            xborder = w_d+int(Acol/2)*(w_d+w_u) + w_u*(Acol%2)
                        elif Acol > Bcol:
                            xborder = w_d+int(Bcol/2)*(w_d+w_u) + w_u*(Bcol%2)
                        else:
                            xborder = -1
                        if Arow < Brow:
                            yborder = w_d+int(Arow/2)*(w_d+w_u) + w_u*(Arow%2)
                        elif Arow > Brow:
                            yborder = w_d+int(Brow/2)*(w_d+w_u) + w_u*(Brow%2)
                        else:
                            yborder = -1

                        if xborder >= 0 and yborder >= 0:
                            Px = np.array([xborder,-100]) 
                            Qx = np.array([xborder,100]) 
                            Py = np.array([-100,yborder]) 
                            Qy = np.array([100,yborder]) 
                            P,Q = (Py,Qy) if distance(A,Py,Qy) < distance(A,Px,Qx) else (Px,Qx)
                        elif xborder >= 0:
                            P = np.array([xborder,-100]) 
                            Q = np.array([xborder,100]) 
                        elif yborder >= 0:
                            P = np.array([-100,yborder]) 
                            Q = np.array([100,yborder]) 
                        else:
                            pass

                        aAB_onPQ = get_angle(A,B,P,Q) 
                        if Ha!=Hb and ( (abs(aAB_onPQ) <= (np.pi/6.)) or (np.random.rand() <= (1-Pb)) ):
                            '''
                            R = line_funcs.ipoint(P,Q,A,B)
                            DlAR = np.sqrt((A[0]-R[0])**2+(A[1]-R[1])**2)
                            phi[i,0] = line_funcs.deflected_angle(P,Q,A,B,max_angle=30)
                            phi[i,0]+= phi_sd*np.random.normal(0,1)
                            xj[i,0] = R[0] + (Dl-DlAR)*np.cos(phi[i,0])
                            yj[i,0] = R[1] + (Dl-DlAR)*np.sin(phi[i,0])
                            '''

                            if  aAB_onPQ < 0:
                                phi[i,0] = get_slope(P,Q)# + phi_sd*abs(np.random.normal(0,1))
                            else:
                                phi[i,0] = get_slope(Q,P)# + phi_sd*abs(np.random.normal(0,1))
                            xj[i,0] = xi[i,0] + Dl*np.cos(phi[i,0])
                            yj[i,0] = yi[i,0] + Dl*np.sin(phi[i,0])
                            re_check = True
                        else:
                            Ln[i,0] -= abs(h)

            P = np.sqrt(np.power(xj[i,0]-X0,2)+np.power(yj[i,0]-Y0,2)) < r_d
            P[i] = 0
            for j in np.arange(M)[P==1]:
                alpha = alphaE if j<Me else alphaI
                if len(H) == 0:
                    W[j,i] = np.random.rand() > (1-alpha)
                    Wp[j,i] = 1
                else:
                    Arow,Acol = get_rect(X0[j],Y0[j],w_d,w_u)
                    Brow,Bcol = get_rect(xj[i,0],yj[i,0],w_d,w_u)

                    Ha = get_H(Arow,Acol,H)
                    Hb = get_H(Brow,Bcol,H)
                    Pa = Pup if Ha == 1 else Pdown
                    Pb = Pup if Hb == 1 else Pdown

                    if (Ha == Hb) or (np.random.rand() > (1-Pb)): # has a problem?
                        W[j,i] = np.random.rand() > (1-alpha)
                        Wp[j,i] = 1

        Xi[:,n] = xj; Yi[:,n] = yj

        phi += phi_sd*np.random.normal(0,1,(M,1))


    W = 1*(W>0); W[:,Me:] *= -1
    Wp[:,Me:] *= -1
    return W,Wp,Xi,Yi
