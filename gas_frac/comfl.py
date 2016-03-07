"""
This base code is from p367
Matches Sod shock tube exact solution

"""

import numpy as np
import pylab as pl

def plot(x):
    pl.plot(x)

#inputs
N  = 1000
dl = 2.0
dt = 4e-7
endtime = 1.75e-3
RG = 287.0
gamma = 1.4

p0l = 150000.0
p0r = 100000.0
it0  = 300.0
iv0  = 0.0

x = np.linspace(0,dl,N)
r0,v0,t0 = np.zeros(N), np.zeros(N), np.zeros(N)

v0.fill(iv0)
t0.fill(it0)
p0 = np.zeros(N)
p0[x<dl/2.0]=p0l
p0[x>=dl/2.0]=p0r

r0=p0/(RG*t0)

time=0
dx=dl/(N-1.0)


p=p0; t=t0; r=r0; v=v0

stop=int(endtime/dt)
i=0
while i<stop:#<endtime:
    i+=1
    #calculate the derivatives.
    # for each interior point
    #calculate pressure from EoS Here.

    rl, rr, rc = np.roll(r,1)[1:-1], np.roll(r,-1)[1:-1],  r[1:-1]
    vl, vr, vc = np.roll(v,1)[1:-1], np.roll(v,-1)[1:-1],  v[1:-1]
    pll, pr    = np.roll(p,1)[1:-1], np.roll(p,-1)[1:-1]
    tl, tr, tc = np.roll(t,1)[1:-1], np.roll(t,-1)[1:-1],  t[1:-1]

    # initial derivatives
    dr0 = -( rr*vr-rl*vl)/2.0/dx
    dv0 = -(pr-pll)/(2.0*rc*dx) - vc*((vc-vl)/dx)
    dt0 = -(gamma-1)*tc*(vr-vl)/2.0/dx - vc*((tc-tl)/dx)

    r1 = r.copy()
    v1 = v.copy()
    t1 = t.copy()

    # prediction step
    r1[1:-1] += dr0*dt
    v1[1:-1] += dv0*dt
    t1[1:-1] += dt0*dt
    p1= r1*RG*t1


    r1l, r1r, r1c = np.roll(r1,1)[1:-1], np.roll(r1,-1)[1:-1],  r1[1:-1]
    v1l, v1r, v1c = np.roll(v1,1)[1:-1], np.roll(v1,-1)[1:-1],  v1[1:-1]
    p1l, p1r      = np.roll(p1,1)[1:-1], np.roll(p1,-1)[1:-1]
    t1l, t1r, t1c = np.roll(t1,1)[1:-1], np.roll(t1,-1)[1:-1],  t1[1:-1]

    # new derivatives
    dr1 = -( r1r*v1r-r1l*v1l)/2.0/dx
    dv1 = -(p1r-p1l)/(2.0*r1c*dx) - v1c*((v1c-v1l)/dx)
    dt1 = -(gamma-1)*t1c*(v1r-v1l)/2.0/dx - v1c*((t1c-t1l)/dx)

    # corrector step
    r[1:-1] += 0.5*(dr0+dr1)*dt
    v[1:-1] += 0.5*(dv0+dv1)*dt
    t[1:-1] += 0.5*(dt0+dt1)*dt
    p=r*RG*t


pl.subplot(411)
pl.plot(x,p)
pl.ylabel('Pa')
pl.subplot(412)
pl.plot(x,t)
pl.ylabel('K')
pl.subplot(413)
pl.plot(x,v)
pl.ylabel('m/s')
pl.subplot(414)
pl.plot(x,r)
pl.ylabel('kg/m^3')
pl.show()
