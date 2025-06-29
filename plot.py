import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# -- SETTINGS --
mp.mp.dps = 50  # high precision
n_phi, n_theta = 80, 80  # mesh resolution

# Ellipsoid parameters
a, b, c = 1.5, 1.0, 0.8
A = np.array([a, 0.0, 0.0])
B = np.array([0.0, b, 0.0])
C = np.array([0.0, 0.0, c])

# -- Numeric integration of geodesic patch volume with error estimate --
def r(phi, theta):
    return [a*mp.sin(phi)*mp.cos(theta),
            b*mp.sin(phi)*mp.sin(theta),
            c*mp.cos(phi)]
def r_phi(phi, theta):
    return [a*mp.cos(phi)*mp.cos(theta),
            b*mp.cos(phi)*mp.sin(theta),
           -c*mp.sin(phi)]
def r_theta(phi, theta):
    return [-a*mp.sin(phi)*mp.sin(theta),
             b*mp.sin(phi)*mp.cos(theta),
             mp.mpf('0')]
def integrand(phi, theta):
    R = r(phi, theta)
    Rp = r_phi(phi, theta); Rt = r_theta(phi, theta)
    cross = [Rp[1]*Rt[2] - Rp[2]*Rt[1],
             Rp[2]*Rt[0] - Rp[0]*Rt[2],
             Rp[0]*Rt[1] - Rp[1]*Rt[0]]
    return (R[0]*cross[0] + R[1]*cross[1] + R[2]*cross[2]) / 3

def inner(theta):
    return mp.quad(lambda phi: integrand(phi, theta), [0, mp.pi/2])

# Outer integral with error estimate
V_patch, err_patch = mp.quad(inner, [0, mp.pi/2], error=True)
V_plane = abs(np.dot(A, np.cross(B, C))) / 6
V_diff = V_patch - V_plane
rel_err_pct = abs(err_patch / V_patch) * 100

# -- Build geodesic patch mesh ----------
phis = np.linspace(0, np.pi/2, n_phi)
thetas = np.linspace(0, np.pi/2, n_theta)
Phi, Theta = np.meshgrid(phis, thetas)
X = a * np.sin(Phi) * np.cos(Theta)
Y = b * np.sin(Phi) * np.sin(Theta)
Z = c * np.cos(Phi)

faces = []
for i in range(n_theta-1):
    for j in range(n_phi-1):
        p00 = (X[i,j], Y[i,j], Z[i,j])
        p10 = (X[i+1,j], Y[i+1,j], Z[i+1,j])
        p01 = (X[i,j+1], Y[i,j+1], Z[i,j+1])
        p11 = (X[i+1,j+1],Y[i+1,j+1],Z[i+1,j+1])
        faces.append([p00,p10,p11]); faces.append([p00,p11,p01])
faces.append([tuple(A),tuple(B),tuple(C)])

def build_wall(P0,P1,edge_pts):
    t = np.linspace(0,1,len(edge_pts))
    P_edge = np.outer(1-t,P0)+np.outer(t,P1)
    for k in range(len(t)-1):
        p0,p1 = tuple(P_edge[k]),tuple(P_edge[k+1])
        r0,r1 = tuple(edge_pts[k]),tuple(edge_pts[k+1])
        faces.append([p0,p1,r1]); faces.append([p0,r1,r0])

AB_edge = list(zip(a*np.cos(thetas), b*np.sin(thetas), np.zeros_like(thetas)))
BC_edge = list(zip(np.zeros_like(phis), b*np.cos(phis), c*np.sin(phis)))
CA_edge = list(zip(a*np.cos(phis), np.zeros_like(phis), c*np.sin(phis)))
build_wall(A,B,AB_edge); build_wall(B,C,BC_edge); build_wall(C,A,CA_edge)

# -- Plot and annotate --------------
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(111, projection='3d')

# ellipsoid wireframe
u_s = np.linspace(0,2*np.pi,80); v_s = np.linspace(0,np.pi,40)
um, vm = np.meshgrid(u_s, v_s)
Xe = a*np.sin(vm)*np.cos(um); Ye = b*np.sin(vm)*np.sin(um); Ze = c*np.cos(vm)
ax.plot_wireframe(Xe,Ye,Ze,rstride=6,cstride=6,color='gray',alpha=0.3)

# filled region polyhedron
poly = Poly3DCollection(faces, facecolors='red', edgecolors='k', linewidths=0.1, alpha=0.5)
ax.add_collection3d(poly)

# annotation moved to bottom-left
label = (
    f"V_patch  = {mp.nstr(V_patch, 12)}\n"
    f"V_base   = {mp.nstr(V_plane, 12)}\n"
    f"ΔV       = {mp.nstr(V_diff, 12)}\n"
    f"RelErr % = {mp.nstr(rel_err_pct, 2)}%"
)
ax.text2D(0.02, 0.05, label, transform=ax.transAxes, fontsize=12,
          bbox=dict(facecolor='white', alpha=0.7))

ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
ax.set_title('Accurate Volume: True Geodesic Patch vs Flat Triangle')
ax.set_box_aspect([a,b,c])
plt.tight_layout()
plt.show()
