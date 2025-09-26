import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Function to cycle 4 coplanar points
# -------------------------
def cycle_quad_coplanar(P):
    P = np.asarray(P, dtype=float)  # shape (4,3)
    c = P.mean(axis=0)

    # plane normal from a non-degenerate triplet
    def plane_normal(P):
        idx = [(0,1,2), (0,1,3), (0,2,3), (1,2,3)]
        for i,j,k in idx:
            n = np.cross(P[j]-P[i], P[k]-P[i])
            if np.linalg.norm(n) > 1e-12:
                return n/np.linalg.norm(n)
        raise ValueError("Points are (nearly) collinear")

    n = plane_normal(P)

    # in-plane basis (u,v)
    u_prime = P[1]-P[0]
    u_prime -= np.dot(u_prime, n)*n
    if np.linalg.norm(u_prime) < 1e-12:
        u_prime = P[2]-P[0]
        u_prime -= np.dot(u_prime, n)*n
    u = u_prime/np.linalg.norm(u_prime)
    v = np.cross(n, u)

    # 2D projection and angle sort
    XY = (P - c) @ np.stack([u, v], axis=1)  # (4,2)
    ang = np.arctan2(XY[:,1], XY[:,0])
    order = np.argsort(ang)

    # optional: enforce right-hand winding with +n
    loop = P[order]
    nx = np.sum((loop[:,1] - np.roll(loop[:,1], -1)) * (loop[:,2] + np.roll(loop[:,2], -1)))
    ny = np.sum((loop[:,2] - np.roll(loop[:,2], -1)) * (loop[:,0] + np.roll(loop[:,0], -1)))
    nz = np.sum((loop[:,0] - np.roll(loop[:,0], -1)) * (loop[:,1] + np.roll(loop[:,1], -1)))
    poly_n = np.array([nx, ny, nz])
    if np.dot(poly_n, n) < 0:
        order = order[::-1]

    return order.tolist(), XY[order]

# -------------------------
# Demo loop
# -------------------------
np.random.seed(1)

fig, axes = plt.subplots(5, 5, figsize=(12, 12))
axes = axes.ravel()

for ax in axes:
    # generate 4 random coplanar points
    # n = np.random.randn(3); n /= np.linalg.norm(n)
    # u = np.random.randn(3); u -= np.dot(u, n) * n; u /= np.linalg.norm(u)
    # v = np.cross(n, u)
    # coeffs = np.random.randn(4, 2) * 5
    # P = coeffs[:,0][:,None]*u + coeffs[:,1][:,None]*
    A = np.random.randn(3); A[2] = 0
    B = np.random.randn(3); B[2] = 0
    C = np.random.randn(3); C[2] = 0
    D = np.random.randn(3); D[2] = 0
    P = np.array([A, B, C, D]) 

    order, XY = cycle_quad_coplanar(P)

    # plot projected 2D quadrilateral
    ax.scatter(XY[:,0], XY[:,1], color='red', zorder=5, s=5)
    # for i, (x, y) in enumerate(XY):
    #     ax.text(x, y, f"P{order[i]}", fontsize=9, ha='center', va='center')

    # close the loop for plotting
    XY_loop = np.vstack([XY, XY[0]])
    ax.plot(XY_loop[:,0], XY_loop[:,1], '-o', color='blue')

    ax.set_aspect('equal')
    # ax.set_title("Random Quad")
    ax.axis('off')

plt.tight_layout()
plt.show()
