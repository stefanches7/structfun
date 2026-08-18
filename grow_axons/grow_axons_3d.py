import numpy as np


def get_cube(x, y, z, w_d, w_u):
    """Return checkerboard-style 3D grid indices for a point."""
    w = w_d + w_u
    col = 2 * int(x / w) + ((w * (x / float(w) - int(x / w))) >= w_d)
    row = 2 * int(y / w) + ((w * (y / float(w) - int(y / w))) >= w_d)
    dep = 2 * int(z / w) + ((w * (z / float(w) - int(z / w))) >= w_d)
    return int(dep), int(row), int(col)


def get_H3(dep, row, col, H):
    if len(H) == 0:
        return 0

    if np.ndim(H) == 2:
        row = min(max(row, 0), np.shape(H)[0] - 1)
        col = min(max(col, 0), np.shape(H)[1] - 1)
        return H[row, col]

    dep = min(max(dep, 0), np.shape(H)[0] - 1)
    row = min(max(row, 0), np.shape(H)[1] - 1)
    col = min(max(col, 0), np.shape(H)[2] - 1)
    return H[dep, row, col]


def _sample_unit_vectors(M):
    phi = 2.0 * np.pi * np.random.rand(M)
    cos_theta = 2.0 * np.random.rand(M) - 1.0
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    return np.column_stack(
        (
            sin_theta * np.cos(phi),
            sin_theta * np.sin(phi),
            cos_theta,
        )
    )


def grow_NC_grid_3d(
    X0,
    Y0,
    Z0,
    H=[],
    h=0.0,
    Pup=1.0,
    Pdown=1.0,
    w_d=300e-3,
    w_u=200e-3,
    Pe=0.8,
    alphaE=0.2,
    alphaI=0.4,
    Dl=10e-3,
    phi_sd=0.1,
    r_d_mu_E=150e-3,
    r_d_sd_E=20e-3,
    r_d_mu_I=150e-3,
    r_d_sd_I=20e-3,
    L_mu_E=1.0 / np.sqrt(np.pi / 2.0),
    L_mu_I=1.0 / np.sqrt(np.pi / 2.0),
    soma_diameter=15e-3,
):
    X0 = np.asarray(X0, dtype=float)
    Y0 = np.asarray(Y0, dtype=float)
    Z0 = np.asarray(Z0, dtype=float)
    M = len(X0)
    Me = int(Pe * M)

    r_d = np.append(
        r_d_mu_E + r_d_sd_E * np.random.normal(0, 1, Me),
        r_d_mu_I + r_d_sd_I * np.random.normal(0, 1, M - Me),
    )
    r_d = np.maximum(r_d, 0.0)

    Ln = np.append(
        np.random.rayleigh(scale=L_mu_E, size=Me),
        np.random.rayleigh(scale=L_mu_I, size=M - Me),
    )

    Nl = max(int(np.max(Ln) / Dl), 1)
    Xi = np.zeros((M, Nl))
    Yi = np.zeros((M, Nl))
    Zi = np.zeros((M, Nl))
    Xi[:, 0] = X0
    Yi[:, 0] = Y0
    Zi[:, 0] = Z0

    dirs = _sample_unit_vectors(M)

    W = np.zeros((M, M))
    Wp = np.zeros((M, M), dtype=int)
    soma_radius = 0.5 * soma_diameter

    for n in np.arange(1, Nl):
        print("%d / %d" % (n, Nl))
        active = (n * Dl) < Ln

        Dxyz = Dl * dirs
        Dxyz *= active[:, None]

        xi = Xi[:, n - 1]
        yi = Yi[:, n - 1]
        zi = Zi[:, n - 1]
        xj = xi + Dxyz[:, 0]
        yj = yi + Dxyz[:, 1]
        zj = zi + Dxyz[:, 2]

        for i in np.where(active)[0]:
            if len(H) > 0:
                Adep, Arow, Acol = get_cube(xi[i], yi[i], zi[i], w_d, w_u)
                Bdep, Brow, Bcol = get_cube(xj[i], yj[i], zj[i], w_d, w_u)

                Ha = get_H3(Adep, Arow, Acol, H)
                Hb = get_H3(Bdep, Brow, Bcol, H)
                Pb = Pup if Hb == 1 else Pdown

                if Ha != Hb and (np.random.rand() <= (1 - Pb)):
                    Ln[i] -= abs(h)
                    xj[i], yj[i], zj[i] = xi[i], yi[i], zi[i]

            P = (
                np.sqrt((xj[i] - X0) ** 2 + (yj[i] - Y0) ** 2 + (zj[i] - Z0) ** 2)
                < (r_d + soma_radius)
            )
            P[i] = 0
            for j in np.where(P)[0]:
                alpha = alphaE if j < Me else alphaI
                if len(H) == 0:
                    W[j, i] = np.random.rand() > (1 - alpha)
                    Wp[j, i] = 1
                else:
                    Adep, Arow, Acol = get_cube(X0[j], Y0[j], Z0[j], w_d, w_u)
                    Bdep, Brow, Bcol = get_cube(xj[i], yj[i], zj[i], w_d, w_u)

                    Ha = get_H3(Adep, Arow, Acol, H)
                    Hb = get_H3(Bdep, Brow, Bcol, H)
                    Pb = Pup if Hb == 1 else Pdown

                    if (Ha == Hb) or (np.random.rand() > (1 - Pb)):
                        W[j, i] = np.random.rand() > (1 - alpha)
                        Wp[j, i] = 1

        Xi[:, n] = xj
        Yi[:, n] = yj
        Zi[:, n] = zj

        dirs += phi_sd * np.random.normal(0, 1, dirs.shape)
        norms = np.linalg.norm(dirs, axis=1)
        norms[norms == 0.0] = 1.0
        dirs /= norms[:, None]

    W = 1 * (W > 0)
    W[:, Me:] *= -1
    Wp[:, Me:] *= -1
    return W, Wp, Xi, Yi, Zi

