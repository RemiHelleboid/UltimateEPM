"""
Generate k-points for a given type of grid (1D, 2D, 3D, BZ, etc.)
"""

import numpy as np
from argparse import ArgumentParser


def IsInIrreducibleWedge(k):
    return (k[2] >= 0.0 and k[2] <= k[1] and k[1] <= k[0] and k[0] <= 1.0) and \
        (np.sum(k) <= 3.0/2.0)


def is_in_fcc_bz(k, tol=1.0e-12):
    kx, ky, kz = k
    cond_1 = abs(kx) <= 1.0 + tol and abs(ky) <= 1.0 + \
        tol and abs(kz) <= 1.0 + tol
    cond_2 = abs(kx) + abs(ky) + abs(kz) <= 3.0 / 2.0 + tol
    return cond_1 and cond_2


def generate_3d_grid_bz_uniform(Nx, Ny, Nz, shift=[0, 0, 0], irreducible_wedge=False):
    list_kpoints = []
    for idx_x in range(Nx):
        for idx_y in range(Ny):
            for idx_z in range(Nz):
                kx = (idx_x + shift[0]) / Nx
                ky = (idx_y + shift[1]) / Ny
                kz = (idx_z + shift[2]) / Nz
                k = np.array([kx, ky, kz])
                if irreducible_wedge:
                    if IsInIrreducibleWedge(k):
                        list_kpoints.append(k)
                else:
                    if is_in_fcc_bz(k):
                        list_kpoints.append(k)
    return list_kpoints


def main():
    parser = ArgumentParser(
        description='Generate k-points for a given type of grid (1D, 2D, 3D, BZ, etc.)')
    parser.add_argument('--grid', type=str, default='3D',
                        choices=['1D', '2D', '3D', 'BZ'],
                        help='Type of grid.')
    parser.add_argument('--Nx', type=int, default=1,
                        help='Number of divisions in x-direction.')
    parser.add_argument('--Ny', type=int, default=1,
                        help='Number of divisions in y-direction.')
    parser.add_argument('--Nz', type=int, default=1,
                        help='Number of divisions in z-direction.')
    parser.add_argument('--min', type=float, default=0.0,
                        help='Minimum value of k-point.')
    parser.add_argument('--max', type=float, default=1.0, required=False,
                        help='Maximum value of k-point.')
    parser.add_argument('--sphere', action='store_true',
                        help='Keep k-points in a sphere, i.e. with a norm less than a given value.')
    parser.add_argument('--sphere_radius', type=float, default=1.0,
                        help='Radius of sphere.')
    parser.add_argument('--sphere_center', type=str, default='0,0,0',
                        help='Center of sphere.')
    parser.add_argument('-dir', '--directions', type=str, default='100',
                        help='Directions in which to generate k-points.', required=False)
    parser.add_argument('--shift', type=str, default='0,0,0',
                        help='Shift in x-, y-, and z-directions.')
    parser.add_argument('--irreducible_wedge', action='store_true',
                        help='Generate only k-points in the irreducible wedge of the BZ.')
    parser.add_argument('--output', type=str, default='kpoints.dat',
                        help='Name of output file.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    Nx, Ny, Nz = args.Nx, args.Ny, args.Nz
    min = args.min
    max = args.max
    
    if args.sphere:
        center = np.array([float(x) for x in args.sphere_center.split(',')])
        radius = args.sphere_radius 
        sphere = True
    else:
        center = np.array([0.0, 0.0, 0.0])
        radius = 1e120
        
    directions = [int(x) for x in args.directions]
    shift = [float(x) for x in args.shift.split(',')]
    print(shift)
    irreducible_wedge = args.irreducible_wedge
    output = args.output
    
    list_kpoints = []
    if args.grid == '1D':
        kx = np.linspace(min, max, Nx)
        if directions == [1, 0, 0]:
            list_kpoints = [[x, 0.0, 0.0] for x in kx if np.linalg.norm(np.array([x, 0.0, 0.0]) - center) <= radius]
        elif directions == [1, 1, 0]:
            list_kpoints = [[x, x, 0.0] for x in kx if np.linalg.norm(np.array([x, x, 0.0]) - center) <= radius]
        elif directions == [1, 1, 1]:
            list_kpoints = [[x, x, x] for x in kx if np.linalg.norm(np.array([x, x, x]) - center) <= radius]
        else:
            raise ValueError('Directions not supported.')
    elif args.grid == '2D':
        kx = np.linspace(min, max, Nx)
        ky = np.linspace(min, max, Ny)
        list_kpoints = [[x, y, 0.0] for x in kx for y in ky if np.linalg.norm(np.array([x, y, 0.0]) - center) <= radius]
        
    elif args.grid == '3D':
        kx = np.linspace(min, max, Nx)
        ky = np.linspace(min, max, Ny)
        kz = np.linspace(min, max, Nz)
        list_kpoints = [[x, y, z] for x in kx for y in ky for z in kz if np.linalg.norm(np.array([x, y, z]) - center) <= radius]
    elif args.grid == 'BZ':
        list_kpoints = generate_3d_grid_bz_uniform(Nx, Ny, Nz, shift=shift, irreducible_wedge=irreducible_wedge)
    else:
        raise ValueError('Unknown grid type: {}'.format(args.grid))
    
    # Apply shift
    list_kpoints = [np.array(k) + np.array(shift) for k in list_kpoints]
    
    print(f"Generated {len(list_kpoints)} k-points.")
    print(f"Writing k-points to {output}.")
    np.savetxt(output, list_kpoints, fmt='%.8e')
    
    
if __name__ == '__main__':
    main()