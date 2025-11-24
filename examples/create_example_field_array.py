#!/usr/bin/env python3
"""
Create example field array HDF5 files for testing field_array_terms loading.

Generates simple test fields:
1. uniform_field.h5 - Uniform electric field (for validation)
2. linear_gradient.h5 - Linear field gradient (tests interpolation)
3. dc_axial_unit.h5 - Unit voltage DC axial field (for scaling test)
"""

import numpy as np
import h5py
import os

def create_uniform_field(filename, E_value=20.0):
    """
    Create uniform electric field normalized to 1V over 50mm domain.
    
    E = 20 V/m (1V / 0.05m) everywhere in z-direction
    This represents a "unit field" that scales linearly with applied voltage:
      - Apply 10V → E = 200 V/m
      - Apply 100V → E = 2000 V/m
    
    Grid: 10x10x10 points, 10mm cubic domain
    Potential: phi = -E*z (linear in z)
    """
    print(f"Creating {filename}...")
    
    # Grid definition
    nx, ny, nz = 10, 10, 10
    xs = np.linspace(-5e-3, 5e-3, nx)  # -5mm to +5mm
    ys = np.linspace(-5e-3, 5e-3, ny)
    zs = np.linspace(0, 50e-3, nz)     # 0 to 50mm (normalized domain)
    
    # Create 3D meshgrid
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    
    # Uniform field in z-direction (axial): 1V / 50mm = 20 V/m
    Ex = np.zeros_like(X)
    Ey = np.zeros_like(X)
    Ez = np.full_like(X, E_value)
    
    # Potential: phi = -E*z (ranges from 0V at z=0 to -1V at z=50mm)
    phi = -E_value * Z
    
    # Save to HDF5
    with h5py.File(filename, 'w') as f:
        # Grid axes (1D arrays)
        f.create_dataset('x', data=xs)
        f.create_dataset('y', data=ys)
        f.create_dataset('z', data=zs)
        
        # Field components (3D arrays, flattened to 1D)
        f.create_dataset('Ex', data=Ex)
        f.create_dataset('Ey', data=Ey)
        f.create_dataset('Ez', data=Ez)
        f.create_dataset('phi', data=phi)
        
        # Metadata
        f.attrs['description'] = f'Uniform electric field E = ({E_value}, 0, 0) V/m'
        f.attrs['nx'] = nx
        f.attrs['ny'] = ny
        f.attrs['nz'] = nz
    
    print(f"  ✓ Created uniform field: E = {E_value} V/m in x-direction")
    print(f"  ✓ Grid: {nx}x{ny}x{nz}, domain: ±5mm")
    print(f"  ✓ Size: {os.path.getsize(filename)/1024:.1f} KB\n")


def create_linear_gradient(filename):
    """
    Create linear field gradient: E_x = 1000 * x V/m
    
    Tests interpolation accuracy for non-uniform fields.
    """
    print(f"Creating {filename}...")
    
    # Grid definition
    nx, ny, nz = 15, 15, 15
    xs = np.linspace(-10e-3, 10e-3, nx)  # -10mm to +10mm
    ys = np.linspace(-10e-3, 10e-3, ny)
    zs = np.linspace(-10e-3, 10e-3, nz)
    
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    
    # Linear gradient: E_x = 1000*x V/m, E_y = 500*y, E_z = 0
    Ex = 1000.0 * X
    Ey = 500.0 * Y
    Ez = np.zeros_like(X)
    
    # Potential: phi = -500*x^2 - 250*y^2
    phi = -500.0 * X**2 - 250.0 * Y**2
    
    with h5py.File(filename, 'w') as f:
        f.create_dataset('x', data=xs)
        f.create_dataset('y', data=ys)
        f.create_dataset('z', data=zs)
        f.create_dataset('Ex', data=Ex)
        f.create_dataset('Ey', data=Ey)
        f.create_dataset('Ez', data=Ez)
        f.create_dataset('phi', data=phi)
        
        f.attrs['description'] = 'Linear field gradient: Ex = 1000*x, Ey = 500*y'
        f.attrs['nx'] = nx
        f.attrs['ny'] = ny
        f.attrs['nz'] = nz
    
    print(f"  ✓ Created gradient field: Ex = 1000*x, Ey = 500*y")
    print(f"  ✓ Grid: {nx}x{ny}x{nz}, domain: ±10mm")
    print(f"  ✓ Size: {os.path.getsize(filename)/1024:.1f} KB\n")


def create_dc_axial_unit(filename):
    """
    Create unit-voltage DC axial field (1 V applied across 50mm)
    
    For testing DC_Axial scaling in field_array_terms.
    E_z = 1 V / 0.05 m = 20 V/m
    """
    print(f"Creating {filename}...")
    
    # Grid definition
    nx, ny, nz = 10, 10, 20  # Elongated in z
    xs = np.linspace(-5e-3, 5e-3, nx)   # ±5mm radial
    ys = np.linspace(-5e-3, 5e-3, ny)
    zs = np.linspace(0, 50e-3, nz)      # 0 to 50mm axial
    
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    
    # Uniform axial field: 1V / 50mm = 20 V/m
    Ex = np.zeros_like(X)
    Ey = np.zeros_like(X)
    Ez = np.full_like(X, 20.0)
    
    # Linear potential: phi = -20*z (from 0V at z=0 to -1V at z=50mm)
    phi = -20.0 * Z
    
    with h5py.File(filename, 'w') as f:
        f.create_dataset('x', data=xs)
        f.create_dataset('y', data=ys)
        f.create_dataset('z', data=zs)
        f.create_dataset('Ex', data=Ex)
        f.create_dataset('Ey', data=Ey)
        f.create_dataset('Ez', data=Ez)
        f.create_dataset('phi', data=phi)
        
        f.attrs['description'] = 'Unit DC axial field: 1V over 50mm → 20 V/m'
        f.attrs['nx'] = nx
        f.attrs['ny'] = ny
        f.attrs['nz'] = nz
        f.attrs['voltage_V'] = 1.0
        f.attrs['length_m'] = 0.05
    
    print(f"  ✓ Created DC axial unit field: Ez = 20 V/m (1V/50mm)")
    print(f"  ✓ Grid: {nx}x{ny}x{nz}, domain: ±5mm radial, 0-50mm axial")
    print(f"  ✓ Size: {os.path.getsize(filename)/1024:.1f} KB\n")


def verify_field(filename):
    """Quick verification that HDF5 file structure is correct"""
    print(f"Verifying {filename}...")
    
    with h5py.File(filename, 'r') as f:
        # Check required datasets
        required = ['xs', 'ys', 'zs', 'Ex', 'Ey', 'Ez', 'phi']
        missing = [name for name in required if name not in f]
        
        if missing:
            print(f"  ✗ Missing datasets: {missing}")
            return False
        
        # Check dimensions
        nx = len(f['xs'][:])
        ny = len(f['ys'][:])
        nz = len(f['zs'][:])
        expected_size = nx * ny * nz
        
        for field in ['Ex', 'Ey', 'Ez', 'phi']:
            actual_size = len(f[field][:])
            if actual_size != expected_size:
                print(f"  ✗ {field}: expected {expected_size} points, got {actual_size}")
                return False
        
        print(f"  ✓ All datasets present and correct size")
        print(f"  ✓ Grid: {nx}x{ny}x{nz} = {expected_size} points")
        
        # Print field statistics
        Ex = f['Ex'][:]
        Ey = f['Ey'][:]
        Ez = f['Ez'][:]
        print(f"  ✓ Ex: [{Ex.min():.2e}, {Ex.max():.2e}] V/m")
        print(f"  ✓ Ey: [{Ey.min():.2e}, {Ey.max():.2e}] V/m")
        print(f"  ✓ Ez: [{Ez.min():.2e}, {Ez.max():.2e}] V/m\n")
        
        return True


if __name__ == '__main__':
    print("=" * 60)
    print("Creating Example Field Array HDF5 Files")
    print("=" * 60 + "\n")
    
    # Create output directory if needed
    output_dir = "field_arrays"
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate test fields
    files = [
        (f"{output_dir}/uniform_field.h5", lambda: create_uniform_field(f"{output_dir}/uniform_field.h5", 1000.0)),
        (f"{output_dir}/linear_gradient.h5", lambda: create_linear_gradient(f"{output_dir}/linear_gradient.h5")),
        (f"{output_dir}/dc_axial_unit.h5", lambda: create_dc_axial_unit(f"{output_dir}/dc_axial_unit.h5")),
    ]
    
    for filename, create_func in files:
        create_func()
        verify_field(filename)
    
    print("=" * 60)
    print("✅ All field arrays created successfully!")
    print("=" * 60)
    print(f"\nFiles created in: {output_dir}/")
    print("\nUsage in config JSON:")
    print("""
{
  "fields": {
    "field_array_terms": [
      {
        "file": "field_arrays/uniform_field.h5",
        "scale_type": "Constant",
        "constant_V": 10.0
      },
      {
        "file": "field_arrays/dc_axial_unit.h5",
        "scale_type": "DC_Axial"
      }
    ]
  }
}
""")
