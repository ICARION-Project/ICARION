## GPU Field Interpolation Plan

### Scope
- Add a GPU-backed field interpolator for grid-based fields (field maps, space-charge grids).
- Keep analytical fields on CPU (overhead not worth it).
- No GPU DomainManager/boundaries; only field lookup is offloaded.

### Objectives
1) Unified interface: `IFieldInterpolator` with CPU and GPU implementations.
2) GPU path uses texture/3D grid sampling; one-time host→device upload.
3) Wire into ForceContext/ForceRegistry so GPU forces can query device-resident fields.
4) Parity/integration tests to ensure CPU/GPU agreement.

### Design
- Interface: `IFieldInterpolator` (method `Vec3 sample(const Vec3& pos_local, double t)` or CUDA kernel equivalent).
- CPU impl wraps existing FieldProviderModel/FieldModel.
- GPU impl (`GPUFieldInterpolator`) in `src/core/gpu/fields/`, owns device grid/texture, spacing, transforms.
- Factory: when `simulation.enable_gpu` and grid data present, create GPU interpolator; otherwise CPU.
- ForceContext: extend with optional `gpu_field_interpolator` pointer for GPU paths; CPU ignores.
- Electric-field GPU kernels use the GPU interpolator; space-charge GPU grid sampling can reuse it if grid-based.

### Data Flow
1) Load field map on CPU (existing loaders).
2) If GPU enabled and grid available:
   - Upload grid to device (cudaArray/texture).
   - Store dims/spacing/origin and local transforms (from DomainManager).
3) Runtime:
   - CPU forces use CPU interpolator.
   - GPU kernels fetch E via GPU interpolator (no per-step host transfer).

### Tests
- Unit parity: CPU vs GPU interpolator on synthetic linear field map (E = ax + by + cz).
-,Integration parity: small ensemble with E-only force CPU vs GPU.
- Space-charge grid parity: reuse existing grid tests comparing CPU vs GPU sampling if grids exist.

### Notes
- RNG/physics unchanged; only field sampling backend differs.
- AnalyticalFieldModel stays CPU.
- Focus is on ROI for grid interpolation; analytical fields remain CPU.
