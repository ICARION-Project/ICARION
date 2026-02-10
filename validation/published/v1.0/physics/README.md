# v1.0 physics validation bundle

Place exported physics validation artifacts here.

Expected contents (small files only):
- manifest.json
- README.md
- logs/
- figures/
- configs/

Export example:

```bash
./validation/scripts/export_run_bundle.sh \
  --run-dir validation/runs/v1.0_physics \
  --out-dir validation/published/v1.0/physics
```
