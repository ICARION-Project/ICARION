# v1.0.0 instrument validation bundle

Place exported instrument validation artifacts here.

Expected contents (small files only):
- manifest.json
- README.md
- logs/
- figures/
- configs/

Export example:

```bash
./validation/scripts/export_run_bundle.sh \
  --run-dir validation/runs/v1.0.0_instruments \
  --out-dir validation/published/v1.0.0/instruments
```
