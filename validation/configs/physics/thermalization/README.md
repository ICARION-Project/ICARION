# Thermalization Test Configs

**Generated:** 90 configs
**Test Matrix:**
- Temperatures: [150, 300, 1000] K
- Pressures: [0.2, 2.0, 20.0, 200.0, 2000.0] Pa
- Ion species: ['H3O+', 'PentanalH+', '2,6-DTBPH+']
- Collision models: ['HSS', 'EHSS']

**Expected Results:**
- Final temperature should match environment temperature
- Thermalization time ~ few collision times
- HSS and EHSS should give similar results
