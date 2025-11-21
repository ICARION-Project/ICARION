
```css
simulation.h5
├── trajectories/ (all dynamic data such as position, velocity and changing ion params)
│   ├── time                    [T]
│   ├── positions               [T × N × 3]
│   ├── velocities              [T × N × 3]
│   ├── active                  [T × N]
│   ├── born                    [T × N]
│   ├── ion_id                  [T × N]
│   ├── species_id              [T × N]
│   └── domain_index            [T × N]
│
├── ions/      (static per-ion metadata)
│   ├── initial_species_id      [N]
│   ├── initial_pos             [N × 3]
│   ├── initial_vel             [N × 3]
│   ├── birth_time_s            [N]
│   └── charge_C                [N] (optional, may change via reactions)
│
├── domains/
│   ├── domain_0/
│   │   ├── instrument
│   │   ├── solver_type
│   │   ├── geometry/…
│   │   ├── env/…
│   │   ├── DC/…
│   │   ├── RF/…
│   │   ├── AC/…
│   │   └── B/…
│   └── domain_1/
│       └── …
│
└── metadata/
    ├── global_params           {json string or scalar datasets}
    ├── reproducibility           
    ├── reactions/           (all reaction metadata for the reactions loaded frmo the database)
    └── species/               (all species metadata for the species loaded frmo the database)
    
```

Reproducibility: 
```css
    ├── reproducibility/
    │   ├── global_seed
    │   ├── rng_algorithm
    │   ├── seed_scheme
    │   ├── git_hash
    │   ├── git_dirty
    │   ├── code_version
    │   ├── compiler_cxx
    │   ├── compiler_flags
    │   ├── cuda_version
    │   ├── gpu_arch
    │   ├── input_hash/
    │   │    ├── config_sha256
    │   │    ├── species_db_sha256
    │   │    ├── reaction_db_sha256
    │   │    ├── geometry_sha256
    │   │    └── ioncloud_sha256
    │   └── execution/
    │        ├── openmp_threads
    │        ├── gpu_threads_per_block
    │        ├── gpu_blocks
    │        └── parallel_scheme
    ├── system/
    │   ├── cpu_model
    │   ├── gpu_model
    │   ├── os
    │   ├── driver_version
    │   └── timestamp
```


maybe also -> run_info(user, timestamp)