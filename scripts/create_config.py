#!/usr/bin/env python3
"""
ICARION Configuration Builder

Interactive tool to create ICARION JSON configuration files.
Supports templates and command-line customization.

Usage:
    python3 create_config.py                           # Interactive mode
    python3 create_config.py --template ims            # Use template
    python3 create_config.py --list-templates          # Show available templates
    python3 create_config.py --help                    # Show help
"""

import json
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Any

# Template definitions
TEMPLATES = {
    "minimal": {
        "description": "Minimal valid configuration",
        "config": {
            "simulation": {
                "total_time_s": 1e-3,
                "dt_s": 1e-9,
                "integrator": "RK4",
                "write_interval": 100
            },
            "ions": {
                "species": [
                    {
                        "id": "H3O+",
                        "count": 100,
                        "position": {
                            "type": "gaussian",
                            "center": [0.0, 0.0, 0.001],
                            "std": [0.001, 0.001, 0.0005]
                        },
                        "velocity": {
                            "type": "thermal",
                            "temperature_K": 300.0
                        }
                    }
                ]
            },
            "physics": {
                "collision_model": "NoCollisions"
            },
            "output": {
                "folder": "./results",
                "trajectory_file": "trajectories.h5",
                "print_progress": True
            },
            "domains": [
                {
                    "name": "test_domain",
                    "instrument": "IMS",
                    "geometry": {
                        "origin_m": [0.0, 0.0, 0.0],
                        "length_m": 0.1,
                        "radius_m": 0.01
                    },
                    "env": {
                        "pressure_Pa": 101325.0,
                        "temperature_K": 300.0,
                        "gas_species": "He",
                        "gas_velocity_m_s": [0.0, 0.0, 0.0]
                    }
                }
            ]
        }
    },
    
    "ims": {
        "description": "Ion Mobility Spectrometry drift tube",
        "config": {
            "simulation": {
                "total_time_s": 5e-4,
                "dt_s": 1e-8,
                "integrator": "RK4",
                "write_interval": 1000,
                "enable_gpu": True,
                "rng_seed": 42
            },
            "ions": {
                "species": [
                    {
                        "id": "H3O+",
                        "count": 1000,
                        "position": {
                            "type": "gaussian",
                            "center": [0.0, 0.0, 0.001],
                            "std": [0.001, 0.001, 0.0005]
                        },
                        "velocity": {
                            "type": "thermal",
                            "temperature_K": 300.0
                        }
                    }
                ]
            },
            "physics": {
                "collision_model": "HSS"
            },
            "output": {
                "folder": "./results/ims",
                "trajectory_file": "ims_trajectories.h5",
                "print_progress": True
            },
            "domains": [
                {
                    "name": "drift_region",
                    "instrument": "IMS",
                    "geometry": {
                        "origin_m": [0.0, 0.0, 0.0],
                        "length_m": 0.05,
                        "radius_m": 0.015
                    },
                    "env": {
                        "pressure_Pa": 200.0,
                        "temperature_K": 300.0,
                        "gas_species": "He",
                        "gas_velocity_m_s": [0.0, 0.0, 0.0]
                    },
                    "fields": {
                        "DC": {
                            "EN_Td": 10.0,
                            "axial_V": 0.0
                        }
                    }
                }
            ]
        }
    },
    
    "tof": {
        "description": "Time-of-Flight mass spectrometer",
        "config": {
            "simulation": {
                "total_time_s": 1e-4,
                "dt_s": 1e-10,
                "integrator": "RK4",
                "write_interval": 100,
                "enable_gpu": True
            },
            "ions": {
                "species": [
                    {
                        "id": "ReserpineH+",
                        "count": 100,
                        "position": {
                            "type": "gaussian",
                            "center": [0.0, 0.0, 0.001],
                            "std": [0.0005, 0.0005, 0.0005]
                        },
                        "velocity": {
                            "type": "thermal",
                            "temperature_K": 300.0
                        }
                    }
                ]
            },
            "physics": {
                "collision_model": "NoCollisions"
            },
            "output": {
                "folder": "./results/tof",
                "trajectory_file": "tof_trajectories.h5",
                "print_progress": True
            },
            "domains": [
                {
                    "name": "acceleration_region",
                    "instrument": "TOF",
                    "geometry": {
                        "origin_m": [0.0, 0.0, 0.0],
                        "length_m": 0.02,
                        "radius_m": 0.01
                    },
                    "env": {
                        "pressure_Pa": 1e-6,
                        "temperature_K": 300.0,
                        "gas_species": "He"
                    },
                    "fields": {
                        "DC": {
                            "axial_V": 5000.0
                        }
                    }
                },
                {
                    "name": "drift_region",
                    "instrument": "TOF",
                    "geometry": {
                        "origin_m": [0.0, 0.0, 0.02],
                        "length_m": 0.5,
                        "radius_m": 0.01
                    },
                    "env": {
                        "pressure_Pa": 1e-6,
                        "temperature_K": 300.0,
                        "gas_species": "He"
                    }
                }
            ]
        }
    },
    
    "lqit": {
        "description": "Linear Quadrupole Ion Trap",
        "config": {
            "simulation": {
                "total_time_s": 1e-3,
                "dt_s": 1e-9,
                "integrator": "RK4",
                "write_interval": 100,
                "enable_gpu": True
            },
            "ions": {
                "species": [
                    {
                        "id": "CaffeineH+",
                        "count": 500,
                        "position": {
                            "type": "gaussian",
                            "center": [0.0, 0.0, 0.0],
                            "std": [0.001, 0.001, 0.005]
                        },
                        "velocity": {
                            "type": "thermal",
                            "temperature_K": 300.0
                        }
                    }
                ]
            },
            "physics": {
                "collision_model": "HSS"
            },
            "output": {
                "folder": "./results/lqit",
                "trajectory_file": "lqit_trajectories.h5",
                "print_progress": True
            },
            "domains": [
                {
                    "name": "trap_region",
                    "instrument": "LQIT",
                    "geometry": {
                        "origin_m": [0.0, 0.0, 0.0],
                        "length_m": 0.05,
                        "radius_m": 0.004
                    },
                    "env": {
                        "pressure_Pa": 0.1,
                        "temperature_K": 300.0,
                        "gas_species": "He"
                    },
                    "fields": {
                        "RF": {
                            "voltage_V": 500.0,
                            "frequency_Hz": 1e6,
                            "phase_rad": 0.0
                        },
                        "DC": {
                            "quad_V": 10.0
                        }
                    }
                }
            ]
        }
    },
    
    "orbitrap": {
        "description": "Orbitrap mass analyzer",
        "config": {
            "simulation": {
                "total_time_s": 1e-3,
                "dt_s": 1e-10,
                "integrator": "RK4",
                "write_interval": 1000,
                "enable_gpu": True
            },
            "ions": {
                "species": [
                    {
                        "id": "ReserpineH+",
                        "count": 100,
                        "position": {
                            "type": "gaussian",
                            "center": [0.0, 0.0, 0.0],
                            "std": [0.002, 0.002, 0.005]
                        },
                        "velocity": {
                            "type": "thermal",
                            "temperature_K": 300.0
                        }
                    }
                ]
            },
            "physics": {
                "collision_model": "NoCollisions"
            },
            "output": {
                "folder": "./results/orbitrap",
                "trajectory_file": "orbitrap_trajectories.h5",
                "print_progress": True
            },
            "domains": [
                {
                    "name": "orbitrap_cell",
                    "instrument": "Orbitrap",
                    "geometry": {
                        "origin_m": [0.0, 0.0, 0.0],
                        "length_m": 0.05,
                        "radius_in_m": 0.006,
                        "radius_out_m": 0.015,
                        "radius_char_m": 0.022
                    },
                    "env": {
                        "pressure_Pa": 1e-8,
                        "temperature_K": 300.0,
                        "gas_species": "He"
                    }
                }
            ]
        }
    },
    
    "quadrupole": {
        "description": "Quadrupole mass filter",
        "config": {
            "simulation": {
                "total_time_s": 5e-5,
                "dt_s": 1e-9,
                "integrator": "RK4",
                "write_interval": 100,
                "enable_gpu": True
            },
            "ions": {
                "species": [
                    {
                        "id": "CaffeineH+",
                        "count": 50,
                        "position": {
                            "type": "gaussian",
                            "center": [0.0, 0.0, 0.001],
                            "std": [0.0005, 0.0005, 0.0005]
                        },
                        "velocity": {
                            "type": "thermal",
                            "temperature_K": 300.0
                        }
                    }
                ]
            },
            "physics": {
                "collision_model": "NoCollisions"
            },
            "output": {
                "folder": "./results/quadrupole",
                "trajectory_file": "quadrupole_trajectories.h5",
                "print_progress": True
            },
            "domains": [
                {
                    "name": "quad_filter",
                    "instrument": "Quadrupole",
                    "geometry": {
                        "origin_m": [0.0, 0.0, 0.0],
                        "length_m": 0.05,
                        "radius_m": 0.005
                    },
                    "env": {
                        "pressure_Pa": 0.001,
                        "temperature_K": 300.0,
                        "gas_species": "He"
                    },
                    "fields": {
                        "RF": {
                            "voltage_V": 100.0,
                            "frequency_Hz": 1e6,
                            "phase_rad": 0.0
                        },
                        "DC": {
                            "quad_V": 10.0
                        }
                    }
                }
            ]
        }
    },
    
    "fticr": {
        "description": "FT-ICR mass spectrometer",
        "config": {
            "simulation": {
                "total_time_s": 1e-3,
                "dt_s": 1e-9,
                "integrator": "Boris",
                "write_interval": 1000,
                "enable_gpu": True
            },
            "ions": {
                "species": [
                    {
                        "id": "ReserpineH+",
                        "count": 50,
                        "position": {
                            "type": "gaussian",
                            "center": [0.0, 0.0, 0.0],
                            "std": [0.001, 0.001, 0.001]
                        },
                        "velocity": {
                            "type": "thermal",
                            "temperature_K": 300.0
                        }
                    }
                ]
            },
            "physics": {
                "collision_model": "NoCollisions"
            },
            "output": {
                "folder": "./results/fticr",
                "trajectory_file": "fticr_trajectories.h5",
                "print_progress": True
            },
            "domains": [
                {
                    "name": "icr_cell",
                    "instrument": "FT-ICR",
                    "geometry": {
                        "origin_m": [0.0, 0.0, 0.0],
                        "length_m": 0.1,
                        "radius_m": 0.02
                    },
                    "env": {
                        "pressure_Pa": 1e-9,
                        "temperature_K": 300.0,
                        "gas_species": "He"
                    },
                    "fields": {
                        "magnetic": {
                            "B_field_T": [0.0, 0.0, 7.0]
                        }
                    }
                }
            ]
        }
    }
}


def list_templates():
    """Print available templates"""
    print("\nAvailable Templates:")
    print("-" * 60)
    for name, tmpl in TEMPLATES.items():
        print(f"  {name:12s} - {tmpl['description']}")
    print("-" * 60)


def load_template(template_name: str) -> Dict[str, Any]:
    """Load a template without modifications"""
    if template_name not in TEMPLATES:
        print(f"Error: Unknown template '{template_name}'")
        list_templates()
        sys.exit(1)
    
    # Deep copy to avoid modifying the original
    return json.loads(json.dumps(TEMPLATES[template_name]["config"]))


def validate_config(config: Dict[str, Any]) -> List[str]:
    """Basic validation of configuration structure"""
    errors = []
    
    # Check required sections
    required_sections = ["simulation", "physics", "output", "domains"]
    for section in required_sections:
        if section not in config:
            errors.append(f"Missing required section: '{section}'")
    
    # Check simulation required fields
    if "simulation" in config:
        if "total_time_s" not in config["simulation"]:
            errors.append("Missing 'total_time_s' in simulation section")
        if "dt_s" not in config["simulation"]:
            errors.append("Missing 'dt_s' in simulation section")
    
    # Check domains
    if "domains" in config:
        if not isinstance(config["domains"], list) or len(config["domains"]) == 0:
            errors.append("'domains' must be a non-empty array")
        else:
            for i, domain in enumerate(config["domains"]):
                if "name" not in domain:
                    errors.append(f"Domain {i}: missing 'name' field")
                if "instrument" not in domain:
                    errors.append(f"Domain {i}: missing 'instrument' field")
                if "geometry" not in domain:
                    errors.append(f"Domain {i}: missing 'geometry' section")
                if "env" not in domain:
                    errors.append(f"Domain {i}: missing 'env' section")
    
    return errors


def main():
    parser = argparse.ArgumentParser(
        description="ICARION Configuration Builder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                    # Interactive mode
  %(prog)s --template ims                     # Use IMS template
  %(prog)s --template tof --time 1e-4         # Customize TOF template
  %(prog)s --list-templates                   # Show available templates
  %(prog)s --template ims --output my.json    # Save to file
        """
    )
    
    parser.add_argument("--template", "-t", help="Use a template (minimal/ims/tof/lqit/orbitrap)")
    parser.add_argument("--list-templates", "-l", action="store_true", help="List available templates")
    parser.add_argument("--output", "-o", default="config.json", help="Output file (default: config.json)")
    
    args = parser.parse_args()
    
    # List templates
    if args.list_templates:
        list_templates()
        sys.exit(0)
    
    # Generate config
    if args.template:
        config = load_template(args.template)
        print(f"\n✓ Using template: {args.template}")
    else:
        parser.print_help()
        sys.exit(1)
    
    # Basic validation
    errors = validate_config(config)
    if errors:
        print("\n⚠ Configuration Validation Warnings:")
        for error in errors:
            print(f"  - {error}")
    
    # Write to file
    output_path = Path(args.output)
    with open(output_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"\n✓ Configuration written to: {output_path.absolute()}")
    print(f"\nValidate with:")
    print(f"  python3 schema/validator.py {output_path}")
    print(f"\nRun simulation:")
    print(f"  ./build/src/icarion_main {output_path}")
    print(f"\nEdit manually if needed, then validate before running.")


if __name__ == "__main__":
    main()
