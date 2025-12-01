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
    }
}


def list_templates():
    """Print available templates"""
    print("\nAvailable Templates:")
    print("-" * 60)
    for name, tmpl in TEMPLATES.items():
        print(f"  {name:12s} - {tmpl['description']}")
    print("-" * 60)


def get_user_input(prompt: str, default: Any = None, value_type: type = str) -> Any:
    """Get user input with default value"""
    if default is not None:
        prompt_str = f"{prompt} [{default}]: "
    else:
        prompt_str = f"{prompt}: "
    
    user_input = input(prompt_str).strip()
    
    if not user_input and default is not None:
        return default
    
    if value_type == bool:
        return user_input.lower() in ('yes', 'y', 'true', '1')
    elif value_type == float:
        return float(user_input)
    elif value_type == int:
        return int(user_input)
    else:
        return user_input


def interactive_mode() -> Dict[str, Any]:
    """Interactive configuration builder"""
    print("\n" + "=" * 60)
    print("ICARION Configuration Builder - Interactive Mode")
    print("=" * 60 + "\n")
    
    # Choose base template
    print("Available templates:")
    for name, tmpl in TEMPLATES.items():
        print(f"  - {name}: {tmpl['description']}")
    
    template_name = get_user_input("\nChoose template (or 'custom')", "minimal")
    
    if template_name in TEMPLATES:
        config = json.loads(json.dumps(TEMPLATES[template_name]["config"]))  # Deep copy
        print(f"\n✓ Loaded template: {template_name}")
    else:
        config = json.loads(json.dumps(TEMPLATES["minimal"]["config"]))
        print(f"\n✓ Starting with minimal template")
    
    # Simulation parameters
    print("\n--- Simulation Parameters ---")
    config["simulation"]["total_time_s"] = get_user_input(
        "Total simulation time [s]", 
        config["simulation"]["total_time_s"], 
        float
    )
    config["simulation"]["dt_s"] = get_user_input(
        "Time step [s]", 
        config["simulation"]["dt_s"], 
        float
    )
    config["simulation"]["integrator"] = get_user_input(
        "Integrator (RK4/RK45/Boris)", 
        config["simulation"]["integrator"]
    )
    config["simulation"]["enable_gpu"] = get_user_input(
        "Enable GPU (yes/no)", 
        config["simulation"].get("enable_gpu", False), 
        bool
    )
    
    # Physics
    print("\n--- Physics Settings ---")
    print("Collision models: NoCollisions, HSD, HSS, EHSS, SDS")
    config["physics"]["collision_model"] = get_user_input(
        "Collision model", 
        config["physics"]["collision_model"]
    )
    
    # Output
    print("\n--- Output Configuration ---")
    config["output"]["folder"] = get_user_input(
        "Output folder", 
        config["output"]["folder"]
    )
    config["output"]["trajectory_file"] = get_user_input(
        "Trajectory file name", 
        config["output"]["trajectory_file"]
    )
    
    # Domains
    print("\n--- Domain Configuration ---")
    num_domains = get_user_input(
        f"Number of domains", 
        len(config["domains"]), 
        int
    )
    
    if num_domains != len(config["domains"]):
        config["domains"] = []
        for i in range(num_domains):
            print(f"\n  Domain {i+1}:")
            domain = {
                "name": get_user_input("    Name", f"domain_{i+1}"),
                "instrument": get_user_input("    Instrument type (IMS/TOF/LQIT/Orbitrap)", "IMS"),
                "geometry": {
                    "origin_m": [0.0, 0.0, float(get_user_input("    Z-origin [m]", 0.0, float))],
                    "length_m": get_user_input("    Length [m]", 0.1, float),
                    "radius_m": get_user_input("    Radius [m]", 0.01, float)
                },
                "env": {
                    "pressure_Pa": get_user_input("    Pressure [Pa]", 101325.0, float),
                    "temperature_K": get_user_input("    Temperature [K]", 300.0, float),
                    "gas_species": get_user_input("    Gas species (He/N2/Ar)", "He")
                }
            }
            config["domains"].append(domain)
    
    return config


def customize_template(template_name: str, args: argparse.Namespace) -> Dict[str, Any]:
    """Customize a template based on command-line args"""
    if template_name not in TEMPLATES:
        print(f"Error: Unknown template '{template_name}'")
        list_templates()
        sys.exit(1)
    
    config = json.loads(json.dumps(TEMPLATES[template_name]["config"]))  # Deep copy
    
    # Apply customizations
    if args.time is not None:
        config["simulation"]["total_time_s"] = args.time
    if args.timestep is not None:
        config["simulation"]["dt_s"] = args.timestep
    if args.integrator is not None:
        config["simulation"]["integrator"] = args.integrator
    if args.gpu is not None:
        config["simulation"]["enable_gpu"] = args.gpu
    
    if args.collision_model is not None:
        config["physics"]["collision_model"] = args.collision_model
    
    if args.output_folder is not None:
        config["output"]["folder"] = args.output_folder
    
    # Customize first domain
    if len(config["domains"]) > 0:
        if args.pressure is not None:
            config["domains"][0]["env"]["pressure_Pa"] = args.pressure
        if args.temperature is not None:
            config["domains"][0]["env"]["temperature_K"] = args.temperature
        if args.gas is not None:
            config["domains"][0]["env"]["gas_species"] = args.gas
    
    return config


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
    parser.add_argument("--interactive", "-i", action="store_true", help="Force interactive mode")
    
    # Customization options
    parser.add_argument("--time", type=float, help="Total simulation time [s]")
    parser.add_argument("--timestep", type=float, help="Time step [s]")
    parser.add_argument("--integrator", choices=["RK4", "RK45", "Boris"], help="Integrator type")
    parser.add_argument("--gpu", type=bool, help="Enable GPU")
    parser.add_argument("--collision-model", help="Collision model (NoCollisions/HSD/HSS/EHSS)")
    parser.add_argument("--pressure", type=float, help="Gas pressure [Pa]")
    parser.add_argument("--temperature", type=float, help="Gas temperature [K]")
    parser.add_argument("--gas", help="Gas species (He/N2/Ar)")
    parser.add_argument("--output-folder", help="Output folder path")
    
    args = parser.parse_args()
    
    # List templates
    if args.list_templates:
        list_templates()
        sys.exit(0)
    
    # Generate config
    if args.interactive or (args.template is None and not args.list_templates):
        config = interactive_mode()
    elif args.template:
        config = customize_template(args.template, args)
        print(f"\n✓ Using template: {args.template}")
    else:
        parser.print_help()
        sys.exit(1)
    
    # Validate
    errors = validate_config(config)
    if errors:
        print("\n⚠ Configuration Validation Warnings:")
        for error in errors:
            print(f"  - {error}")
        
        if not args.interactive:
            proceed = input("\nProceed anyway? (yes/no): ").strip().lower()
            if proceed not in ('yes', 'y'):
                print("Aborted.")
                sys.exit(1)
    
    # Write to file
    output_path = Path(args.output)
    with open(output_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"\n✓ Configuration written to: {output_path.absolute()}")
    print(f"\nValidate with:")
    print(f"  python3 src/core/config/schema/validator.py {output_path}")
    print(f"  ./build/src/icarion_main --validate-config {output_path}")
    print(f"\nRun simulation:")
    print(f"  ./build/src/icarion_main {output_path}")


if __name__ == "__main__":
    main()
