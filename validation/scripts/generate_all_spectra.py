#!/usr/bin/env python3
"""
Master script to generate all publication-quality instrument spectra
Creates scientific visualizations for TOF, Orbitrap, and FTICR instruments
"""

import subprocess
import sys
from pathlib import Path
import time

def run_script(script_path, description):
    """Run a Python script and handle errors"""
    print(f"\n{'='*60}")
    print(f"Generating {description}...")
    print(f"{'='*60}")
    
    try:
        start_time = time.time()
        result = subprocess.run([sys.executable, script_path], 
                              capture_output=True, text=True, cwd=script_path.parent)
        end_time = time.time()
        
        if result.returncode == 0:
            print(f"✅ {description} completed successfully!")
            print(f"   Execution time: {end_time - start_time:.2f} seconds")
            if result.stdout:
                print(f"   Output: {result.stdout.strip()}")
        else:
            print(f"❌ {description} failed!")
            print(f"   Error: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"❌ Failed to run {description}: {e}")
        return False
    
    return True

def main():
    """Generate all instrument spectra"""
    print("ICARION Validation: Publication-Quality Spectrum Generation")
    print("=" * 60)
    print("Creating scientific visualizations for mass spectrometry instruments")
    print("Replacing tabular presentations with spectral analysis plots")
    
    scripts_dir = Path('/home/chsch95/ICARION/validation/scripts')
    
    # Define scripts to run
    spectrum_scripts = [
        {
            'script': scripts_dir / 'create_tof_spectrum.py',
            'description': 'TOF Mass Spectrum (Time-of-Flight Analysis)'
        },
        {
            'script': scripts_dir / 'create_orbitrap_spectrum.py', 
            'description': 'Orbitrap FFT Spectrum (Frequency Domain Analysis)'
        },
        {
            'script': scripts_dir / 'create_fticr_spectrum.py',
            'description': 'FTICR Cyclotron Spectrum (Ultra-High Resolution)'
        }
    ]
    
    # Track success/failure
    results = []
    total_start = time.time()
    
    # Run each spectrum generation script
    for spectrum in spectrum_scripts:
        if spectrum['script'].exists():
            success = run_script(spectrum['script'], spectrum['description'])
            results.append({
                'name': spectrum['description'],
                'success': success
            })
        else:
            print(f"⚠️  Script not found: {spectrum['script']}")
            results.append({
                'name': spectrum['description'],
                'success': False
            })
    
    total_end = time.time()
    
    # Summary report
    print(f"\n{'='*60}")
    print("SPECTRUM GENERATION SUMMARY")
    print(f"{'='*60}")
    print(f"Total execution time: {total_end - total_start:.2f} seconds")
    print(f"")
    
    successful = sum(1 for r in results if r['success'])
    total = len(results)
    
    for result in results:
        status = "✅ SUCCESS" if result['success'] else "❌ FAILED"
        print(f"{status}: {result['name']}")
    
    print(f"\nOverall Status: {successful}/{total} spectra generated successfully")
    
    if successful == total:
        print("\n🎉 All spectra generated! Publication-quality plots are ready.")
        print("\nGenerated files:")
        print("  📊 validation/figures/tof_mass_spectrum_validation.png")
        print("  📊 validation/figures/orbitrap_fft_spectrum_validation.png") 
        print("  📊 validation/figures/fticr_cyclotron_spectrum_validation.png")
        print("\nDetailed analysis logs:")
        print("  📋 validation/logs/TOF_SPECTRUM_ANALYSIS.txt")
        print("  📋 validation/logs/ORBITRAP_FFT_ANALYSIS.txt")
        print("  📋 validation/logs/FTICR_CYCLOTRON_ANALYSIS.txt")
        
        print(f"\n{'='*60}")
        print("VISUALIZATION QUALITY UPGRADE COMPLETE")
        print(f"{'='*60}")
        print("• TOF: Tabular → Mass Spectrum Analysis")
        print("• Orbitrap: Simple plot → FFT Frequency Spectrum")  
        print("• FTICR: Basic validation → Cyclotron Spectrum")
        print("• All plots now match the quality of thermalization/IMS/quadrupole")
        print("• Ready for publication in validation documentation")
        
    else:
        failed_count = total - successful
        print(f"\n⚠️  {failed_count} spectrum generation(s) failed. Check error messages above.")
        
    return successful == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)