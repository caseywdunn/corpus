#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path

def run_snakemake_with_fallback():
    """Run snakemake with fallback options if CBC solver is not available"""
    
    # Try different scheduler options
    scheduler_options = [
        [],  # Default scheduler
        ["--scheduler", "greedy"],  # Greedy scheduler
        ["--scheduler", "ilp", "--scheduler-solver-path", ""],  # Try to skip solver
    ]
    
    base_cmd = ["snakemake", "--cores", "1"]
    
    for scheduler_opts in scheduler_options:
        cmd = base_cmd + scheduler_opts
        
        print(f"Trying: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print("SUCCESS!")
                print(result.stdout)
                return True
            else:
                print(f"Failed with return code {result.returncode}")
                if "PulpSolverError" not in result.stderr:
                    print("STDERR:", result.stderr)
                    return False
                
        except Exception as e:
            print(f"Exception: {e}")
            continue
    
    print("\nAll scheduler options failed. You may need to install CBC solver:")
    print("conda install -c conda-forge coincbc")
    print("or")
    print("pip install coincbc")
    return False

if __name__ == "__main__":
    success = run_snakemake_with_fallback()
    sys.exit(0 if success else 1)