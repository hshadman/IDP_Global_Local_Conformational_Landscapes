import os, sys
# Add repo/src to sys.path so imports work without needing pip install -e .
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
