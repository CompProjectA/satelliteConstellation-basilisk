#!/usr/bin/env python
"""Test DRL Tab"""
import sys
import os
import tkinter as tk
from tkinter import ttk

# Add paths
sys.path.insert(0, os.path.dirname(__file__))

print("=" * 60)
print("DRL TAB DIAGNOSTIC")
print("=" * 60)

# Test 1: Import
print("\n1. Testing import...")
try:
    from gui_tab.drl_tab import DRLTab
    print("   ✓ DRLTab imported successfully")
except Exception as e:
    print(f"   ✗ Import failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 2: Check class
print("\n2. Checking class structure...")
print(f"   Class: {DRLTab}")
print(f"   Has __init__: {hasattr(DRLTab, '__init__')}")
print(f"   Module: {DRLTab.__module__}")

# Test 3: Create test window
print("\n3. Creating test window...")
try:
    root = tk.Tk()
    root.title("DRL Tab Test")
    root.geometry("800x600")
    
    print("   ✓ Root window created")
    
    # Create DRL tab
    print("\n4. Creating DRL Tab...")
    drl_tab = DRLTab(root, parent_app=None, default_script="", default_results_dir="")
    drl_tab.pack(fill="both", expand=True)
    print("   ✓ DRL Tab created and packed")
    
    # Check children
    print("\n5. Checking tab structure...")
    children = drl_tab.winfo_children()
    print(f"   Children count: {len(children)}")
    for i, child in enumerate(children):
        print(f"   - Child {i}: {type(child).__name__}")
    
    if len(children) > 0:
        print("\n   ✓ DRL Tab has content!")
        # Check notebook
        for child in children:
            if isinstance(child, ttk.Notebook):
                tabs = child.tabs()
                print(f"\n   Found notebook with {len(tabs)} tabs:")
                for tab in tabs:
                    print(f"   - {child.tab(tab, 'text')}")
    else:
        print("\n   ✗ DRL Tab is EMPTY - no children widgets!")
    
    print("\n6. Showing window for 5 seconds...")
    print("   (Check if you can see the DRL tab content)")
    root.after(5000, root.destroy)
    root.mainloop()
    
    print("\n✓ TEST COMPLETED SUCCESSFULLY")
    
except Exception as e:
    print(f"\n✗ Error creating tab: {e}")
    import traceback
    traceback.print_exc()

print("=" * 60)