#!/usr/bin/env python
"""
style.py

Defines the visual styling for the Spacecraft Fault Simulator GUI.
This module separates styling concerns from the main application code.
"""
import tkinter as tk
from tkinter import ttk

def setup_style(root):
    """
    Set up the application's visual style
    
    Parameters:
    root (tk.Tk): The root window of the application
    
    Returns:
    ttk.Style: The configured style object
    """
    # Create a modern style for the application
    style = ttk.Style(root)
    
    # Try to use a modern theme if available
    available_themes = style.theme_names()
    if 'clam' in available_themes:
        style.theme_use('clam')
    elif 'vista' in available_themes:
        style.theme_use('vista')
    
    # Configure styles for different widgets
    style.configure('TFrame', background='#f5f5f5')
    style.configure('TLabel', background='#f5f5f5', font=('Segoe UI', 10))
    style.configure('TButton', font=('Segoe UI', 10))
    style.configure('Header.TLabel', font=('Segoe UI', 12, 'bold'))
    style.configure('Title.TLabel', font=('Segoe UI', 14, 'bold'), padding=5)
    style.configure('Section.TLabel', font=('Segoe UI', 11, 'bold'), padding=3)
    style.configure('Info.TLabel', font=('Segoe UI', 9), foreground='#555555')
    
    # Configure Notebook tab style
    style.configure('TNotebook.Tab', padding=(12, 8), font=('Segoe UI', 10, 'bold'))
    style.map('TNotebook.Tab', background=[('selected', '#4a6ea9'), ('!selected', '#e1e1e1')],
                foreground=[('selected', '#ffffff'), ('!selected', '#333333')])
    
    # Configure Treeview
    style.configure('Treeview.Heading', font=('Segoe UI', 10, 'bold'))
    style.configure('Treeview', rowheight=25)

    # Status bar style
    style.configure('Status.TFrame', background='#e1e1e1')
    style.configure('Status.TLabel', background='#e1e1e1', font=('Segoe UI', 10))
    
    # Help button style
    style.configure('Help.TButton', font=('Segoe UI', 9), padding=2)
    
    # Run button style - make it MUCH more prominent
    style.configure('Run.TButton', 
                   font=('Segoe UI', 11, 'bold'),
                   padding=(10, 5),
                   background='#2d8659',  # Green background
                   foreground='white')    # White text
                   
    style.map('Run.TButton', 
              background=[('active', '#1d6849')],  # Darker green when active
              foreground=[('active', 'white')])
              
    return style