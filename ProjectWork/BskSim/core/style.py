#!/usr/bin/env python
"""
style.py

Enhanced styling module for the Spacecraft Fault Simulator GUI
Integrates styling features from the original design with modern enhancements
"""
import tkinter as tk
from tkinter import ttk

def setup_style(root):
    """Set up enhanced styling for the application"""
    
    # Create enhanced style object
    style = ttk.Style()
    
    # Set enhanced theme
    try:
        # Try to use a modern theme if available
        available_themes = style.theme_names()
        if 'clam' in available_themes:
            style.theme_use('clam')
        elif 'alt' in available_themes:
            style.theme_use('alt')
        else:
            style.theme_use('default')
    except:
        style.theme_use('default')
    
    # Enhanced color scheme (from old file inspiration)
    colors = {
        'bg_primary': '#f5f5f5',
        'bg_secondary': '#ffffff',
        'accent_blue': '#2E86C1',
        'accent_green': '#28B463',
        'accent_red': '#E74C3C',
        'accent_orange': '#F39C12',
        'text_primary': '#2C3E50',
        'text_secondary': '#7F8C8D',
        'border': '#BDC3C7'
    }
    
    # Configure enhanced styles
    
    # Enhanced Title style
    style.configure('Title.TLabel',
                   font=('Segoe UI', 16, 'bold'),
                   foreground=colors['text_primary'],
                   background=colors['bg_primary'])
    
    # Enhanced Run button style (prominent like old file)
    style.configure('Run.TButton',
                   font=('Segoe UI', 10, 'bold'),
                   foreground='white',
                   background=colors['accent_green'],
                   borderwidth=2,
                   focuscolor='none',
                   padding=(20, 10))
    
    # Enhanced Run button hover effect
    style.map('Run.TButton',
             background=[('active', '#239B56'),
                        ('pressed', '#1E8449')])
    
    # Enhanced Help button style (from old file)
    style.configure('Help.TButton',
                   font=('Segoe UI', 8),
                   foreground=colors['text_primary'],
                   background=colors['bg_secondary'],
                   borderwidth=1,
                   width=2)
    
    # Enhanced Info label style
    style.configure('Info.TLabel',
                   font=('Segoe UI', 9, 'italic'),
                   foreground=colors['text_secondary'],
                   background=colors['bg_primary'])
    
    # Enhanced Status frame and label
    style.configure('Status.TFrame',
                   background=colors['bg_secondary'],
                   relief='sunken',
                   borderwidth=1)
    
    style.configure('Status.TLabel',
                   font=('Segoe UI', 9),
                   foreground=colors['text_primary'],
                   background=colors['bg_secondary'],
                   padding=(5, 2))
    
    # Enhanced Notebook styling
    style.configure('TNotebook',
                   background=colors['bg_primary'],
                   borderwidth=1,
                   tabmargins=[2, 5, 2, 0])
    
    style.configure('TNotebook.Tab',
                   font=('Segoe UI', 9),
                   foreground=colors['text_primary'],
                   background=colors['bg_secondary'],
                   padding=[20, 8],
                   borderwidth=1)
    
    style.map('TNotebook.Tab',
             background=[('selected', colors['accent_blue']),
                        ('active', '#AED6F1')],
             foreground=[('selected', 'white'),
                        ('active', colors['text_primary'])])
    
    # Enhanced LabelFrame styling
    style.configure('TLabelframe',
                   background=colors['bg_primary'],
                   borderwidth=2,
                   relief='groove')
    
    style.configure('TLabelframe.Label',
                   font=('Segoe UI', 10, 'bold'),
                   foreground=colors['text_primary'],
                   background=colors['bg_primary'])
    
    # Enhanced Button styling
    style.configure('TButton',
                   font=('Segoe UI', 9),
                   foreground=colors['text_primary'],
                   background=colors['bg_secondary'],
                   borderwidth=1,
                   focuscolor='none',
                   padding=(10, 5))
    
    style.map('TButton',
             background=[('active', colors['accent_blue']),
                        ('pressed', '#2874A6')],
             foreground=[('active', 'white'),
                        ('pressed', 'white')])
    
    # Enhanced Entry styling
    style.configure('TEntry',
                   font=('Segoe UI', 9),
                   foreground=colors['text_primary'],
                   borderwidth=1,
                   insertcolor=colors['text_primary'])
    
    # Enhanced Combobox styling
    style.configure('TCombobox',
                   font=('Segoe UI', 9),
                   foreground=colors['text_primary'],
                   borderwidth=1)
    
    # Enhanced Checkbutton styling
    style.configure('TCheckbutton',
                   font=('Segoe UI', 9),
                   foreground=colors['text_primary'],
                   background=colors['bg_primary'],
                   focuscolor='none')
    
    # Enhanced Radiobutton styling
    style.configure('TRadiobutton',
                   font=('Segoe UI', 9),
                   foreground=colors['text_primary'],
                   background=colors['bg_primary'],
                   focuscolor='none')
    
    # Enhanced Spinbox styling
    style.configure('TSpinbox',
                   font=('Segoe UI', 9),
                   foreground=colors['text_primary'],
                   borderwidth=1)
    
    # Enhanced Treeview styling
    style.configure('Treeview',
                   font=('Segoe UI', 9),
                   foreground=colors['text_primary'],
                   background=colors['bg_secondary'],
                   borderwidth=1,
                   rowheight=25)
    
    style.configure('Treeview.Heading',
                   font=('Segoe UI', 9, 'bold'),
                   foreground=colors['text_primary'],
                   background=colors['bg_primary'],
                   borderwidth=1)
    
    style.map('Treeview',
             background=[('selected', colors['accent_blue'])],
             foreground=[('selected', 'white')])
    
    # Enhanced Scale styling
    style.configure('TScale',
                   background=colors['bg_primary'],
                   borderwidth=1,
                   lightcolor=colors['accent_blue'],
                   darkcolor=colors['accent_blue'])
    
    # Enhanced Progressbar styling
    style.configure('TProgressbar',
                   background=colors['accent_green'],
                   borderwidth=1,
                   lightcolor=colors['accent_green'],
                   darkcolor=colors['accent_green'])
    
    # Enhanced Separator styling
    style.configure('TSeparator',
                   background=colors['border'])
    
    # Enhanced Scrollbar styling
    style.configure('TScrollbar',
                   background=colors['bg_secondary'],
                   borderwidth=1,
                   arrowcolor=colors['text_primary'])
    
    # Set enhanced root window properties
    root.configure(bg=colors['bg_primary'])
    
    # Enhanced window icon (if available)
    try:
        import os
        assets_dir = os.path.join(os.path.dirname(__file__), '..', 'assets')
        icon_path = os.path.join(assets_dir, "satellite_icon.ico")
        if os.path.exists(icon_path):
            root.iconbitmap(icon_path)
    except:
        pass  # Skip if icon not available
    
    return style


def create_enhanced_tooltip(widget, text):
    
    def on_enter(event):
        tooltip = tk.Toplevel()
        tooltip.wm_overrideredirect(True)
        tooltip.wm_geometry(f"+{event.x_root+10}+{event.y_root+10}")
        
        label = tk.Label(tooltip, text=text,
                        font=('Segoe UI', 8),
                        background='#FFFFCC',
                        foreground='#000000',
                        relief='solid',
                        borderwidth=1,
                        padx=5, pady=3)
        label.pack()
        
        widget.tooltip = tooltip
    
    def on_leave(event):
        if hasattr(widget, 'tooltip'):
            widget.tooltip.destroy()
            del widget.tooltip
    
    widget.bind('<Enter>', on_enter)
    widget.bind('<Leave>', on_leave)


def apply_enhanced_grid_weights(parent, rows=None, cols=None):
    """Apply enhanced grid weights for better layout (from old file concept)"""
    if rows:
        for row in rows:
            parent.grid_rowconfigure(row, weight=1)
    
    if cols:
        for col in cols:
            parent.grid_columnconfigure(col, weight=1)


# Enhanced color constants for consistency
ENHANCED_COLORS = {
    'RED': '#E74C3C',
    'BLUE': '#3498DB', 
    'GREEN': '#2ECC71',
    'YELLOW': '#F1C40F',
    'PURPLE': '#9B59B6',
    'ORANGE': '#E67E22',
    'CYAN': '#1ABC9C',
    'MAGENTA': '#E91E63',
    'DARK_BLUE': '#2C3E50',
    'LIGHT_GRAY': '#ECF0F1',
    'DARK_GRAY': '#7F8C8D'
}


def get_enhanced_color(color_name):
    """Get enhanced color value by name"""
    return ENHANCED_COLORS.get(color_name.upper(), '#000000')


def create_enhanced_separator(parent, orient='horizontal', **kwargs):
    """Create enhanced separator with better styling"""
    separator = ttk.Separator(parent, orient=orient)
    if kwargs:
        separator.configure(**kwargs)
    return separator


def create_enhanced_labelframe(parent, text, **kwargs):
    """Create enhanced label frame with consistent styling"""
    frame = ttk.LabelFrame(parent, text=text, **kwargs)
    frame.configure(padding=10)
    return frame