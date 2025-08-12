#!/usr/bin/env python
"""
base_tab.py

Base class for all tab components in the Spacecraft Constellation Fault Simulator.
Provides common functionality for all tabs.
"""
import tkinter as tk
from tkinter import ttk, scrolledtext
import sys
import os

# Add parent directory to path to import help content
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from rw_fault_help import  get_general_help
except ImportError:
    print("Warning: Could not import help content module")
    get_general_help = lambda section: f"No help content available for: {section}"

class BaseTab:
    """Base class for all tab components"""
    
    def __init__(self, parent_app, parent_frame):
        """
        Initialize the base tab
        
        Parameters:
        parent_app (SatelliteSimulatorApp): The parent application instance
        parent_frame (ttk.Frame): The parent frame to build the tab in
        """
        self.parent_app = parent_app
        self.parent_frame = parent_frame
        
    def add_help_button(self, parent, title, topic=None, command=None):
        """
        Add a help button to a frame with a popup dialog for the given topic
        
        Parameters:
        parent (ttk.Frame): Parent frame to add the button to
        title (str): Title for the help dialog
        topic (str, optional): Topic identifier for looking up help content
        command (callable, optional): Custom command to run, overrides default help dialog
        
        Returns:
        ttk.Button: The created help button
        """
        if command is None:
            command = lambda: self.show_help_content(title, topic)
            
        try:
            help_button = ttk.Button(parent, text="?", width=2, style='Help.TButton',
                                   command=command)
            help_button.pack(side=tk.RIGHT, padx=5, pady=5, anchor='ne')
            return help_button
        except:
            # Fallback if style doesn't exist
            help_button = ttk.Button(parent, text="?", width=2, command=command)
            help_button.pack(side=tk.RIGHT, padx=5, pady=5, anchor='ne')
            return help_button
        
    def show_help_content(self, title, topic=None, content=None):
        """
        Show help content in a popup dialog
        
        Parameters:
        title (str): Title for the help dialog
        topic (str, optional): Topic identifier for looking up help content
        content (str, optional): Explicit content to show, overrides topic lookup
        """
        try:
            help_window = tk.Toplevel(self.parent_app.root)
            help_window.title(f"Help: {title}")
            help_window.geometry("600x400")
            help_window.transient(self.parent_app.root)  # Make window modal
            help_window.grab_set()
            
            # Create a frame for the help content
            frame = ttk.Frame(help_window, padding=10)
            frame.pack(fill=tk.BOTH, expand=True)
            
            # Add title
            try:
                ttk.Label(frame, text=title, style='Title.TLabel').pack(fill=tk.X, pady=(0, 10))
            except:
                # Fallback if style doesn't exist
                ttk.Label(frame, text=title, font=('Segoe UI', 12, 'bold')).pack(fill=tk.X, pady=(0, 10))
            
            # Create scrolled text widget for content
            help_text = scrolledtext.ScrolledText(frame, wrap=tk.WORD, width=70, height=20)
            help_text.pack(fill=tk.BOTH, expand=True, pady=5)
            
            # Get the content
            if content is None and topic:
                # Try to get help content directly from the help module
                try:
                    # Special topics that use get_general_help
                    if topic in ["overview", "simulation", "visualization"]:
                        content = get_general_help(topic)
                    else:
                        content = get_help_content(topic)
                except Exception as e:
                    print(f"Error getting help content for {topic}: {e}")
                    content = f"No help content available for topic: {topic}"
            elif content is None:
                content = f"No help content available for: {title}"
                    
            help_text.insert(tk.END, content)
            
            # Make it read-only
            help_text.config(state=tk.DISABLED)
            
            # Add close button
            ttk.Button(frame, text="Close", command=help_window.destroy).pack(pady=10)
        except Exception as e:
            print(f"Error showing help content: {e}")
        
    def show_message(self, title, message, message_type="info"):
        """
        Show a message dialog
        
        Parameters:
        title (str): Dialog title
        message (str): Message to display
        message_type (str): Type of message - 'info', 'warning', or 'error'
        """
        try:
            from tkinter import messagebox
            
            if message_type == "info":
                messagebox.showinfo(title, message)
            elif message_type == "warning":
                messagebox.showwarning(title, message)
            elif message_type == "error":
                messagebox.showerror(title, message)
            else:
                messagebox.showinfo(title, message)
        except Exception as e:
            print(f"Error showing message: {e}")
            
    def add_log(self, message):
        """
        Add a message to the application log
        
        Parameters:
        message (str): Message to add to the log
        """
        try:
            if hasattr(self.parent_app, 'add_log'):
                self.parent_app.add_log(message)
            else:
                print(f"Log: {message}")  # Fallback to console
        except Exception as e:
            print(f"Log: {message}")  # Fallback to console if error