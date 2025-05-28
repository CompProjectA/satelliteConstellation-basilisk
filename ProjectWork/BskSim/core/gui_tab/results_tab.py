import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from datetime import datetime
import matplotlib.pyplot as plt

class ResultsTab:
    def __init__(self, app, parent):
        self.app = app
        self.parent = parent
        self.plot_figure = None
        self.plot_canvas = None
        self.plot_toolbar = None
        self.plot_listbox = None
        self.plot_info_text = None
        self.current_plot_path = None

        self.create_ui()

    def create_ui(self):
        main_frame = ttk.Frame(self.parent, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        controls_frame = ttk.Frame(main_frame)
        controls_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(controls_frame, text="Simulation Results", style="Title.TLabel").pack(side=tk.LEFT)

        ttk.Button(controls_frame, text="Refresh", command=self.refresh_plot_list).pack(side=tk.RIGHT, padx=5)
        ttk.Button(controls_frame, text="Open Results Folder", command=self.app.open_results_folder).pack(side=tk.RIGHT, padx=5)
        ttk.Button(controls_frame, text="Export Plot", command=self.export_current_plot).pack(side=tk.RIGHT, padx=5)

        split_frame = ttk.Frame(main_frame)
        split_frame.pack(fill=tk.BOTH, expand=True)

        list_frame = ttk.LabelFrame(split_frame, text="Available Plots", padding=10, width=320)
        list_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))

        list_container = ttk.Frame(list_frame)
        list_container.pack(fill=tk.BOTH, expand=True)

        scrollbar = ttk.Scrollbar(list_container)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.plot_listbox = tk.Listbox(list_container, selectmode=tk.SINGLE, yscrollcommand=scrollbar.set, font=('Segoe UI', 10), exportselection=False)
        self.plot_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.plot_listbox.yview)

        self.plot_listbox.bind('<<ListboxSelect>>', self.on_plot_selected)

        display_frame = ttk.LabelFrame(split_frame, text="Plot View", padding=10)
        display_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0))

        self.plot_figure = Figure(figsize=(12, 8), dpi=100)
        self.plot_canvas = FigureCanvasTkAgg(self.plot_figure, display_frame)
        self.plot_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar_frame = ttk.Frame(display_frame)
        toolbar_frame.pack(fill=tk.X, pady=(5, 0))
        self.plot_toolbar = NavigationToolbar2Tk(self.plot_canvas, toolbar_frame)
        self.plot_toolbar.update()

        info_frame = ttk.LabelFrame(main_frame, text="Plot Information", padding=10)
        info_frame.pack(fill=tk.X, pady=(10, 0))

        self.plot_info_text = scrolledtext.ScrolledText(info_frame, height=4, wrap=tk.WORD)
        self.plot_info_text.pack(fill=tk.BOTH, expand=True)
        self.plot_info_text.config(state=tk.DISABLED)

        self.refresh_plot_list()

    def refresh_plot_list(self):
        self.plot_listbox.delete(0, tk.END)

        self.plot_info_text.config(state=tk.NORMAL)
        self.plot_info_text.delete(1.0, tk.END)
        self.plot_info_text.insert(tk.END, "Select a plot to view details and analysis")
        self.plot_info_text.config(state=tk.DISABLED)

        self.plot_figure.clear()
        ax = self.plot_figure.add_subplot(111)
        ax.text(0.5, 0.5, "No plot selected\nSelect a plot from the list to view", ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.axis('off')
        self.plot_canvas.draw()

        if not os.path.exists(self.app.plots_dir):
            return

        plot_files = [f for f in os.listdir(self.app.plots_dir) if f.endswith(('.png', '.jpg', '.jpeg', '.pdf', '.svg'))]
        plot_files.sort(key=lambda x: os.path.getmtime(os.path.join(self.app.plots_dir, x)), reverse=True)

        for filename in plot_files:
            self.plot_listbox.insert(tk.END, filename)

        if plot_files:
            self.plot_listbox.selection_set(0)
            self.on_plot_selected(None)

    def on_plot_selected(self, event):
        selection = self.plot_listbox.curselection()
        if not selection:
            return

        index = selection[0]
        filename = self.plot_listbox.get(index)
        full_path = os.path.join(self.app.plots_dir, filename)
        self.current_plot_path = full_path

        self.plot_info_text.config(state=tk.NORMAL)
        self.plot_info_text.delete(1.0, tk.END)

        try:
            file_size = os.path.getsize(full_path) / 1024
            mod_time = datetime.fromtimestamp(os.path.getmtime(full_path))
            info_text = f"File: {filename}\nSize: {file_size:.1f} KB\nCreated: {mod_time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            self.plot_info_text.insert(tk.END, info_text)
        except Exception as e:
            self.plot_info_text.insert(tk.END, f"Error reading file info: {e}")
        self.plot_info_text.config(state=tk.DISABLED)

        try:
            self.plot_figure.clear()
            if filename.endswith(('.png', '.jpg', '.jpeg')):
                img = plt.imread(full_path)
                ax = self.plot_figure.add_subplot(111)
                ax.imshow(img)
                ax.axis('off')
                ax.set_title(f"Plot: {filename}", fontsize=14, pad=20)
            else:
                ax = self.plot_figure.add_subplot(111)
                ax.text(0.5, 0.5, f"Plot: {filename}\nUse 'Export Plot' to save", ha='center', va='center', transform=ax.transAxes, fontsize=12)
                ax.axis('off')
            self.plot_canvas.draw()
        except Exception as e:
            self.plot_figure.clear()
            ax = self.plot_figure.add_subplot(111)
            ax.text(0.5, 0.5, f"Error displaying plot:\n{e}", ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.axis('off')
            self.plot_canvas.draw()

    def export_current_plot(self):
        if not self.current_plot_path:
            messagebox.showwarning("No Plot Selected", "Please select a plot first.")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("JPEG", "*.jpg"), ("PDF", "*.pdf"), ("SVG", "*.svg")],
            title="Export Plot"
        )
        if file_path:
            try:
                self.plot_figure.savefig(file_path, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Export Successful", f"Plot exported to:\n{file_path}")
            except Exception as e:
                messagebox.showerror("Export Failed", f"Could not export plot:\n{e}")
