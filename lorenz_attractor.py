#!/usr/bin/env python3
"""
Lorenz Attractor 3D Visualization with GUI
Allows users to modify initial conditions and parameters through input fields
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
import tkinter as tk
from tkinter import ttk


class LorenzAttractor:
    """Class to compute and visualize the Lorenz Attractor"""
    
    def __init__(self, sigma=10.0, rho=28.0, beta=8.0/3.0):
        """
        Initialize Lorenz Attractor with parameters
        
        Parameters:
        -----------
        sigma : float
            Prandtl number (default: 10.0)
        rho : float
            Rayleigh number (default: 28.0)
        beta : float
            Geometric factor (default: 8.0/3.0)
        """
        self.sigma = sigma
        self.rho = rho
        self.beta = beta
        
    def lorenz_system(self, t, state):
        """
        Lorenz system differential equations
        
        Parameters:
        -----------
        t : float
            Time (not used but required by solve_ivp)
        state : array-like
            Current state [x, y, z]
            
        Returns:
        --------
        list
            Derivatives [dx/dt, dy/dt, dz/dt]
        """
        x, y, z = state
        dx = self.sigma * (y - x)
        dy = x * (self.rho - z) - y
        dz = x * y - self.beta * z
        return [dx, dy, dz]
    
    def solve(self, initial_state, t_span=(0, 40), t_points=10000):
        """
        Solve the Lorenz system
        
        Parameters:
        -----------
        initial_state : array-like
            Initial conditions [x0, y0, z0]
        t_span : tuple
            Time span (start, end)
        t_points : int
            Number of time points to evaluate
            
        Returns:
        --------
        tuple
            (t, solution) where solution has shape (t_points, 3)
        """
        t_eval = np.linspace(t_span[0], t_span[1], t_points)
        solution = solve_ivp(
            self.lorenz_system, 
            t_span, 
            initial_state, 
            t_eval=t_eval,
            method='RK45'
        )
        return solution.t, solution.y.T


class LorenzGUI:
    """GUI for Lorenz Attractor visualization with input controls"""
    
    def __init__(self, root):
        """Initialize the GUI"""
        self.root = root
        self.root.title("Lorenz Attractor - 3D Visualization")
        
        # Create main container
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Create control panel
        control_frame = ttk.LabelFrame(main_frame, text="Parameters", padding="10")
        control_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        
        # Initial conditions
        ttk.Label(control_frame, text="Initial Conditions:", font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan=2, pady=5)
        
        ttk.Label(control_frame, text="x0:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.x0_var = tk.StringVar(value="1.0")
        ttk.Entry(control_frame, textvariable=self.x0_var, width=15).grid(row=1, column=1, padx=5, pady=2)
        
        ttk.Label(control_frame, text="y0:").grid(row=2, column=0, sticky=tk.W, padx=5, pady=2)
        self.y0_var = tk.StringVar(value="1.0")
        ttk.Entry(control_frame, textvariable=self.y0_var, width=15).grid(row=2, column=1, padx=5, pady=2)
        
        ttk.Label(control_frame, text="z0:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=2)
        self.z0_var = tk.StringVar(value="1.0")
        ttk.Entry(control_frame, textvariable=self.z0_var, width=15).grid(row=3, column=1, padx=5, pady=2)
        
        # Lorenz parameters
        ttk.Label(control_frame, text="Lorenz Parameters:", font=('Arial', 10, 'bold')).grid(row=4, column=0, columnspan=2, pady=(10, 5))
        
        ttk.Label(control_frame, text="σ (sigma):").grid(row=5, column=0, sticky=tk.W, padx=5, pady=2)
        self.sigma_var = tk.StringVar(value="10.0")
        ttk.Entry(control_frame, textvariable=self.sigma_var, width=15).grid(row=5, column=1, padx=5, pady=2)
        
        ttk.Label(control_frame, text="ρ (rho):").grid(row=6, column=0, sticky=tk.W, padx=5, pady=2)
        self.rho_var = tk.StringVar(value="28.0")
        ttk.Entry(control_frame, textvariable=self.rho_var, width=15).grid(row=6, column=1, padx=5, pady=2)
        
        ttk.Label(control_frame, text="β (beta):").grid(row=7, column=0, sticky=tk.W, padx=5, pady=2)
        self.beta_var = tk.StringVar(value="2.667")
        ttk.Entry(control_frame, textvariable=self.beta_var, width=15).grid(row=7, column=1, padx=5, pady=2)
        
        # Update button
        ttk.Button(control_frame, text="Update Plot", command=self.update_plot).grid(row=8, column=0, columnspan=2, pady=10)
        
        # Create matplotlib figure
        self.fig = plt.Figure(figsize=(10, 8))
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        # Create canvas
        canvas_frame = ttk.Frame(main_frame)
        canvas_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=canvas_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Configure grid weights
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(0, weight=1)
        
        # Initial plot
        self.update_plot()
        
    def update_plot(self):
        """Update the plot with current parameter values"""
        try:
            # Get initial conditions
            x0 = float(self.x0_var.get())
            y0 = float(self.y0_var.get())
            z0 = float(self.z0_var.get())
            
            # Get Lorenz parameters
            sigma = float(self.sigma_var.get())
            rho = float(self.rho_var.get())
            beta = float(self.beta_var.get())
            
            # Create Lorenz attractor
            lorenz = LorenzAttractor(sigma=sigma, rho=rho, beta=beta)
            
            # Solve the system
            t, solution = lorenz.solve([x0, y0, z0])
            
            # Clear previous plot
            self.ax.clear()
            
            # Plot the attractor
            self.ax.plot(solution[:, 0], solution[:, 1], solution[:, 2], 
                        linewidth=0.5, color='blue', alpha=0.7)
            
            # Mark the initial point
            self.ax.scatter([x0], [y0], [z0], color='red', s=50, label='Initial Point')
            
            # Set labels and title
            self.ax.set_xlabel('X', fontsize=10)
            self.ax.set_ylabel('Y', fontsize=10)
            self.ax.set_zlabel('Z', fontsize=10)
            self.ax.set_title(f'Lorenz Attractor\nσ={sigma}, ρ={rho}, β={beta:.3f}\nInitial: ({x0}, {y0}, {z0})', 
                            fontsize=12)
            self.ax.legend()
            
            # Refresh canvas
            self.canvas.draw()
            
        except ValueError as e:
            print(f"Error: Invalid input - {e}")
        except Exception as e:
            print(f"Error updating plot: {e}")


def main():
    """Main function to run the application"""
    root = tk.Tk()
    app = LorenzGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
