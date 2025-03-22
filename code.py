import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib.widgets import Slider

def update_plot(val):
    theta = theta_slider.val  # Get the current slider value
    ax1.cla()
    ax2.cla()
    ax3.cla()

    # Define the sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 50)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    ax1.plot_surface(x, y, z, color='c', alpha=0.1, edgecolor='gray')
    ax1.plot_wireframe(x, y, z, color="gray", alpha=0.2, linewidth=0.5)  # Improved visibility

    # Axes
    ax1.quiver(0, 0, 0, 1, 0, 0, color="red", linewidth=2)    # x-axis
    ax1.quiver(0, 0, 0, 0, 1, 0, color="green", linewidth=2)  # y-axis
    ax1.quiver(0, 0, 0, 0, 0, 1, color="blue", linewidth=2)   # z-axis

    # Compute qubit position
    x_q = np.sin(theta) * np.cos(phi)
    y_q = np.sin(theta) * np.sin(phi)
    z_q = np.cos(theta)

    # Plot qubit state vector
    ax1.quiver(0, 0, 0, x_q, y_q, z_q, color="black", linewidth=4, arrow_length_ratio=0.2)
    
    # Labels
    ax1.text(1.2, 0, 0, '|+⟩', color='red', fontsize=12)
    ax1.text(0, 1.2, 0, '|i⟩', color='green', fontsize=12)
    ax1.text(0, 0, 1.2, '|0⟩', color='blue', fontsize=12)
    ax1.text(0, 0, -1.2, '|1⟩', color='blue', fontsize=12)

    ax1.set_xlim([-1, 1])
    ax1.set_ylim([-1, 1])
    ax1.set_zlim([-1, 1])
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_zticks([])
    ax1.set_title("Bloch Sphere Representation")

    # Probability Graph (2D)
    P_0 = np.cos(theta / 2) ** 2  # Probability of |0>
    P_1 = np.sin(theta / 2) ** 2  # Probability of |1>
    states = ['|0>', '|1>']
    probabilities = [P_0, P_1]
    bars = ax2.bar(states, probabilities, color=['blue', 'red'])
    ax2.set_ylim(0, 1)
    ax2.set_ylabel('Probability')
    ax2.set_title("Measurement Probabilities")
    
    # Add labels to bars
    for bar, prob in zip(bars, probabilities):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05, f"{prob:.2f}", ha='center', fontsize=12, fontweight='bold')

    # Calculation Details
    ax3.axis('off')
    calculation_text = (
        f"Theta: {theta:.2f}\n"
        f"Phi: {phi:.2f}\n"
        f"\nCalculations:\n"
        f"P(|0⟩) = cos²(θ/2) = cos²({theta/2:.2f}) = {P_0:.2f}\n"
        f"P(|1⟩) = sin²(θ/2) = sin²({theta/2:.2f}) = {P_1:.2f}\n"
        
        
        f"\nQubit Position Matrix:\n"
        f"[{x_q:.2f}]\n"
        f"[{y_q:.2f}]\n"
        f"[{z_q:.2f}]\n"
    )
    ax3.text(0, 0.6, calculation_text, fontsize=12, verticalalignment='top', horizontalalignment='left')
    
    # Display quantum state equation
    state_text = (
        r"$|\psi\rangle = \cos(\frac{\theta}{2}) |0\rangle + e^{i\phi} \sin(\frac{\theta}{2}) |1\rangle$"
    )
    ax3.text(0, 0.2, state_text, fontsize=14, verticalalignment='center', horizontalalignment='left', color='purple')
    
    plt.draw()

def rotate_view(frame):
    ax1.view_init(elev=30, azim=frame % 360)  # Continuous rotation
    plt.draw()

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
ax1 = fig.add_subplot(131, projection='3d')
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

theta = np.pi / 4
phi = np.pi / 4

# Add interactive theta slider
fig.subplots_adjust(bottom=0.2)
ax_slider = plt.axes([0.2, 0.05, 0.65, 0.03])
theta_slider = Slider(ax_slider, 'Theta', 0, np.pi, valinit=theta)

theta_slider.on_changed(update_plot)

ani = animation.FuncAnimation(fig, rotate_view, frames=360, interval=50, repeat=True)  # Smooth rotation
plt.show()
