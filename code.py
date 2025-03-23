import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

def calculate_phi(alpha, beta):
    """Calculate the azimuthal angle phi from a qubit state."""
    if np.isclose(beta.real, 0) and np.isclose(beta.imag, 0):
        return 0  # Avoid division by zero
    return np.arctan2(beta.imag, beta.real)

# Get user input for alpha and beta
try:
    alpha_real = float(input("Enter real part of alpha: "))
    alpha_imag = float(input("Enter imaginary part of alpha: "))
    beta_real = float(input("Enter real part of beta: "))
    beta_imag = float(input("Enter imaginary part of beta: "))
    alpha = complex(alpha_real, alpha_imag)
    beta = complex(beta_real, beta_imag)
    
    # Normalize alpha and beta
    norm = np.sqrt(abs(alpha)**2 + abs(beta)**2)
    alpha /= norm
    beta /= norm
except ValueError:
    print("Invalid input, using default values for alpha and beta.")
    alpha = np.cos(np.pi / 8)
    beta = np.sin(np.pi / 8) * np.exp(1j * (np.pi / 4))

# Compute theta and phi dynamically
theta = 2 * np.arccos(abs(alpha))
phi = calculate_phi(alpha, beta)

# Create figure
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
ax1 = fig.add_subplot(131, projection='3d')
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

def plot_bloch_sphere(ax, angle):
    ax.clear()
    # Define the sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(np.pi, 0, 50)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='c', alpha=0.1, edgecolor='gray')
    ax.plot_wireframe(x, y, z, color="gray", alpha=0.2, linewidth=0.5)

    # Axes
    ax.quiver(0, 0, 0, 1, 0, 0, color="red", linewidth=2)
    ax.quiver(0, 0, 0, 0, 1, 0, color="green", linewidth=2)
    ax.quiver(0, 0, 0, 0, 0, 1, color="blue", linewidth=2)

    # Compute qubit position
    x_q = np.sin(theta) * np.cos(phi)
    y_q = np.sin(theta) * np.sin(phi)
    z_q = np.cos(theta)

    # Plot qubit state vector
    ax.quiver(0, 0, 0, x_q, y_q, z_q, color="black", linewidth=4, arrow_length_ratio=0.2)

    # Labels
    ax.text(1.2, 0, 0, '|+⟩', color='red', fontsize=12)
    ax.text(0, 1.2, 0, '|i⟩', color='green', fontsize=12)
    ax.text(0, 0, 1.2, '|0⟩', color='blue', fontsize=12)
    ax.text(0, 0, -1.2, '|1⟩', color='blue', fontsize=12)

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_title("Bloch Sphere Representation")
    ax.view_init(elev=20, azim=angle)

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
    f"|Ψ⟩ = [{alpha.real:.2f} + {alpha.imag:.2f}i] |0⟩ + [{beta.real:.2f} + {beta.imag:.2f}i] |1⟩\n"
    f"Normalization: |α|² + |β|² = {abs(alpha)**2:.2f} + {abs(beta)**2:.2f} = 1\n"
    f"Theta: θ = 2 arccos(|α|) = 2 arccos({abs(alpha):.2f}) = {theta:.2f} rad\n"
    f"Phi: φ = arctan2(Im(β), Re(β)) = arctan2({beta.imag:.2f}, {beta.real:.2f}) = {phi:.2f} rad\n"
    f"P(|0⟩) = cos²(θ/2) = {P_0:.2f}\n"
    f"P(|1⟩) = sin²(θ/2) = {P_1:.2f}"
)
ax3.text(0, 0.5, calculation_text, fontsize=12, verticalalignment='center', horizontalalignment='left')

# Rotate and update plot
for angle in range(0, 360, 2):
    plot_bloch_sphere(ax1, angle)
    plt.pause(0.05)  # Adjust speed

plt.show()
