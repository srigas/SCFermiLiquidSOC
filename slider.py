from numpy import pi, sin
import numpy as np
import cmath
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

mycol = (0.13333, 0.35294, 0.38824) # Teal color

plt.rc('font', size=12)
plt.rc('axes', titlesize=18)
plt.rc('axes', labelsize=20)
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('legend', fontsize=12)
plt.rc('figure', titlesize=18)

# Constants
HBAR = 1
M = 0.5
IMAG = 1j
PI = cmath.pi

def EPSONE(KFERMI,BETA,DELTA,RASHBA):
    singlen = (HBAR**2/(2*M))*(KRANGE**2-KFERMI**2)
    rasterm = RASHBA*KRANGE
    subsqrt = np.sqrt((singlen*rasterm)**2 + (singlen**2 + DELTA**2)*BETA**2)
    return np.sqrt(singlen**2 + rasterm**2 + DELTA**2 + BETA**2 + 2*subsqrt)

def EPSTWO(KFERMI,BETA,DELTA,RASHBA):
    singlen = (HBAR**2/(2*M))*(KRANGE**2-KFERMI**2)
    rasterm = (RASHBA)*(KRANGE)
    subsqrt = np.sqrt((singlen*rasterm)**2 + (singlen**2 + DELTA**2)*BETA**2)
    return np.sqrt(singlen**2 + rasterm**2 + DELTA**2 + BETA**2 - 2*subsqrt)

axis_color = 'white'

fig = plt.figure(figsize=(10, 8))
fig.suptitle('Energy eigenvalues for superconducting rashba liquid', fontsize=16)
ax = fig.add_subplot(111)
ax.set_xlabel(r'$k$', labelpad=14)
ax.set_ylabel(r'$E$', labelpad=12)
ax.axhline(y=0.0, color='black', linestyle='-', linewidth=0.5)
ax.axvline(x=0.0, color='black', linestyle='-', linewidth=0.5)

# leave some space for the sliders and buttons
fig.subplots_adjust(bottom=0.4)

KRANGE = np.arange(-4.0, 4.0, 0.005)
KFERMI_0 = 1.0
DELTA_0 = 0.0
RASHBA_0 = 0.0
BETA_0 = 0.0

# Draw the plot for the default values
# The lines are four: 2 negative and 2 positive energy ones (particle-hole symmetry)
[lineone] = ax.plot(KRANGE, EPSONE(KFERMI_0, BETA_0, DELTA_0, RASHBA_0), linewidth=1, color=mycol)
[linetwo] = ax.plot(KRANGE, EPSTWO(KFERMI_0, BETA_0, DELTA_0, RASHBA_0), linewidth=1, color=mycol)
[linethree] = ax.plot(KRANGE, -EPSONE(KFERMI_0, BETA_0, DELTA_0, RASHBA_0), linewidth=1, color=mycol)
[linefour] = ax.plot(KRANGE, -EPSTWO(KFERMI_0, BETA_0, DELTA_0, RASHBA_0), linewidth=1, color=mycol)

vone = ax.axvline(KFERMI_0, color='black', linestyle='--', linewidth=0.5)
vtwo = ax.axvline(-KFERMI_0, color='black', linestyle='--', linewidth=0.5)

ax.set_xlim([-2, 2])
ax.set_ylim([-5, 5])

# Add four sliders for tweaking the parameters

kfermi_slider_ax  = fig.add_axes([0.25, 0.25, 0.65, 0.03], facecolor=axis_color)
kfermi_slider = Slider(kfermi_slider_ax, r'$k_F$', 0.0, 4.0, valinit=KFERMI_0)

beta_slider_ax = fig.add_axes([0.25, 0.2, 0.65, 0.03], facecolor=axis_color)
beta_slider = Slider(beta_slider_ax, r'$B$', 0.0, 10.0, valinit=BETA_0)

delta_slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axis_color)
delta_slider = Slider(delta_slider_ax, r'$\Delta$', 0.0, 10.0, valinit=DELTA_0)

rashba_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axis_color)
rashba_slider = Slider(rashba_slider_ax, r'$\alpha_R$', 0.0, 10.0, valinit=RASHBA_0)

# Upon the change of any slider's value
def sliders_on_changed(val):
    lineone.set_ydata(EPSONE(kfermi_slider.val, beta_slider.val, delta_slider.val, rashba_slider.val))
    linetwo.set_ydata(EPSTWO(kfermi_slider.val, beta_slider.val, delta_slider.val, rashba_slider.val))
    linethree.set_ydata(-EPSONE(kfermi_slider.val, beta_slider.val, delta_slider.val, rashba_slider.val))
    linefour.set_ydata(-EPSTWO(kfermi_slider.val, beta_slider.val, delta_slider.val, rashba_slider.val))

    vone.set_xdata(kfermi_slider.val)
    vtwo.set_xdata(-kfermi_slider.val)

    fig.canvas.draw_idle()
kfermi_slider.on_changed(sliders_on_changed)
beta_slider.on_changed(sliders_on_changed)
delta_slider.on_changed(sliders_on_changed)
rashba_slider.on_changed(sliders_on_changed)

# Button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    kfermi_slider.reset()
    beta_slider.reset()
    delta_slider.reset()
    rashba_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

plt.show()