#
#  ISC License
#
#  Copyright (c) 2016, Autonomous Vehicle Systems Lab, University of Colorado at Boulder
#
#  Permission to use, copy, modify, and/or distribute this software for any
#  purpose with or without fee is hereby granted, provided that the above
#  copyright notice and this permission notice appear in all copies.
#
#  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
#  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
#  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
#  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
#  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
#  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
#  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
import matplotlib.pyplot as plt
import numpy as np
from Basilisk.utilities import RigidBodyKinematics
from Basilisk.utilities import macros as mc
from Basilisk.utilities import orbitalMotion
from Basilisk.utilities import unitTestSupport

# --------------------------------- COMPONENTS & SUBPLOT HANDLING ----------------------------------------------- #
color_x = 'dodgerblue'
color_y = 'salmon'
color_z = 'lightgreen'
m2km = 1.0 / 1000.0

def show_all_plots():
    plt.show()

def clear_all_plots():
    plt.close("all")

def save_all_plots(fileName, figureNames):
    figureList = {}
    numFigures = len(figureNames)
    for i in range(0, numFigures):
        pltName = fileName + "_" + figureNames[i]
        figureList[pltName] = plt.figure(i+1)
    return figureList


def plot3components(timeAxis, vec, id=None):
    plt.figure(id)
    time = timeAxis * mc.NANO2MIN
    plt.xlabel('Time, min')
    plt.plot(time, vec[:, 0], color_x)
    plt.plot(time, vec[:, 1], color_y)
    plt.plot(time, vec[:, 2], color_z)


def plot_sigma(timeAxis, sigma, id=None):
    plot3components(timeAxis, sigma, id)
    plt.legend([r'$\sigma_1$', r'$\sigma_2$', r'$\sigma_3$'])
    plt.ylabel('MRP')


def plot_omega(timeAxis, omega, id=None):
    plot3components(timeAxis, omega, id)
    plt.ylabel('Angular Rate, rad/s')
    plt.legend([r'$\omega_1$', r'$\omega_2$', r'$\omega_3$'])

def subplot_sigma(subplot, timeAxis, sigma, id=None):
    plot3components(timeAxis, sigma, id)
    plt.legend([r'$\sigma_1$', r'$\sigma_2$', r'$\sigma_3$'])
    plt.ylabel('MRP')


def subplot_omega(subplot, timeAxis, omega, id=None):
    plot3components(timeAxis, omega, id)
    plt.ylabel('Angular Rate, rad/s')
    plt.legend([r'$\omega_1$', r'$\omega_2$', r'$\omega_3$'])


# ------------------------------------- MAIN PLOT HANDLING ------------------------------------------------------ #

# def plot_bodyTorque(Lb):
#     plt.figure()
#     plot3components(Lb)
#     plt.title('Body Torque $L_b$')
#     plt.ylabel('Body Torque, $N \cdot m$')
#     plt.legend(['$L_{b,1}$', '$L_{b,2}$', '$L_{b,3}$'])


def plot_controlTorque(timeAxis, Lr, id=None):
    plot3components(timeAxis, Lr, id)
    plt.ylabel(r'Control Torque, $N \cdot m$')
    plt.legend(['$L_{r,1}$', '$L_{r,2}$', '$L_{r,3}$'])
    plt.title('Control Torque $L_r$')
    return


def plot_trackingError(timeAxis, sigma_BR, omega_BR_B, id=None):
    # plt.figure(id)
    plt.subplot(211)
    plot_sigma(timeAxis, sigma_BR, id)
    plt.title(r'Att Error: $\sigma_{BR}$')

    plt.subplot(212)
    #plt.figure(id)
    plot_omega(timeAxis, omega_BR_B, id)
    plt.title(r'Rate Error: $^B{\omega_{BR}}$')
    return


def plot_attitudeGuidance(timeAxis, sigma_RN, omega_RN_N, id=None):
    plot_sigma(timeAxis, sigma_RN, id)
    plt.ylim([-1.0, 1.0])
    plt.title(r'Ref Att: $\sigma_{RN}$')

    plot_omega(timeAxis, omega_RN_N, id)
    plt.title(r'Ref Rate: $^N{\omega_{RN}}$')
    return


def plot_rotationalNav(timeAxis, sigma_BN, omega_BN_B, id=None):
    plt.figure()
    plot_sigma(timeAxis, sigma_BN, id)
    plt.title(r'Sc Att: $\sigma_{BN}$')

    plot_omega(timeAxis, omega_BN_B, id)
    plt.title(r'Sc Rate: $^B{\omega_{BN}}$')
    return


def plot_shadow_fraction(timeAxis, shadow_factor, id=None):
    plt.figure(id)
    plt.plot(timeAxis, shadow_factor)
    plt.xlabel('Time min')
    plt.ylabel('Shadow Fraction')
    return


def plot_sun_point(timeAxis, sunPoint, id=None):
    plot3components(timeAxis, sunPoint, id)
    plt.xlabel('Time')
    plt.ylabel('Sun Point Vec')
    return


def plot_orbit(r_BN, id=None):
    plt.figure(id)
    plt.xlabel('$R_x$, km')
    plt.ylabel('$R_y$, km')
    plt.plot(r_BN[:, 0] * m2km, r_BN[:, 1] * m2km, color_x)
    plt.scatter(0, 0, c=color_x)
    plt.title('Spacecraft Orbit')
    return


def plot_attitude_error(timeLineSet, dataSigmaBR, id=None):
    plt.figure(id)
    fig = plt.gcf()
    ax = fig.gca()
    vectorData = unitTestSupport.pullVectorSetFromData(dataSigmaBR)
    sNorm = np.array([np.linalg.norm(v) for v in vectorData])
    plt.plot(timeLineSet, sNorm,
             color=unitTestSupport.getLineColor(1, 3),
             )
    plt.xlabel('Time [min]')
    plt.ylabel(r'Attitude Error Norm $|\sigma_{B/R}|$')
    ax.set_yscale('log')


def plot_control_torque(timeLineSet, dataLr, id=None, livePlot=False):
    plt.figure(id)
    for idx in range(3):
        plt.plot(timeLineSet, dataLr[:, idx],
                 color=unitTestSupport.getLineColor(idx, 3),
                 label='$L_{r,' + str(idx) + '}$')
    if not livePlot:
        plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('Control Torque $L_r$ [Nm]')


def plot_rate_error(timeLineSet, dataOmegaBR, id=None, livePlot=False):
    plt.figure(id)
    for idx in range(3):
        plt.plot(timeLineSet, dataOmegaBR[:, idx],
                 color=unitTestSupport.getLineColor(idx, 3),
                 label=r'$\omega_{BR,' + str(idx) + '}$')
    if not livePlot:
        plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('Rate Tracking Error [rad/s] ')
    return


def plot_orientation(timeLineSet, vectorPosData, vectorVelData, vectorMRPData, id=None, livePlot=False):
    data = np.empty([len(vectorPosData), 3])
    for idx in range(0, len(vectorPosData)):
        ir = vectorPosData[idx] / np.linalg.norm(vectorPosData[idx])
        hv = np.cross(vectorPosData[idx], vectorVelData[idx])
        ih = hv / np.linalg.norm(hv)
        itheta = np.cross(ih, ir)
        dcmBN = RigidBodyKinematics.MRP2C(vectorMRPData[idx])
        data[idx] = [np.dot(ir, dcmBN[0]), np.dot(itheta, dcmBN[1]), np.dot(ih, dcmBN[2])]
    plt.figure(id)
    labelStrings = (r'$\hat\imath_r\cdot \hat b_1$'
                    , r'${\hat\imath}_{\theta}\cdot \hat b_2$'
                    , r'$\hat\imath_h\cdot \hat b_3$')
    for idx in range(0, 3):
        plt.plot(timeLineSet, data[:, idx],
                 color=unitTestSupport.getLineColor(idx, 3),
                 label=labelStrings[idx])
    if not livePlot:
        plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('Orientation Illustration')


def plot_rw_cmd_torque(timeData, dataUsReq, numRW, id=None, livePlot=False):
    plt.figure(id)
    for idx in range(3):
        plt.plot(timeData, dataUsReq[:, idx],
                 '--',
                 color=unitTestSupport.getLineColor(idx, numRW),
                 label=r'$\hat u_{s,' + str(idx) + '}$')
    if not livePlot:
        plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('RW Motor Torque (Nm)')


def plot_rw_cmd_actual_torque(timeData, dataUsReq, dataRW, numRW, id=None, livePlot=False):
    """compare commanded and actual RW motor torques"""
    plt.figure(id)
    for idx in range(numRW):
        plt.plot(timeData, dataUsReq[:, idx],
                 '--',
                 color=unitTestSupport.getLineColor(idx, numRW),
                 label=r'$\hat u_{s,' + str(idx) + '}$')
        plt.plot(timeData, dataRW[idx],
                 color=unitTestSupport.getLineColor(idx, numRW),
                 label='$u_{s,' + str(idx) + '}$')
    plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('RW Motor Torque (Nm)')


def plot_rw_speeds(timeData, dataOmegaRW, numRW, id=None, livePlot=False):
    plt.figure(id)
    for idx in range(numRW):
        plt.plot(timeData, dataOmegaRW[:, idx] / mc.RPM,
                 color=unitTestSupport.getLineColor(idx, numRW),
                 label=r'$\Omega_{' + str(idx) + '}$')
    if not livePlot:
        plt.legend(loc='upper right')
    plt.xlabel('Time [min]')
    plt.ylabel('RW Speed (RPM) ')

def plot_rw_friction(timeData, dataFrictionRW, numRW, dataFaultLog=[],  id=None, livePlot=False):
    plt.figure(id)
    for idx in range(numRW):
        plt.plot(timeData, dataFrictionRW[idx],
                 color=unitTestSupport.getLineColor(idx, numRW),
                 label=r'$RW_{' + str(idx+1) + '} Friction$')
    if dataFaultLog:
        # fourth column of dataFaultLog is the fault times
        plt.scatter([row[3] for row in dataFaultLog], np.zeros(len(dataFaultLog)), marker="x", color=(1,0,0),
            label='Faults')
    if not livePlot:
        plt.legend(loc='upper right')
    plt.xlabel('Time [min]')
    plt.ylabel('RW Static Friction ')


def plot_planet(oe, planet):
    b = oe.a * np.sqrt(1 - oe.e * oe.e)
    plt.figure(figsize=np.array((1.0, b / oe.a)) * 4.75, dpi=100)
    plt.axis(np.array([-oe.a, oe.a, -b, b]) / 1000 * 1.75)
    # draw the planet
    fig = plt.gcf()
    ax = fig.gca()
    planetColor = '#008800'
    planetRadius = planet.radEquator / 1000
    ax.add_artist(plt.Circle((0, 0), planetRadius, color=planetColor))
    plt.xlabel('$i_e$ Cord. [km]')
    plt.ylabel('$i_p$ Cord. [km]')
    plt.grid()


def plot_peri_and_orbit(oe, mu, r_BN_N, v_BN_N, id=None):
    # draw the actual orbit
    rData = []
    fData = []
    p = oe.a * (1 - oe.e * oe.e)
    for idx in range(0, len(r_BN_N)):
        oeData = orbitalMotion.rv2elem(mu, r_BN_N[idx], v_BN_N[idx])
        rData.append(oeData.rmag)
        fData.append(oeData.f + oeData.omega - oe.omega)
    plt.figure(id)
    plt.plot(rData * np.cos(fData) / 1000, rData * np.sin(fData) / 1000, color='#aa0000', linewidth=3.0)
    # draw the full osculating orbit from the initial conditions
    fData = np.linspace(0, 2 * np.pi, 100)
    rData = []
    for idx in range(0, len(fData)):
        rData.append(p / (1 + oe.e * np.cos(fData[idx])))
    plt.plot(rData * np.cos(fData) / 1000, rData * np.sin(fData) / 1000, '--', color='#555555')


def plot_rel_orbit(timeData, r_chief, r_deputy, id=None, livePlot=False):
    plt.figure(id)
    x = np.array(r_chief[:, 0]) - np.array(r_deputy[:, 0])
    y = np.array(r_chief[:, 1]) - np.array(r_deputy[:, 1])
    z = np.array(r_chief[:, 2]) - np.array(r_deputy[:, 2])
    plt.plot(timeData, x, label="x")
    plt.plot(timeData, y, label="y")
    plt.plot(timeData, z, label="z")
    if not livePlot:
        plt.legend()
    plt.grid()


def plot_target_visibility(timeData, r_BN_N, targets):
    """
    Plot target visibility information

    Args:
        timeData: Array of simulation times in minutes
        r_BN_N: Array of spacecraft positions in inertial frame
        targets: List of target dictionaries with 'name', 'lat', 'lon' keys
        
    Returns:
        fig: The matplotlib figure object for display or saving
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from Basilisk.utilities import macros
    
    # Calculate visibility for each target over time
    visibilityData = []
    targetNames = []
    
    for target in targets:
        targetNames.append(target["name"])
        visibility = []
        
        for i, t in enumerate(timeData):
            # Convert time to hours for Earth rotation calculation
            time_hours = t / 60.0  # minutes to hours
            earth_rotation = time_hours * 15.0 * macros.D2R  # Earth rotates ~15 degrees per hour
            
            # Create rotation matrix
            c_rot = np.cos(earth_rotation)
            s_rot = np.sin(earth_rotation)
            R_N_E = np.array([
                [c_rot, -s_rot, 0],
                [s_rot, c_rot, 0],
                [0, 0, 1]
            ])
            
            # Convert lat/lon to ECEF coordinates
            lat_rad = target["lat"] * macros.D2R
            lon_rad = target["lon"] * macros.D2R
            r_Earth = 6371000.0  # Earth radius in meters
            
            # Target position in Earth-fixed frame
            r_Target_E = np.array([
                (r_Earth) * np.cos(lat_rad) * np.cos(lon_rad),
                (r_Earth) * np.cos(lat_rad) * np.sin(lon_rad),
                (r_Earth) * np.sin(lat_rad)
            ])
            
            # Convert to inertial frame
            r_Target_N = R_N_E @ r_Target_E
            
            # Check if target is visible (above horizon)
            sc_pos = r_BN_N[i]
            sc_to_target = r_Target_N - sc_pos
            sc_to_center = -sc_pos
            
            # If angle between these vectors is less than 90 degrees, target is visible
            cos_angle = np.dot(sc_to_target, sc_to_center) / (np.linalg.norm(sc_to_target) * np.linalg.norm(sc_to_center))
            
            # Set a value that represents visibility (1 for visible, 0 for not visible)
            if cos_angle < 0:  # Angle > 90 degrees, target is visible
                visibility.append(1)
            else:
                visibility.append(0)
        
        visibilityData.append(visibility)
    
    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    
    for i, vis in enumerate(visibilityData):
        ax.plot(timeData, [i+1]*len(timeData), 'k-', alpha=0.2)  # baseline
        visible_segments = []
        start_idx = None
        
        # Find continuous visible segments
        for j, v in enumerate(vis):
            if v == 1 and start_idx is None:
                start_idx = j
            elif v == 0 and start_idx is not None:
                visible_segments.append((timeData[start_idx], timeData[j-1], i+1))
                start_idx = None
        
        # Handle case where visibility extends to the end
        if start_idx is not None:
            visible_segments.append((timeData[start_idx], timeData[-1], i+1))
        
        # Plot visible segments
        for start, end, level in visible_segments:
            ax.plot([start, end], [level, level], 'g-', linewidth=3)
    
    ax.set_yticks(range(1, len(targetNames)+1))
    ax.set_yticklabels(targetNames)
    ax.set_xlabel('Time [min]')
    ax.set_title('Target Visibility Timeline')
    ax.grid(True)
    
    return fig

def pull_outputs(self, showPlots):
    """Process and plot simulation outputs"""
    # FSW process outputs, remove first data point as it is before FSW is called
    attErrRec = self.msgRecList[self.attGuidName]

    sigma_BR = np.delete(attErrRec.sigma_BR, 0, 0)
    omega_BR_B = np.delete(attErrRec.omega_BR_B, 0, 0)
    
    num_RW = 4
    RW_speeds = np.delete(self.rwSpeedRec.wheelSpeeds[:, range(num_RW)], 0, 0)
    RW_friction = []
    for i in range(num_RW):
        RW_friction.append(np.delete(self.rwLogs[i].u_f, 0, 0))

    # Plot results
    BSK_plt.clear_all_plots()
    timeData = np.delete(attErrRec.times(), 0, 0) * macros.NANO2MIN
    BSK_plt.plot_attitude_error(timeData, sigma_BR, id=1)
    BSK_plt.plot_rate_error(timeData, omega_BR_B, id=2)
    BSK_plt.plot_rw_speeds(timeData, RW_speeds, num_RW, id=3)
    BSK_plt.plot_rw_friction(timeData, RW_friction, num_RW, self.DynModels.RWFaultLog, id=4)
    
    # Create the target visibility figure (using an explicit figure number)
    vis_fig = self.plot_target_visibility(timeData)
    
    # Make sure we have a local figureList dict
    figureList = {}
    
    if showPlots:
        BSK_plt.show_all_plots()
        plt.figure(vis_fig.number)  # Ensure visibility figure is shown
        plt.show()
    else:
        fileName = os.path.basename(os.path.splitext(__file__)[0])
        figureNames = ["attitudeErrorNorm", "rateError", "RWSpeeds", "RWFriction"]
        figureList = BSK_plt.save_all_plots(fileName, figureNames)
        
        # Manually add the visibility figure
        pltName = fileName + "_targetVisibility"
        figureList[pltName] = vis_fig
        vis_fig.savefig(pltName + ".png")

    return figureList