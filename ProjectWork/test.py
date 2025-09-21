import numpy as np
import matplotlib.pyplot as plt
from Basilisk.simulation import spacecraft, simpleBattery, simpleSolarPanel, simplePowerSink
from Basilisk.utilities import macros, SimulationBaseClass
from Basilisk.architecture import messaging

# --- Create simulation ---
scSim = SimulationBaseClass.SimBaseClass()
simTaskName = "simTask"
simProcessName = "simProcess"
dynamicsProcess = scSim.CreateNewProcess(simProcessName)
simulationTimeStep = macros.sec2nano(1.0)  # 1 second step
dynamicsProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))

numSats = 16
scObjects = []
batteries = []
sinks = []
solars = []

for i in range(numSats):
    scObject = spacecraft.Spacecraft()
    scObject.ModelTag = f"sat-{i}"
    scSim.AddModelToTask(simTaskName, scObject)
    scObjects.append(scObject)

    # Battery
    battery = simpleBattery.SimpleBattery()
    battery.ModelTag = f"battery-{i}"
    battery.storageCapacity = 100.0   # Wh
    battery.storedCharge_Init = 80.0  # start at 80%
    scSim.AddModelToTask(simTaskName, battery)
    batteries.append(battery)

    # Solar panel (just keeps adding some power)
    solarPanel = simpleSolarPanel.SimpleSolarPanel()
    solarPanel.ModelTag = f"solarPanel-{i}"
    solarPanel.setPanelParameters([1.0, 0.0, 0.0], 0.001, 0.0001)
    solarPanel.stateInMsg.subscribeTo(scObject.scStateOutMsg)
    scSim.AddModelToTask(simTaskName, solarPanel)
    solars.append(solarPanel)

    # Sun direction
    rawSun = np.array([1.0, 0.2, 0.1])
    sunDir = (rawSun / np.linalg.norm(rawSun)).tolist()
    sunMsgData = messaging.SpicePlanetStateMsgPayload()
    sunMsgData.PositionVector = sunDir
    sunMsg = messaging.SpicePlanetStateMsg().write(sunMsgData)
    solarPanel.sunInMsg.subscribeTo(sunMsg)

    # Power sink (10W baseline load)
    sink = simplePowerSink.SimplePowerSink()
    sink.ModelTag = f"powerSink{i}"
    sink.nodePowerOut = -0.01  # 10 W
    scSim.AddModelToTask(simTaskName, sink)
    sinks.append(sink)

    # Connect
    battery.addPowerNodeToModel(solarPanel.nodePowerOutMsg)
    battery.addPowerNodeToModel(sink.nodePowerOutMsg)

# --- Fault injection (extra 10W for sat-5 after 15 min) ---
fault_sat = 5
faultTime = macros.min2nano(15)
fault_power_kw = 10.0 / 1000.0
scSim.powerSink5=sinks[fault_sat]
scSim.createNewEvent(
    "batteryFaultEvent",
    simulationTimeStep,
    True,
    [f"self.TotalSim.CurrentNanos >= {faultTime}"],
    [
        f"print('*** Fault injected on sat-{fault_sat}: extra 10W drain ***')",
        f"self.{sinks[fault_sat].ModelTag}.nodePowerOut = -{0.01 + fault_power_kw}",
        "self.setEventActivity('batteryFaultEvent', False)"
    ]
)

# --- Set up data logging ---
dataRecorders = []
for bat in batteries:
    dataRec = bat.batPowerOutMsg.recorder()
    scSim.AddModelToTask(simTaskName, dataRec)
    dataRecorders.append(dataRec)

# --- Run simulation ---
simTime = macros.min2nano(30)  # 30 min
scSim.InitializeSimulation()
scSim.ConfigureStopTime(simTime)
scSim.ExecuteSimulation()

# --- Plot battery SOC ---
plt.figure(figsize=(10,6))
timeAxis = dataRecorders[0].times() * macros.NANO2MIN  # convert to minutes

for i, rec in enumerate(dataRecorders):
    soc = rec.storageLevel / rec.storageCapacity * 100  # percent
    if i == fault_sat:
        plt.plot(timeAxis, soc, label=f"Sat-{i} (FAULT)", linewidth=2.5, color="red")
    else:
        plt.plot(timeAxis, soc, label=f"Sat-{i}")

plt.xlabel("Time [min]")
plt.ylabel("Battery SOC [%]")
plt.title("Battery SOC for 16 Satellites (Fault at 15 min on Sat-5)")
plt.legend()
plt.grid(True)
plt.show()
