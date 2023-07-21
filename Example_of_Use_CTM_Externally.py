import Cell_Transmission_Model as CTM
from datetime import datetime

start = datetime.now()
# Initialize simulation instance
sim = CTM.initializeCTM()
for t in range(sim.total_steps):
    # Run the simulation step by step
    CTM.simulation_run_step(sim)

    # Dynamiclly and temporally change parameters of simulation
    if t >= 100 and t <= 300:
        CTM.Cell.idcase['A0.100000001.C3'].qmax = 600
end = datetime.now()
print("Elapsed Time:", end - start)

# Write results to csv files
for key in CTM.Corridor.idcase:
    CTM.Corridor.idcase[key].printResults()
