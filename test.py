from mycosimulator import MyceliumSimulator, section

simulator = MyceliumSimulator()

simulator.add_sec(section((0,0,0),(0,0,0),(1,0,0)))
simulator.add_sec(section((0,0,0),(0,0,0),(-1,0,0)))

simulator.run_tropisms(30, a=5, limit_f1=0.1, Ic=0.01, Id=0.01, k=0.5)