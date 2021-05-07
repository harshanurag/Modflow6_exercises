import os
import numpy
import matplotlib.pyplot as plt
import flopy

os.chdir('/home/harsh/Desktop/modflow_testcase/test_1')
name = "mf6_test_1"
h1 = 100 #head at the side
h2 = 90 #head at the center, top layer
Nlay = 10 #number of layers
N = 101 #number of rows and columns
L = 400 #length of the side of the model
H = 50 #aquifer thickness
k = 1 #hydraulic conductivity

'''
Simulations in mf6 made up of:
- Temporal discretization - here TDIS packaege
- one or more models - GWF model
- zero or more exchanges
- solutions - here iterative model solution (IMS)

'''
mf6_exe = '/home/harsh/work/modflow/mf6.2.1/bin/mf6'
#creating the simulation object
sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6_exe, version="mf6", sim_ws=".")

#creating the TDIS object
tdis = flopy.mf6.ModflowTdis(sim, pname="tdis", time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])

#creating the IMS package object
ims = flopy.mf6.ModflowIms(sim, pname="ims", complexity="SIMPLE")

#creating the groundwater flow model object
model_nam_file = "{}.nam".format(name)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, model_nam_file=model_nam_file)

###
#After creating the gwf model object, can now build the gwf model
#by adding packages to it that would define the model's characteristics
###

#creating the discretization package
#using H and Nlay for indicating top and bottom value
#using L and N for calculating number of rows and columns
bot = np.linspace(-H/Nlay, -H, Nlay)
delrow = delcol = L/(N-1)
dis = flopy.mf6.ModflowGwfdis(
    gwf,
    nlay=Nlay,
    nrow=N,
    ncol=N,
    delr=delrow,
    delc=delcol,
    top=0.0,
    botm=bot,
)

#create the initial condition package
start = h1 * np.ones((Nlay, N, N))
ic = flopy.mf6.ModflowGwfic(gwf, pname="ic", strt=start)

#create the node property flow package
npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k, save_flows=True)

#create the constant head (chd) package
chd_rec = []
chd_rec.append(((0, int(N / 4), int(N / 4)), h2))
for layer in range(0, Nlay):
    for row_col in range(0, N):
        chd_rec.append(((layer, row_col, 0), h1))
        chd_rec.append(((layer, row_col, N - 1), h1))
        if row_col != 0 and row_col != N - 1:
            chd_rec.append(((layer, 0, row_col), h1))
            chd_rec.append(((layer, N - 1, row_col), h1))
chd = flopy.mf6.ModflowGwfchd(
    gwf,
    maxbound=len(chd_rec),
    stress_period_data=chd_rec,
    save_flows=True,
)


#CHD package stored the constant head in a structured array;
#can get a pointer to the recarray for the first stress period (iper=0) as follows:
iper = 0
ra = chd.stress_period_data.get_data(key=iper)

#create the output control package
headfile = "{}.hds".format(name)
head_filerecord = [headfile]
budgetfile = "{}.cbb".format(name)
budget_filerecord = [budgetfile]
saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")]
printrecord = [("HEAD", "LAST")]
oc = flopy.mf6.ModflowGwfoc(
    gwf,
    saverecord=saverecord,
    head_filerecord=head_filerecord,
    budget_filerecord=budget_filerecord,
    printrecord=printrecord,
)


#after the model is created, you can create the input dataset and run the simulations

#writing the datasets
#output is currently in home folder
sim.write_simulation()


#run the simulation
success, buff = sim.run_simulation()
if not success:
    raise Exception("MODFLOW 6 did not terminate normally.")
    
    

#TAKING A PEAK AT THE OUTPUT FILES
#map of layer 1
hds = flopy.utils.binaryfile.HeadFile(headfile)
h = hds.get_data(kstpkper=(0,0))
x = y = np.linspace(0, L, N)
y = y[::-1]
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1, 1, 1, aspect="equal")
c = ax.contour(x, y, h[0], np.arange(90,100.1, 0.2), colors="black")
plt.clabel(c, fmt="%2.1f")

#map of layer 10
x = y = np.linspace(0, L, N)
y = y[::-1]
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
c = ax.contour(x, y, h[-1], np.arange(90,100.1,0.2), colors='black')
plt.clabel(c, fmt="%1.1f")

#plot cross section along row 51
z = np.linspace(-H/Nlay/2, -H + H/Nlay/2, Nlay)
fig = plt.figure(figsize=(5,2.5))
ax = fig.add_subplot(1,1,1,aspect='equal')
c = ax.contour(x, z, h[:, 50, :], np.arange(90,100.1,0.2), colors='black')
plt.clabel(c, fmt="%1.1f")


#USING FLOPY'S PlotMapView() capabilities in MF6

#creating MF5 ibound array for plotting the locations of the constant head
ibd = np.ones((Nlay, N, N), dtype=int)
for k, i, j in ra['cellid']:
    ibd[k, i, j] = -1

#plotting map of layers 1 and 10
fig, axes = plt.subplots(2, 1, figsize=(6, 12), constrained_layout=True)

ax = axes[0]
ax.set_title("Model layer 1")
modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax)
quadmesh = modelmap.plot_ibound(ibound=ibd)
linecollection = modelmap.plot_grid(lw=0.5, color="0.5")
contours = modelmap.contour_array(h[0], levels=np.arange(90,100,0.2), colors='black')
ax.clabel(contours, fmt="%2.1f")
# second subplot
ax = axes[1]
ax.set_title("Model Layer {}".format(Nlay))
modelmap = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=Nlay - 1)
quadmesh = modelmap.plot_ibound(ibound=ibd)
linecollection = modelmap.plot_grid(lw=0.5, color="0.5")
pa = modelmap.plot_array(h[0])
contours = modelmap.contour_array(
    h[0], levels=np.arange(90, 100.1, 0.2), colors="black"
)
cb = plt.colorbar(pa, shrink=0.5, ax=ax)
ax.clabel(contours, fmt="%2.1f")