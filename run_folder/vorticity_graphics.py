from PyQt5 import QtGui  # (the example applies equally well to PySide)
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtWidgets
from PyQt5.QtCore import Qt
import pyqtgraph.opengl as gl

animation_on = False
def stop_visual():
    global animation_on
    animation_on = False

def start_visual():
    global animation_on
    animation_on = True

def reset_visual():
    global i, animation_on
    i = 0
    animation_on = False
    iterationvalue.setText('T= {}'.format(i))
    val1 = 0
    colors  = np.genfromtxt("./data/colours.txt", skip_header=i*len(faces), max_rows=len(faces), usecols = (val1, val1+1, val1+2, val1+3))
    m1.setMeshData(vertexes=verts, faces=faces, faceColors=colors, smooth=False, drawEdges=True, edgeColor=(1,1,1,1)) 
    vels = np.genfromtxt("./data/velocities.txt", skip_header=i*len(faces), max_rows=len(faces), usecols = (0, 1, 2)) 
    for l in range(int(len(vels))):
        lines[l].setData(pos=np.array([cell_centers[l], vels[l]]))


def change_particle_size():
    global i, colors, verts, faces
    val1 = 0
    colors  = np.genfromtxt("./data/colours.txt", skip_header=i*len(faces), max_rows=len(faces), usecols = (val1, val1+1, val1+2, val1+3))
    m1.setMeshData(vertexes=verts, faces=faces, faceColors=colors, smooth=False, drawEdges=True, edgeColor=(1,1,1,1))
    
## Always start by initializing Qt (only once per application)
app = QtWidgets.QApplication([])

## Define a top-level widget to hold everything
w = QtWidgets.QWidget()
w.setWindowTitle("vorticity")
## Create some widgets to be placed inside
start_btn = QtWidgets.QPushButton('start')
start_btn.clicked.connect(start_visual)
stop_btn = QtWidgets.QPushButton('stop')
stop_btn.clicked.connect(stop_visual)
reset_btn = QtWidgets.QPushButton('reset')
reset_btn.clicked.connect(reset_visual)
## Create a grid layout to manage the widgets size and position
layout = QtWidgets.QGridLayout()
w.setLayout(layout)
w.resize(875, 800)
w.show()
## Add widgets to the layout in their proper positions
i = 0
iterationvalue = QtWidgets.QLabel('T= {}'.format(i))
layout.addWidget(iterationvalue, 1, 1)
layout.addWidget(start_btn, 2, 1)
layout.addWidget(stop_btn, 3, 1)
layout.addWidget(reset_btn, 4, 1)

grid = gl.GLViewWidget()
grid.setCameraPosition(elevation=90, azimuth=0) ##, distance=0)

verts = np.loadtxt("./data/vertices.txt")
faces = np.loadtxt("./data/faces.txt", dtype=int)
colors = np.genfromtxt("./data/colours.txt", max_rows=len(faces), usecols = (0, 1, 2, 3))
## Mesh item will automatically compute face normals.
m1 = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=colors, smooth=False, drawEdges=True, edgeColor=(1,1,1,1))
m1.translate(-0.75, -0.25, 0.0)
m1.setGLOptions('additive')
m1.rotate(90, 0,0,1)

lines = []
vels = np.genfromtxt("./data/velocities.txt", max_rows=len(faces), usecols = (0, 1, 2))
cell_centers = np.loadtxt("./data/cellcenters.txt")
for l in range(int(len(vels))):
    lines.append(gl.GLLinePlotItem(pos=np.array([cell_centers[l], vels[l]]), color=(255,255,255,255), antialias=True))
    lines[l].translate(0.25, -0.75, 0.0)
    grid.addItem(lines[l])

grid.addItem(m1)

gw = pg.GradientEditorItem(orientation='right')
# load predefined color gradient
gw.restoreState({'ticks': [(0.0, (0, 0, 255, 255)), (0.5, (0, 255, 0, 255)), (1.0, (255, 0, 0, 255))], 'mode': 'rgb'})

ax = pg.AxisItem('left')
ax.setRange(0.0, 1.0)
cb = pg.GraphicsLayoutWidget()

cb.addItem(ax)
cb.addItem(gw)
cb.resize(100, 700)
cb.show()
layout.addWidget(cb, 0, 1)
layout.addWidget(grid, 0, 0)
layout.setColumnStretch(0, 1)

def updateData():
    global colors, i, animation_on, verts, faces, total_frames

    QtCore.QTimer.singleShot(100, updateData)

    if (animation_on):
        iterationvalue.setText('T= {}'.format(i))
        val1 = 0
        colors  = np.genfromtxt("./data/colours.txt", skip_header=i*len(faces), max_rows=len(faces), usecols = (val1, val1+1, val1+2, val1+3))
        i = (i+1)%total_frames
        m1.setMeshData(vertexes=verts, faces=faces, faceColors=colors, smooth=False, drawEdges=True, edgeColor=(1,1,1,1))   
        m1.setMeshData(vertexes=verts, faces=faces, faceColors=colors, smooth=False, drawEdges=False)
        vels = np.genfromtxt("./data/velocities.txt", skip_header=i*len(faces), max_rows=len(faces), usecols = (0, 1, 2))
        for l in range(int(len(vels))):
            lines[l].setData(pos=np.array([cell_centers[l], vels[l]]))
        grid.grabFramebuffer().save('./images/fileName_%d.png' % i)

updateData()

## Start Qt event loop unless running in interactive mode.

if __name__ == '__main__':
    import sys
    total_frames = int(sys.argv[1])
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtWidgets.QApplication.instance().exec_()
