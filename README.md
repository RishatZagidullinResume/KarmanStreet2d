# Vorticity equation solver in 2D using finite volume approach.

## Code description

The code solves the voritcity equation in 2D. The problem set is the one that generates the Karman street. The space domain is a triangulated mesh. We use finite volume approach to find the numerical solution.

* Parameters:
`dr` - approximate space step of the triangulation, 
`dt` - time step, 
`TIME` - total number of iterations, 
`TOTAL_FRAMES` - number of iterations that will be dumped to a file,
`length_x` - length of the computational domain, 
`length_y` - height of the computational domain, 
`init velocity_x` - velocity of the flux.

* To generate the mesh domain we use **Fade2D** library.

* The numerical solver takes a config file as a sole command line argument. See __input_2d__ for reference.

* The __Makefile__ outputs the numerical solver executable (__2D_CPU_CPU.exe__).

## Visualization and demo

* For visualization purposes see __vorticity_graphics.py__ script. It uses **PyQt5** for a simple GUI interactable and **pyopengl** for visualization. The script accepts 1 command line argument. For example, `python3 2d_py.py 100` meaning: 100 time iterations are going to be visualized.

* Below is an example of the numerical solution.

![image.gif](/image.gif)


