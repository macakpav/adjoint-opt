import numpy as np
#import matplotlib.pyplot as pyplot

if __name__ == "__main__":
#-----------------------------------------------------------------------------------------
    filename = "GHS1_4760_Camber_Line.csv"
    output_filename = "airfoilcpsBsplines0"
    param_t = 0.0476

    #nCPsU=14 # number of coordinates in the file minus 2 has to be divisible by this number
    nCPsV=6
    nCPsW=3
#-----------------------------------------------------------------------------------------

    print("Generating control points file...")
    file = open(filename)
    camber_pts = np.loadtxt(file, delimiter=",")
    
    print("len(camber_pts)")
    print(len(camber_pts))
    print("camber_pts")
    print(camber_pts)
   
    lead_pt = camber_pts[0]+(camber_pts[0] - camber_pts[1])/np.linalg.norm(camber_pts[0]  - camber_pts[1])*0.016
    lead_pt2 = camber_pts[0]+(camber_pts[0] - camber_pts[1])/np.linalg.norm(camber_pts[0]  - camber_pts[1])*0.008
    trail_pt = camber_pts[-1]+(camber_pts[-1] - camber_pts[-2])/np.linalg.norm(camber_pts[-1] - camber_pts[-2])*0.016
    trail_pt2 = camber_pts[-1]+(camber_pts[-1] - camber_pts[-2])/np.linalg.norm(camber_pts[-1] - camber_pts[-2])*0.008

    cp_row = np.concatenate(([lead_pt, lead_pt2],camber_pts,[trail_pt2, trail_pt]),axis=0)
    print("len(cp_row)")
    print(len(cp_row))
    print("cp_row")
    print(cp_row)
    # pyplot.scatter(cp_row[:,0],cp_row[:,1])
    # pyplot.show()



    u_coords = cp_row[:,0]
    y_coords = cp_row[:,1]
    z_coords = cp_row[:,2]
    
    print("u_coords")    
    print(u_coords)

    v_translations = np.linspace(-0.4 * param_t, 0.4 * param_t, num=nCPsV)
    w_translations = np.linspace(-.05, .05, num=nCPsW)

    uvw_coords =[]

    for z in range(nCPsW):
        for y in range(nCPsV):
            for x in range(len(u_coords)):
                if z == 0 and y == 0:
                    print(u_coords[x])
                uvw_coords.append([ u_coords[x], y_coords[x]+v_translations[y], z_coords[x]+w_translations[z]])

    #print(uvw_coords)
    uvw_coords = np.array(uvw_coords)


    # fig = pyplot.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(uvw_coords[:,0], uvw_coords[:,1], uvw_coords[:,2])
    # pyplot.show()

    of_header = """
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2006                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "../constant/controlPoints";
    object      airfoilcpsBsplines0;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

    f = open(output_filename, 'w')
    f.write(of_header)
    f.write("\n\ncontrolPoints " + str(len(uvw_coords)) + "\n")
    f.write("(\n")

    form="f"
    for cp in uvw_coords:
        f.write("( " + format(cp[0], form) + " " + format(cp[1], form) + " " + format(cp[2], form) + " )\n")

    f.write(");\n")
    f.close()

    print("Control points file generated.")
