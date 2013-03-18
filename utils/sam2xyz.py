#!pvpython

# export three .xyz files for displacement on a horizontal slice out of gamra HDF5 output

import getopt
import glob
import os
import sys
from paraview.simple import *

def usage():
    print 'sam2xyz.py exports a slice of samrai hdf5 export file to .xyz format'
    print ''
    print 'usage: sam2xyz.py --depth=0 summary.samrai'
    print ''
    print 'options:'
    print '  -d --depth:  the depth of the horizontal slice'
    print ''
    print 'example:'
    print './sam2xyz.py --depth=0 ./Landers.visit/visit_dump.00000/summary.samrai'
    print ''
    print 'or'
    print ''
    print './sam2xyz.py --depth=0 */*/summary.samrai'
    print ''


def main():

    # default parameters
    depth = 0

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:", ["help","depth="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print >> sys.stderr, 'sam2xyz.py:', str(err) # will print something like "option -a not recognized"
        print ''
        usage()
        sys.exit(2)

    if 0 == len(args):
	usage()
	sys.exit(2)

    for o, a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-d", "--depth"):
            depth = float(a)
        else:
            print >> sys.stderr, 'sam2xyz.py: unhandled option:', o, a
            assert False, "unhandled option"

    print '# sam2xyz.py '+" ".join(sys.argv[1:])

    # loop over the models
    for i in xrange(len(args)):
        ifile = os.path.abspath(args[i])
        ofile = os.path.abspath(os.path.dirname(args[i])+'/.'+os.path.basename(args[i])+'.csv')

        summary_samrai = VisItSAMRAIReader( FileName = ifile )
        summary_samrai.CellArrays = []
        summary_samrai.Materials = []
        summary_samrai.Meshes = ['amr_mesh']
        summary_samrai.PointArrays = []
        summary_samrai.CellArrays = ['Displacement']

        CellDatatoPointData2 = CellDatatoPointData()

        Slice2 = Slice( SliceType="Plane" )
        Slice2.SliceType.Origin = [-20.0, -20.0, depth+0.011]
        Slice2.SliceType = "Plane"
        Slice2.SliceType.Normal = [0,0,1]
        SetActiveSource(Slice2)
        UpdatePipeline()
        writer = CreateWriter(ofile, Slice2)
        writer.FieldAssociation = "Points" # or "Cells"
        writer.UpdatePipeline()
        del writer

        ofiles = glob.glob(os.path.dirname(ifile)+'/.'+os.path.basename(ifile)+'*.csv')

        fname_north = os.path.abspath(args[i]+'-north.xyz')
        fname_east  = os.path.abspath(args[i]+'-east.xyz')
        fname_up    = os.path.abspath(args[i]+'-up.xyz')

        of_north = open(fname_north, 'wb')
        of_east  = open(fname_east , 'wb')
        of_up    = open(fname_up   , 'wb')
        #of_north.write('> x y u1(north)\n')
        #of_east.write( '> x y u2(east)\n')
        #of_up.write(   '> x y u3(down)\n')

        for fname in ofiles:
            f = open(fname,'r')

            # "Displacement:0","Displacement:1","Displacement:2","Points:0","Points:1","Points:2"
            header = f.readline()

            for line in f:
                line = line.strip()
                cols = map(float,line.split(','))

                if len(cols)>3:
                    of_north.write("%+05.3e %+05.3e %+05.3e\n" % (cols[4],cols[3],cols[0]))
                    of_east.write("%+05.3e %+05.3e %+05.3e\n" % (cols[4],cols[3],cols[1]))
                    of_up.write("%+05.3e %+05.3e %+05.3e\n" % (cols[4],cols[3],-cols[2]))
            f.close()

            os.remove(fname)

        of_north.close()
        of_east.close()
        of_up.close()

if __name__ == "__main__":
    main()


