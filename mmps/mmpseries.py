from rdkit import Chem
import sys
from collections import namedtuple

Frag = namedtuple( 'Frag', [ 'id', 'scaffold', 'rgroup' ] )

class Series():
    def __init__( self ):
        self.rgroups = []
        self.scaffold = ""

def getFrags( filename ):
    frags = []
    for line in open( filename ):
        broken = line.rstrip().split( "," )
        if broken[ 2 ]: # single cut
            continue
        smiles = broken[ -1 ].split( "." )

        mols = [ Chem.MolFromSmiles( smi ) for smi in smiles ]
        numAtoms = [ mol.GetNumAtoms() for mol in mols ]

        if numAtoms[ 0 ] > 5 and numAtoms[ 1 ] < 12:
            frags.append( Frag( broken[1], smiles[0], smiles[1] ) )
        if numAtoms[ 1 ] > 5 and numAtoms[ 0 ] < 12:
            frags.append( Frag( broken[1], smiles[1], smiles[0] )  )
    frags.sort( key=lambda x:( x.scaffold, x.rgroup ) )
    return frags

def getSeries( frags ):
    oldfrag = Frag( None, None, None )
    series = Series()
    for frag in frags:
        if frag.scaffold != oldfrag.scaffold:
            if len( series.rgroups ) >= 2:
                series.scaffold = oldfrag.scaffold
                yield series
            series = Series()
        series.rgroups.append( ( frag.rgroup, frag.id ) )
        oldfrag = frag
    if len( series.rgroups ) >= 2:
        series.scaffold = oldfrag.scaffold
        yield series

if __name__ == "__main__":
    filename = sys.argv[1]

    frags = getFrags( filename )
    it = getSeries( frags )
    for series in it:

        print( "# %s" % series.scaffold )
        for rgroup in sorted( series.rgroups ):
            print( "%s %s" % ( rgroup[0], rgroup[1] ) )
