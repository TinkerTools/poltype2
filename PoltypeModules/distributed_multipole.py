from collections import OrderedDict
import sys

# This class handles I/O of GDMA input file
# It is designed to be flexible but not very smart.
# The keywords are stored in OrderedDict.
# The Radius can be added/updated by any of the following methods:
#     inp = GdmaInput(Radius_H=0.35, Radius_O=0.65)
#     inp = GdmaInput(Radius_H=0.35); inp.update(Radius_O=0.65)
#     inp = GdmaInput(Radius_H=0.35); inp.update({'Radius_O':0.65})
#     inp = GdmaInput(); inp.update({'Radius H':0.35, 'Radius O':0.65})
#     inp = GdmaInput(); inp.update([('Radius H',0.35), ('Radius O',0.65)])
# The following way also works but is not recommended
#     inp = GdmaInput(Radius="H 0.35")
# And below is a wrong example, as the option for Radius H will be lost
#     inp = GdmaInput(Radius="H 0.35"); inp.update(Radius="O 0.65")
class GdmaInput:
    """ Store, update and write GDMA options
    """
    def __init__(self, **kwargs):
        self._header = OrderedDict({"Title":"poltype gdma", "File":"dma.fchk density MP2", "Names":None, "Angstrom":"", "AU":""})
        self._footer = OrderedDict({"":"", "Start":"", "Finish":""})
        self._multipole = OrderedDict({"Multipoles":"", "Switch":0, "Limit":2, "Punch":"dma.punch"})

        self.update(**kwargs)

    def update(self, keywords=None, **kwargs):
        '''Update GDMA keywords
        '''
        if keywords is None:
            keywords = {}
        kws = OrderedDict(keywords)
        kws.update(kwargs)
        for kw0, val in kws.items():
            kw = ' '.join(kw0.replace('_', ' ').split())
            for _dict in (self._header, self._footer, self._multipole):
                if kw in _dict:
                    break
            _dict[kw] = val

    def write_file(self, outfile):
        '''Write GMDA input file
        '''
        _data = list(self._header.items()) + list(self._multipole.items()) + list(self._footer.items())
        outp = ''
        for kw, val in _data:
            if val is None:
                continue
            outp += ('%s %s'%(kw, val)).strip() + '\n'
        with open(outfile, 'w') as fh:
            fh.write(outp)

def get_dma_default(method='dma0'):
    if method == 'dma4':
        inp = GdmaInput(Switch=4, 
            #GDMA suggested
            Radius_H=0.325, Radius_C=0.65, Radius_N=0.65, Radius_O=0.65, Radius_F=0.65, Radius_Cl=1.0, 
            #Scaled based on Cl, vdw-radii
            Radius_Br=1.05, Radius_I=1.10, Radius_S=0.90, Radius_P=0.90)
    else:
        inp = GdmaInput(Switch=0, Radius_H=0.65, Radius_S=0.80, Radius_P=0.75, 
            Radius_Cl=1.00, Radius_Br=1.1, Radius_I=1.3)
    return inp

if __name__ == '__main__':
    method = 'dma0'
    if len(sys.argv) > 2:
        method = sys.argv[2]
    if len(sys.argv) > 1:
        outfile = sys.argv[1]
        inp = get_dma_default(method)
        inp.write_file(outfile)


