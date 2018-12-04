'''
cell parser and synapse mapper 
for single cell simulations
with NeuroNet subcellular synapse distributions
'''

from writer import write_cell_simulation
from writer import write_landmark_file
from writer import write_cell_synapse_locations
from writer import write_synapse_activation_file
from writer import write_synapse_weight_file
from writer import write_sim_results
from writer import write_all_traces
from writer import write_PSTH
from writer import write_presynaptic_spike_times
from writer import write_spike_times_file
from reader import read_scalar_field
from reader import read_synapse_realization
from reader import read_synapse_activation_file
from reader import read_complete_synapse_activation_file
from reader import read_spike_times_file
from reader import read_synapse_weight_file
from reader import read_landmark_file
from synapse_mapper import SynapseMapper
from cell import Cell, PySection, PointCell
from cell import SynParameterChanger
from cell_parser import CellParser
#from synapse import activate_functional_synapse
from network import NetworkMapper
from network_realizations import create_synapse_realization
from network_realizations import create_functional_network
#from sim_control import SimControl
import neuron
from sumatra.parameters import build_parameters
from sumatra.parameters import NTParameterSet
import numpy as np

#------------------------------------------------------------------------------ 
# commonly used functions required for running single neuron simulations
#------------------------------------------------------------------------------ 
def load_NMODL_parameters(parameters):
    '''
    automatically loads NMODL mechanisms from paths in parameter file
    '''
    for mech in parameters.NMODL_mechanisms.values():
        neuron.load_mechanisms(mech)
    try:
        for mech in parameters.mech_globals.keys():
                for param in parameters.mech_globals[mech]:
                    paramStr = param + '_' + mech + '='
                    paramStr += str(parameters.mech_globals[mech][param])
                    print 'Setting global parameter', paramStr
                    neuron.h(paramStr)
    except AttributeError:
        pass

def create_cell(parameters, scaleFunc=None, allPoints=False):
    '''
    default way of creating NEURON cell models;
    includes spatial discretization and inserts
    biophysical mechanisms according to parameter file
    '''
    print '-------------------------------'
    print 'Starting setup of cell model...'
    axon = False
    if 'AIS' in parameters.keys():
        axon = True
    print 'Loading cell morphology...'
    parser = CellParser(parameters.filename)
    parser.spatialgraph_to_cell(axon, scaleFunc)
    print 'Setting up biophysical model...'
    parser.set_up_biophysics(parameters, allPoints)
    print '-------------------------------'
    return parser.cell

def init_neuron_run(simparam, vardt=False, *events):
    '''
    Default NEURON run with inital parameters
    according to parameter file.
    Optional parameters: callable "events" that are
    passed to Event objects holding a FInitializeHandler.
    This can be used to implement changes of parameters during
    the course of the simulation using h.cvode.event(t, "statement")
    in the supplied callable, where "statement" is another
    Python callable which may be used to change parameters.
    '''
#    use fixed time step for now
    neuron.h.load_file('stdrun.hoc')
    if vardt:
        cvode = neuron.h.CVode()
        cvode.active(1)
        # minimum tolerance: heuristically
        # tested with BAC firing
        # to give good tradeoff accuracy/speed
#        cvode.atol(1e-2)
#        cvode.rtol(2e-3)
#    neuron.h('using_cvode_=1')
#    neuron.h('cvode_active(1)')
#    cvode.use_local_dt(1)
#    cvode.condition_order(2)
#    cvode.atol(1e-3)
#    cvode.rtol(1e-12)
    eventList = []
    for event in events:
        e = Event(event)
        eventList.append(e)
#        print 'added cvode event to EventList'
    neuron.h.dt = simparam.dt
    neuron.h.celsius = simparam.T
    vInitStr = 'v_init=' + str(simparam.Vinit)
    neuron.h(vInitStr)
    neuron.h('init()')
#    neuron.h('run()')
#    neuron.h.finitialize(simparam.Vinit)
    neuron.run(simparam.tStop)

def sec_distance_to_soma(currentSec):
    '''compute path length from sec(x=0) to soma'''
    parentSec = currentSec.parent
    dist = 0.0
    parentLabel = parentSec.label
    while parentLabel != 'Soma':
        dist += parentSec.L
        currentSec = parentSec
        parentSec = currentSec.parent
        parentLabel = parentSec.label
    return dist

class Event():
    def __init__(self, func):
        self.callback = func
        self.fih = neuron.h.FInitializeHandler(1, self.callback)


