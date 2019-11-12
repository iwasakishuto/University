#coding: utf-8
import nest
import matplotlib.pyplot as plt

def symbol2units(symbol):
    initial = symbol[0]
    if initial == "C": return "pF"
    elif initial in ["V","E"]: return "mV"
    elif initial == "t": return "ms"
    elif initial == "g": return "nS"
    elif initial == "I": return "pA"
    else: return ""

class SingleNeuralModel():
    def __init__(self, model_name, neuron_verbose="", display=True):
        self.display = display
        nest.ResetKernel() # reset simulation kernel
        #=== These Nodes only return the tuple of the ids. ===
        self.neuron = nest.Create(model_name) # id=(1,)
        if self.display: print(neuron_verbose)
        self.spikegenerator = nest.Create('spike_generator') # id=(2,)
        if self.display: print("Created a spike generator.")
        self.voltmeter = nest.Create('voltmeter') # id=(3,)
        if self.display: print("Created a voltmeter.")
        self.spikedetector = nest.Create('spike_detector')
        if self.display: print("Created a spike detector.")
        #=====================================================
        self.built = False

    def set_params(self, node, **kwargs):
        nest.SetStatus(self.__dict__[node], kwargs)

    def build(self, weight=1e3):
        """ Connect the spike generator and voltmeter to the neuron. """
        if self.built:
            print("Requirement already satisfied.")
        else:
            nest.Connect(self.spikegenerator, self.neuron, syn_spec={"weight": weight})
            if self.display: print(f"Connect spike generator with a given synaptic specification (weight={weight})")
            nest.Connect(self.voltmeter, self.neuron)
            if self.display: print("Connected voltmeter to the neuron for measurements.")
            nest.Connect(self.neuron, self.spikedetector)
            if self.display: print("Connected the neuron to spikedetector for measurements.")
            self.built=True

    def spike_params(self, *keys):
        if keys: return nest.GetStatus(self.spikegenerator, keys)
        else: return nest.GetStatus(self.spikegenerator)

    def voltage_params(self, *key):
        if keys: return nest.GetStatus(self.voltmeter, keys)
        else: return nest.GetStatus(self.voltmeter)

    def simulate(self, ms=100., ax=None, **plotkwargs):
        if ax==None: fig, ax = plt.subplots()
        nest.Simulate(float(ms))
        # Read out recording time and voltage from voltmeter (check the parameter list!)
        results = nest.GetStatus(self.voltmeter)[0]['events']
        self.spikes = nest.GetStatus(self.spikedetector, "n_events")[0]
        times = results['times']
        voltage = results['V_m']
        #=== Plot the Results. ===
        ax.plot(times,voltage,**plotkwargs)
        ax.set_xlabel('Time (ms)'), ax.set_ylabel('Membrance potential (mV)'), ax.grid(linestyle='-', linewidth=.25, alpha=.7)
        return ax

    def neuron_params(self, *keys):
        if keys:
            # retrieve a particular set of parameters
            return nest.GetStatus(self.neuron, keys)
        else:
            # get the parameter list and values of a node
            parameters = nest.GetStatus(self.neuron)[0]
            width = max([len(key) for key in self.neuron_disparams])
            for key in self.neuron_disparams:
                print(f"{key:<{width}}: {parameters[key]} {symbol2units(key)}")


class LIFmodel(SingleNeuralModel):
    def __init__(self, display=True):
        super().__init__(model_name='iaf_psc_alpha',
                         neuron_verbose="Created LIF neuron with alpha-function shaped synaptic currents.",
                         display=display)
        self.neuron_disparams = ["C_m", "E_L", "tau_m", "V_m", "V_reset", "V_th", "I_e"]


class HHmodel(SingleNeuralModel):
    def __init__(self, display=True):
        super().__init__(model_name="hh_psc_alpha",
                         neuron_verbose="Create Hodgkin-Huxley neuron with delta-shaped synaptic currents.",
                         display=display)
        self.neuron_disparams = ["C_m", "E_K", "E_L", "E_Na", "g_K", "g_L", "g_Na", "I_e", "V_m", "t_ref"]
