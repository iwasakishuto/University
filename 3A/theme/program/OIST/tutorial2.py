#coding: utf-8
import nest
import numpy as np
import matplotlib.pyplot as plt

Abbre2Formal = {
    "hh"   : "Hodgkin-Huxley",
    "iaf"  : "LIF",
    "alpha": "alpha-function shaped",
    "delta": "delta-shaped",
}

def symbol2units(symbol):
    initial = symbol[0]
    if initial == "C": return "pF"
    elif initial in ["V","E"]: return "mV"
    elif initial == "t": return "ms"
    elif initial == "g": return "nS"
    elif initial == "I": return "pA"
    else: return ""

def model2verbose(name):
    network, _, func = name.split("_")
    return f"Create {Abbre2Formal[network]} neuron with {Abbre2Formal[func]} synaptic currents."

class SingleNeuralModel():
    def __init__(self, model_name, spike_generator="spike_generator",
                 voltmeter=True, spike_detector=True, display=True, reset=True):
        if reset: nest.ResetKernel()
        #=== These Nodes only return the tuple of the ids. ===
        self.neuron          = nest.Create(model_name) # id=(1,)
        self.spike_generator = nest.Create(spike_generator)  if spike_generator else None
        self.voltmeter       = nest.Create('voltmeter')      if voltmeter       else None
        self.spike_detector  = nest.Create('spike_detector') if spike_detector  else None
        if display:
            print(model2verbose(model_name))
            if spike_generator: print(f"Created a {spike_generator}.")
            if voltmeter:       print("Created a voltmeter.")
            if spike_detector:  print("Created a spike detector.")
        #=====================================================
        self.built = False

    def set_params(self, node, **kwargs):
        nest.SetStatus(self.__dict__[node], kwargs)

    def build(self, display=True, **kwargs):
        """ Connect the spike generator and voltmeter to the neuron. """
        if self.built: print("Requirement already satisfied.")
        else:
            # spike_generator to Neuron.
            if self.spike_generator is not None:
                nest.Connect(self.spike_generator, self.neuron, syn_spec=kwargs)
                if display: print(f"Connect spike generator with a given synaptic specification ({kwargs})")
            # Voltmeter to Neuron.
            if self.voltmeter is not None:
                nest.Connect(self.voltmeter, self.neuron)
                if display: print("Connected voltmeter to the neuron for measurements.")
            # Neuron to Spike detector.
            if self.spike_detector is not None:
                nest.Connect(self.neuron, self.spike_detector)
                if display: print("Connected the neuron to spike detector for measurements.")
            self.built=True

    def spike_params(self, *keys):
        if keys: return nest.GetStatus(self.spike_generator, keys)
        else: return nest.GetStatus(self.spike_generator)

    def voltage_params(self, *key):
        if keys: return nest.GetStatus(self.voltmeter, keys)
        else: return nest.GetStatus(self.voltmeter)

    def simulate(self, ms=100., ax=None, plot=True, **plotkwargs):
        if ax==None and plot==True: fig, ax = plt.subplots()
        nest.Simulate(float(ms))
        # Read out recording time and voltage from voltmeter (check the parameter list!)
        self.spikes = nest.GetStatus(self.spike_detector, "n_events")[0]
        if self.voltmeter is not None:
            results = nest.GetStatus(self.voltmeter)[0]['events']
            times = results['times']
            voltage = results['V_m']
            #=== Plot the Results. ===
            if plot:
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
    def __init__(self, model_name='iaf_psc_alpha', spike_generator="spike_generator", display=True, reset=True):
        super().__init__(model_name=model_name,
                         spike_generator=spike_generator,
                         display=display,
                         reset=reset)
        self.neuron_disparams = ["C_m", "E_L", "I_e", "V_m", "V_reset", "V_th", "t_ref", "tau_m",]

class HHmodel(SingleNeuralModel):
    def __init__(self, model_name="hh_psc_alpha", spike_generator="spike_generator", display=True, reset=True):
        super().__init__(model_name=model_name,
                         spike_generator=spike_generator,
                         display=display,
                         reset=reset)
        self.neuron_disparams = ["C_m", "E_K", "E_L", "E_Na", "g_K", "g_L", "g_Na", "I_e", "V_m", "t_ref"]

def mkMultiDetectors(num, display=True):
    multimeters = nest.Create('multimeter', num)
    nest.SetStatus(multimeters, {'record_from': ['V_m']})
    if display: print(f"Created {num} multimeter to record membrane potential of the neurons.")
    spike_detectors = nest.Create('spike_detector', num)
    nest.SetStatus(spike_detectors, [{'withtime': True,'withgid': True,'to_file': False}])
    if display: print(f"Created {num} spike detectors to record spikes times and neuron identifiers, but not record from file.")
    return multimeters,spike_detectors

def plotMultResults(multimeters, spike_detectors, rate=False, **plotkwargs):
    num = len(multimeters)
    if num != len(spike_detectors): raise ValueError(f"'multimeters' and 'spike detectors' should be the same lengths ({num}!={len(spike_detectors)})")

    fig = plt.figure(figsize=(14,4))
    axL = fig.add_subplot(1,2,1)
    axR = fig.add_subplot(1,2,2)
    for i in range(num):
        data   = nest.GetStatus([multimeters[i]])[0]['events']
        V_mem  = data["V_m"]
        times  = data["times"]
        spikes = nest.GetStatus([spike_detectors[i]])[0]['events']['times']
        plot_info = dict([(k,vals[i]) for k,vals in plotkwargs.items()])
        if rate:
            plot_info["label"] += f" (Rate of neuron stimulated: {float(len(spikes))/max(times)*1e3:.3f})"
        else:
            print(f"Rate of neuron stimulated with {plot_info['label']} input: {float(len(spikes))/max(times)*1e3:.3f}")
        hist, _ = np.histogram(V_mem)
        if max(hist)/len(V_mem) < 0.8: axL.hist(V_mem, 100, alpha=0.7, **plot_info)
        axR.plot(times, V_mem, **plot_info)
    #=== Design ===
    axL.set_xlabel(r'$V_m$ (mV)'), axL.grid(linestyle='-', linewidth=.25, alpha=.7)
    axR.set_xlabel('Time (ms)'), axR.set_ylabel(r'$V_m$ (mV)')
    axR.set_xlim([0., max(times)]), axR.set_ylim([-5., 20.])
    axR.grid(linestyle='-', linewidth=.25, alpha=.7)
    return (axL, axR)

def FiringRate(model_name, I_e, T_sim=1e3, **kwargs):
    model = SingleNeuralModel(model_name=model_name, voltmeter=False, display=False)
    model.set_params("neuron", **kwargs)
    model.set_params("neuron", I_e=I_e)
    model.set_params("spike_detector", withtime=True, withgid=True, to_file=False)
    model.build(display=False)
    model.simulate(T_sim, plot=False)
    return model.spikes
