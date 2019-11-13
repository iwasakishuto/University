#coding: utf-8
import nest
import time
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# reset kernel for new example
class neuron():
    def __init__(self, model_name, multimeter=True, spike_detector=True,
                 label=None, color=None):
        self.neuron         = nest.Create(model_name)
        self.multimeter     = nest.Create("multimeter", params={"withtime": True, "record_from":["V_m"]})
        self.spike_detector = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
        self.label = label
        self.color = color
        self.built = False
    def set_params(self, node, **kwargs):
        nest.SetStatus(self.__dict__[node], kwargs)
    def build(self, display=True, **kwargs):
        if self.built: print("Requirement already satisfied.")
        else:
            nest.Connect(self.multimeter, self.neuron)
            nest.Connect(self.neuron, self.spike_detector)
            self.built=True

class NeuronA(neuron):
    def __init__(self, model_name="iaf_psc_delta"):
        super().__init__(model_name=model_name, label="A", color="tab:red")
        self.set_params("neuron", I_e=376.0)

class NeuronB(neuron):
    def __init__(self, model_name="iaf_psc_delta"):
        super().__init__(model_name=model_name, label="B", color="tab:blue")
        self.set_params("neuron", I_e=0.0)

def plotMultipleNeurons(*neurons, T=300):
    nest.Simulate(float(T))

    fig, ax = plt.subplots(2, 1, sharex=True, sharey=False)
    y_ticks=()
    for neuron in neurons:
        #=== multimeter ===
        multimeter_readout = nest.GetStatus(neuron.multimeter)[0]
        V = multimeter_readout["events"]["V_m"]
        t = multimeter_readout["events"]["times"]
        #=== spike detector ===
        spike_detector_readout = nest.GetStatus(neuron.spike_detector, keys="events")[0]
        event = spike_detector_readout["senders"]
        te    = spike_detector_readout["times"]

        ax[0].plot(t, V, color=neuron.color, label=neuron.label)
        ax[1].plot(te, event, ".", color=neuron.color)
        y_ticks+=neuron.neuron
    ax[0].set_ylabel('V (mV)')
    ax[0].legend()
    ax[0].grid(linestyle='-', linewidth=.25, alpha=.7)
    ax[1].set_ylim(0,5)
    ax[1].set_yticks(list(y_ticks))
    ax[1].set_yticklabels([neuron.label for neuron in neurons])
    ax[1].set_xlabel('time (ms)'), ax[1].set_ylabel('Spike times')
    ax[1].grid(linestyle='-', linewidth=.25, alpha=.7)
    return ax

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
Default_Neuron_Params = {"C_m": 1.0,
                         "tau_m": 20.0,
                         "t_ref": 2.0,
                         "E_L": 0.0,
                         "V_reset": 0.0,
                         "V_m": 0.0,
                         "V_th": 20.0}

class BrunelNodes():
    syn_spec=None
    def __init__(self, N_neurons, N_rec_volt, N_rec, noise):
        """Arguments
        @param N_neurons  : (int)
        @param N_rec_volt : (int) the amount of neurons whose voltages are recorded.
        @param N_rec,     : (int)
        @param noise      : (tuple)
        """
        self.neuron          = nest.Create("iaf_psc_delta" , N_neurons)
        self.spikes_detector = nest.Create("spike_detector", params={"withtime": True, "withgid": True, "to_file": False})
        self.voltmeter       = nest.Create("voltmeter"     , N_rec_volt)
        nest.Connect(noise, self.neuron, syn_spec=BrunelNodes.syn_spec)
        nest.Connect(self.neuron[:N_rec], self.spikes_detector, syn_spec=BrunelNodes.syn_spec)
        for i in range(N_rec_volt):
            nest.Connect([self.voltmeter[i]], [self.neuron[i]])

class BrunelNetwork():
    """Random balanced network (delta synapses)
    ----------------------------------------------
    This script simulates an excitatory and an inhibitory population on
    the basis of the network used in [1].

    > When connecting the network customary synapse models are used, which allow
      for querying the number of created synapses.
    > Using spike detectors the average firing rates of the neurons in the
      populations are established.
    > The building, as well as the simulation time of the network, are recorded.

    References
    ~~~~~~~~~~
    .. [Github] (Skill Pills +): https://github.com/oist-ncbc/skill-pill-plus/
    .. [1] Brunel N (2000). Dynamics of sparsely connected networks of
           excitatory and inhibitory spiking neurons.
           Journal of Computational Neuroscience 8, 183-208.
    """
    def __init__(self, g=5., eta=2., delay=1.5, epsilon=0.1, dt=0.1):
        """Arguments
        @param g      : (float) ratio of inhibitory / excitatory weight
        @param eta    : (float) external noise rate relative to the threshold rate
        @param delay  : (float) synaptic delay [ms]
        @param dt     : (float) the resolution [ms]
        @param epsilon: (float) connection probability
        """
        self.g = g
        self.eta = eta
        self.delay = delay
        self.dt = dt
        self.epsilon = epsilon

    def calcurate_Poisson_params(self, J, tauMem, theta, ex_synapses_per_neuron):
        self.J_ex = J                    # amplitude of excitatory postsynaptic potential
        self.J_in = -self.g * self.J_ex  # amplitude of inhibitory postsynaptic potential
        nu_th = theta / (J * ex_synapses_per_neuron * tauMem)
        nu_ex = self.eta * nu_th
        p_rate = 1000.0 * nu_ex * ex_synapses_per_neuron
        return float(p_rate)

    def build(self, N_total_neuron=12500, N_rec_volt=3, N_rec=50, J=0.1, **neuron_params):
        """Arguments
        @param N_total_neuron : (int) the amount population
        @param N_rec_volt     : (int) the amount of neurons whose voltages are recorded.
        @param N_rec          : (int)
        @param J              : (float) postsynaptic amplitude in mV.
        """
        self.n_rec_volt = N_rec_volt
        self.n_rec = N_rec
        if len(neuron_params)==0:
            print(f"Using default neuron parameters for `iaf_psc_delta` model. If you want to know more, please look at 'https://github.com/iwasakishuto/University/blob/gh-pages/3A/theme/program/OIST/tutorial3.py#L67'")
            neuron_params = Default_Neuron_Params

        #=== Initialization of the parameters =========
        nest.ResetKernel()
        start = time.time()
        n_excitatory = int(N_total_neuron*0.8)
        n_inhibitory = N_total_neuron-n_excitatory
        neuron_digit = len(str(N_total_neuron))
        print(f"Number of neurons : {N_total_neuron:>{neuron_digit}}")
        print(f"        Excitatory: {n_excitatory:>{neuron_digit}}")
        print(f"        Inhibitory: {n_inhibitory:>{neuron_digit}}")
        #----------------------------------------------
        ex_synapses_per_neuron  = int(self.epsilon * n_excitatory) # num of excitatory synapses per neuron.
        in_synapses_per_neuron  = int(self.epsilon * n_inhibitory) # num of inhibitory synapses per neuron.
        tot_synapses_per_neuron = ex_synapses_per_neuron + in_synapses_per_neuron # total number of synapse per neuron.
        synapses_digit = len(str(int(tot_synapses_per_neuron * N_total_neuron)))
        print(f"Number of synapses: {int(tot_synapses_per_neuron * N_total_neuron):>{synapses_digit}}")
        print(f"        Excitatory: {int(ex_synapses_per_neuron * N_total_neuron):>{synapses_digit}}")
        print(f"        Inhibitory: {int(in_synapses_per_neuron * N_total_neuron):>{synapses_digit}}")
        #----------------------------------------------
        tauMem = neuron_params.get("tau_m") # Time constant of membrane potential in ms
        theta = neuron_params.get("V_th")   # Membrane threshold potential in mV
        p_rate = self.calcurate_Poisson_params(J,tauMem,theta,ex_synapses_per_neuron)

        ###############################################
        #=== Configuration of the simulation kernel ===
        nest.SetKernelStatus({"resolution": self.dt, "print_time": True, "overwrite_files": True})
        nest.SetDefaults("iaf_psc_delta", neuron_params)
        nest.CopyModel("static_synapse", "excitatory", {"weight": self.J_ex, "delay": self.delay})
        nest.CopyModel("static_synapse", "inhibitory", {"weight": self.J_in, "delay": self.delay})
        BrunelNodes.syn_spec = "excitatory"

        poisson_noise = nest.Create("poisson_generator", params={"rate": p_rate})
        nodes_ex = BrunelNodes(n_excitatory, N_rec_volt, N_rec, poisson_noise)
        nodes_in = BrunelNodes(n_inhibitory, N_rec_volt, N_rec, poisson_noise)

        nest.Connect(nodes_ex.neuron,
                     nodes_ex.neuron + nodes_in.neuron,
                     {'rule': 'fixed_indegree', 'indegree': ex_synapses_per_neuron},
                     "excitatory")
        nest.Connect(nodes_in.neuron,
                     nodes_ex.neuron + nodes_in.neuron,
                     {'rule': 'fixed_indegree', 'indegree': in_synapses_per_neuron},
                     "inhibitory")
        self.nodes_ex = nodes_ex
        self.nodes_in = nodes_in
        print(f"Building time     : {time.time()-start:.2f} s")

    def simulate(self, T=300., memorize=True):
        self.T = T
        start = time.time()
        nest.Simulate(float(T))
        print(f"Simulation time   : {T} ms")
        print(f"Processing time   : {time.time()-start:.2f} s")

        events_ex = nest.GetStatus(self.nodes_ex.spikes_detector, "n_events")[0]
        events_in = nest.GetStatus(self.nodes_in.spikes_detector, "n_events")[0]
        rate_ex = events_ex / T * 1000.0 / self.n_rec
        rate_in = events_in / T * 1000.0 / self.n_rec
        # num_synapses = (nest.GetDefaults("excitatory")["num_connections"] +
        #                 nest.GetDefaults("inhibitory")["num_connections"])
        allevents_E = nest.GetStatus(self.nodes_ex.spikes_detector, "events")[0]
        allevents_I = nest.GetStatus(self.nodes_in.spikes_detector, "events")[0]

        memb_potential_ex = [nest.GetStatus([self.nodes_ex.voltmeter[i]])[0]["events"]["V_m"]   for i in range(self.n_rec_volt)]
        memb_times_ex     = [nest.GetStatus([self.nodes_ex.voltmeter[i]])[0]["events"]["times"] for i in range(self.n_rec_volt)]
        memb_potential_in = [nest.GetStatus([self.nodes_in.voltmeter[i]])[0]["events"]["V_m"]   for i in range(self.n_rec_volt)]
        memb_times_in     = [nest.GetStatus([self.nodes_in.voltmeter[i]])[0]["events"]["times"] for i in range(self.n_rec_volt)]

        if memorize:
            self.allevents_E = allevents_E
            self.allevents_I = allevents_I
            self.memb_potential_ex = memb_potential_ex
            self.memb_times_ex = memb_times_ex
            self.memb_potential_in = memb_potential_in
            self.memb_times_in = memb_times_in
            print("Completed. Please use following methods to plot results.",
                  "\n- `plotVoltageEx(n)`, `plotVoltageIn(n)`",
                  "\n- `plotVoltageCompare(idx=0)`",
                  "\n- `plotRaster()`")

    def plotVoltageEx(self, n=None):
        n = self.n_rec_volt if n is None else n
        print("\033[07m\033 Membrane potential plots [Excitatory Neurons] \033[0m")
        self._plotVoltage(n, self.memb_potential_ex[:n], self.memb_times_ex[:n])

    def plotVoltageIn(self, n=None):
        n = self.n_rec_volt if n is None else n
        print("\033[07m\033 Membrane potential plots [Inhibitory Neurons] \033[0m")
        self._plotVoltage(n, self.memb_potential_in[:n], self.memb_times_in[:n])

    def _plotVoltage(self, n, memb_potential, times):
        fig = plt.figure(figsize=(10,2*n))
        for i, (t_E,v_E) in enumerate(zip(times, memb_potential)):
            ax = fig.add_subplot(n, 1, i+1)
            ax.axhline(0, color='tab:grey', linestyle='--'), ax.axhline(20, color='tab:grey', linestyle='--')
            ax.text(305, -2, "V$_{reset}$", color='tab:grey', fontsize=12), ax.text(305, 18, "V$_{th}$",    color='tab:grey', fontsize=12)
            plt.plot(t_E, v_E, color=cm.hsv(float(i) / n))
            ax.set_xlim(0, 300), ax.set_ylim(-5, 25), ax.set_ylabel('V$_{E,0}$ (mV)')
            if i+1<n: plt.setp(ax.get_xticklabels(), visible=False)
        plt.tight_layout()
        plt.show()

    def plotVoltageCompare(self, idx=0):
        print("\033[07m\033 Membrane potential Compaison. \033[0m")
        if idx>=self.n_rec_volt:
            raise IndexError(f"list index out of range. Please specify idx < {self.n_rec_volt}")
        fig = plt.figure(figsize=(10,6))
        lst = [
            (self.memb_times_ex[idx], self.memb_potential_ex[idx], "Excitatory", "blue"),
            (self.memb_times_in[idx], self.memb_potential_in[idx], "Inhibitory", "red"),
        ]
        for i,(t,v,title,color) in enumerate(lst):
            ax = fig.add_subplot(2,1,i+1)
            ax.axhline(0, color='tab:grey', linestyle='--'), ax.axhline(20, color='tab:grey', linestyle='--')
            ax.plot(t, v, color=color, label="")
            ax.set_xlim(0,self.T), ax.set_ylim(-5, 25), ax.set_ylabel('V$_{E,0}$ (mV)'), ax.set_title(title)
            ax.text(self.T+10, -2, "V$_{reset}$", color='tab:grey', fontsize=12), ax.text(self.T+10, 18, "V$_{th}$", color='tab:grey', fontsize=12)
            if i==0: plt.setp(ax.get_xticklabels(), visible=False)
        plt.tight_layout()
        plt.show()

    def plotRaster(self):
        print(f"\033[07m\033 Raster Plot.(N={self.n_rec}) \033[0m")
        fig, ax = plt.subplots(2, 2, sharex=True, sharey=False, figsize=(14,6))
        lst = [
            (self.allevents_E, "Excitatory", "blue"),
            (self.allevents_I, "Inhibitory", "red"),
        ]
        for i,(e,types,color) in enumerate(lst):
            ax[0,i].set_title(f"{types} Population"), ax[0,i].set_xlim(0, self.T), ax[0,i].set_ylabel('Neuron ID')
            ax[0,i].plot(e['times'], e['senders'], '.', color=color, markersize=2)
            ax[1,i].set_xlim(0, self.T), ax[1,i].set_xlabel('time (ms)'), ax[1,i].set_ylabel('Frequency (Hz)')
            count, bins = np.histogram(e['times'], bins=100)
            binwidth = bins[1]-bins[0]
            ax[1,i].plot(bins[:-1], count*1000./self.n_rec/binwidth, color=color)
        plt.tight_layout()
        plt.show()
