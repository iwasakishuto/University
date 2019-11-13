#coding: utf-8
import nest
import time
import random
import numpy as np
import matplotlib.pyplot as plt

class AddSTDP():
    """ Additive STDP rule """
    def __init__(self, **params):
        if len(params) > 0:
            self.load_params(**params)

    def load_params(self, **params):
        self.__dict__.update(params)

    def simulate(self, time_differences, ax=None):
        if ax is None: fig,ax = plt.subplots()
        weight_changes = [self._add_STDP(t_diff) for t_diff in time_differences]
        ax.plot(time_differences, weight_changes, color="black"), ax.grid(linestyle='-', linewidth=.25, alpha=.7)
        ax.set_title("Additive STDP rule curve"), ax.set_xlabel("$\Delta t\ (ms)$"), ax.set_ylabel("$\Delta w$")

    def _add_STDP(self, t_diff):
        return self.a_LTP * np.exp(t_diff/self.tau_LTP) if t_diff <= 0 else -self.a_LTD * np.exp(-t_diff/self.tau_LTD)

#-------------------------------------------------------------------------------
#---------------------------------- [ Step2 ] ----------------------------------
def plotPatternRaster(spikes, ax=None, display=True, cmap="binary_r", **kwargs):
    if ax is None: fig,ax=plt.subplots()
    if display: print("White plots mean Firing.")
    H,W = spikes.shape
    ax.matshow(spikes, aspect=W/H, cmap=cmap, **kwargs)
    ax.set_title("Raster plot of presynaptic neurons"), ax.set_xlabel("time (ms)"), ax.set_ylabel("neuron ID")
    return ax

def plotFiringRates(spikes, fig=None):
    if fig is None: fig = plt.figure()
    num_neurons, _ = spikes.shape
    firing_rates = np.mean(spikes, axis=1)*100
    neurons_IDs  = np.arange(num_neurons)
    #=== Plot ===
    ax1 = fig.add_subplot(2,1,1)
    ax1.barh(neurons_IDs, firing_rates, height=1)
    ax1.set_title("Average firing rates of presynaptic neurons"), ax1.set_ylabel("Neuron ID")

    ax2 = fig.add_subplot(2,1,2)
    ax2.hist(firing_rates)
    ax2.set_xlabel("Average firing rate (Hz)"), ax2.set_ylabel("Frequency")
    return fig

class PoissonNoiseGenerator():
    def __init__(self, **params):
        if len(params) > 0: self.load_params(**params)
        self.repeated=False
        self.pattern_repeated_spikes=None

    def load_params(self, **params):
        self.__dict__.update(params)

    def generate_spikes(self, simulation_length=1000, num_presynaptic_neurons=2000, seed=None):
        np.random.seed(seed)
        start = time.time()

        self.simulation_length = simulation_length
        self.num_presynaptic_neurons = num_presynaptic_neurons

        spikes = np.zeros(shape=(num_presynaptic_neurons,simulation_length))
        spikes[:,0] = np.clip(np.random.poisson(self.starting_Hz/1000), 0, 1)
        # ds is randomly picked from a uniform distribution over [−360,+360] Hz/s
        ds = np.random.uniform(self.min_ds, self.max_ds, size=(num_presynaptic_neurons, simulation_length))
        dr = np.zeros_like(ds)
        for t in range(simulation_length-1):
            # dr is clipped to within [−1800,1800] (Hz/s)
            dr[:,t+1] = np.clip(dr[:,t] + ds[:,t], self.min_dr, self.max_dr)
        # r varies beteen [0,90] (Hz)
        r = np.insert(dr[:,:-1], 0, self.starting_Hz, axis=1)
        for t in range(simulation_length-1):
            r[:,t+1] = np.clip(r[:,t] + dr[:,t], self.min_r, self.max_r)
            # Manually add some additional spikes to guarantee that
            # "in every 50ms time bin, each neuron spikes at least once".
            spikes[:,t+1] = np.clip(np.random.poisson(r[:,t+1]/1000), 0, 1)
            if t+1>=self.min_span:
                mask    = np.where(np.sum(spikes[:,t+1-self.min_span:t+1], axis=1)==0)[0]
                add_pos = np.random.randint(t+1-self.min_span, t+1, size=len(mask))
                for idx,pos in zip(mask,add_pos):
                    spikes[idx,pos] = 1
        self.spikes = spikes
        print(f"Processing Time: {time.time()-start:.3f}[s]")
        return spikes

    def select_segment(self, seed=None):
        self.num_segments = int(self.simulation_length/self.segment_length)
        selection = np.random.RandomState(seed).randint(1, self.num_segments)
        pattern_spikes = self.spikes[:,selection*self.segment_length:(selection+1)*self.segment_length]
        self.pattern_spikes = pattern_spikes
        self.selection = selection
        print(f"Spikes divided into {self.num_segments} segments.")
        print(f"{selection}-th segment was selected as a Pattern.")

    def repeat_segment(self, proportion_of_neurons=0.5, frequency_of_patterns=0.25):
        if self.repeated:
            print("Requirement already satisfied.")
            return (None,None)
        else:
            self.pattern_repeated_spikes = np.copy(self.spikes)
            self.repeated_position = np.zeros(self.spikes.shape)
            num_repeats = round(frequency_of_patterns*self.num_segments)
            idx_repeats = random.sample(range(self.num_segments), num_repeats)
            for idx in idx_repeats:
                self.pattern_repeated_spikes[0:int(self.num_presynaptic_neurons*proportion_of_neurons), idx*self.segment_length:(idx+1)*self.segment_length] = self.pattern_spikes[0:int(self.num_presynaptic_neurons*proportion_of_neurons),:]
                self.repeated_position[0:int(self.num_presynaptic_neurons*proportion_of_neurons),idx*self.segment_length:(idx+1)*self.segment_length] = 1
            print(f"Number of spikes added: {np.sum(self.repeated_position)}")
            print(f"Spikes are saved as `self.pattern_repeated_spikes`.")
            self.repeated = True
            return self.pattern_repeated_spikes, self.repeated_position

    def generatePoissonNoise(self, frequency=10):
        noise = np.random.poisson(frequency/1000, size=(self.spikes.shape))
        return noise

    def addPoisonNoise(self, frequency=10):
        spikes = self.pattern_repeated_spikes if self.pattern_repeated_spikes is not None else self.spikes
        noise = self.generatePoissonNoise(frequency=frequency)
        self.noise_added_spikes = np.clip(spikes + noise, 0, 1)
        print(f"Spikes are saved as `self.noise_added_spikes`.")
        return self.noise_added_spikes
