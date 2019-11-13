#coding: utf-8
import nest
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
def poisson_train_generator(simulation_length, num_presynaptic_neurons, min_span,
                            starting_Hz, min_r, max_r,
                            min_ds, max_ds, min_dr, max_dr, **kwargs):
    spikes = np.zeros(shape=(num_presynaptic_neurons,simulation_length))
    spikes[:,0] = np.clip(np.random.poisson(starting_Hz/1000), 0, 1)

    # ds is randomly picked from a uniform distribution over [−360,+360] Hz/s
    ds = np.random.uniform(min_ds, max_ds, size=(num_presynaptic_neurons,simulation_length))

    dr = np.zeros_like(ds)
    for t in range(simulation_length-1):
        # dr is clipped to within [−1800,1800] (Hz/s)
        dr[:,t+1] = np.clip(dr[:,t] + ds[:,t], min_dr, max_dr)

    # r varies beteen [0,90] (Hz)
    r = np.insert(dr[:,:-1], 0, starting_Hz, axis=1)
    for t in range(simulation_length-1):
        r[:,t+1] = np.clip(r[:,t] + dr[:,t], min_r, max_r)

        # Manually add some additional spikes to guarantee that
        # "in every 50ms time bin, each neuron spikes at least once".
        spikes[:,t+1] = np.clip(np.random.poisson(r[:,t+1]/1000), 0, 1)
        if t+1>=min_span:
            mask    = np.where(np.sum(spikes[:,t+1-min_span:t+1], axis=1)==0)[0]
            add_pos = np.random.randint(t+1-min_span, t+1, size=len(mask))
            for idx,pos in zip(mask,add_pos):
                spikes[idx,pos] = 1
    return spikes
