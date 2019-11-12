#coding: utf-8
import numpy as np

def LIFsimulation(V_rest=-65, V_thre=-40, V_reset=-70, V_fire=20, C=2e-4, R=7.5e4, I=5e-4, dt=0.1, T=100):
    """Leaky Integrate-and-Fire(LIF) model.
    I-(V-Vrest)/Rm = Cm*(dV/dt)
    """
    Vm = V_rest
    Time = np.arange(0,T+dt,dt)
    V = np.zeros_like(Time)
    for i in range(len(V)):
        V[i] = Vm
        # dV/dt=(IR-(Vm-Vrest))/RC
        dVdt = (I*R-(Vm-V_rest))/(R*C)
        V_plus_dV = Vm+dt*dVdt
        Vm = V_reset if Vm>V_thre else V_fire if V_plus_dV > V_thre else V_plus_dV
    return Time,V

def LITSpikeFrequency(Is, **params):
    freqs=[]
    for I in Is:
        params["I"]=I
        Time,V = LIFsimulation(**params)
        first, second = np.where(V==params.get("V_fire", 20))[0][:2]
        span = (second-first)*params.get("dt", 0.1)*1e-3
        freqs.append(1/span)
    return np.asarray(freqs)

def HHsimulation(V_rest=-65, V_thre=-55, V_reset=-70, V_fire=20,
                 gl=3.0e1, gK=3.6e3, gNa=1.2e4,
                 Vl=-54.402, VK=-77.0, VNa=50.0,
                 Cm=1e2, I=0, dt=0.1, T=100):
    """Hodgkin-Huxley model. (`hh_psc_alpha`)
    I=Cm*(dVm/dt) + gl*(Vm-Vl) + gK*n^4*()
    """
    alpha_n = lambda Vm: (1e-2*(Vm+55)) / (1-np.exp(-(Vm+55)/10))
    beta_n  = lambda Vm: 0.125*np.exp(-(Vm+65)/80)
    alpha_m = lambda Vm: (0.1*(Vm+40)) / (1-np.exp(-(Vm+40)/10))
    beta_m  = lambda Vm: 4*np.exp(-(Vm+65)/18)
    alpha_h = lambda Vm: 7e-2*np.exp(-(Vm+65)/20)
    beta_h  = lambda Vm: 1 / (1+np.exp(-(Vm+35)/10))

    # Initialization
    Vm = V_rest
    Time = np.arange(0,T+dt,dt)
    V = np.zeros_like(Time)

    n = alpha_n(Vm) / (alpha_n(Vm)+beta_n(Vm))
    m = alpha_m(Vm) / (alpha_m(Vm)+beta_m(Vm))
    h = alpha_h(Vm) / (alpha_h(Vm)+beta_h(Vm))
    for i,t in enumerate(Time):
        V[i] = Vm
        dVdt = (I-gl*(Vm-Vl)-gK*(n**4)*(Vm-VK)-gNa*(m**3)*h*(Vm-VNa))/Cm
        V_plus_dV = Vm+dt*dVdt
        n += (alpha_n(Vm)*(1-n) - beta_n(Vm)*n)*dt
        m += (alpha_m(Vm)*(1-m) - beta_m(Vm)*m)*dt
        h += (alpha_h(Vm)*(1-h) - beta_h(Vm)*h)*dt
        Vm = V_reset if Vm>V_thre else V_fire if V_plus_dV > V_thre else V_plus_dV
    return Time,V

def STDPsimulation(pre, post, dt=0.1, T=100):
    Time = np.arange(0,T+dt,dt)
    X  = np.zeros_like(Time); Y  = np.zeros_like(Time)
    Pre_Neuron  = np.zeros_like(Time); Post_Neuron = np.zeros_like(Time)
    Pre_Neuron[Time == pre] = 1; Post_Neuron[Time == post] = 1
    for i in range(len(Time)-1):
        X[i+1] = X[i]*np.exp(-0.1/20)+Pre_Neuron[i]
        Y[i+1] = Y[i]*np.exp(-0.1/20)+Post_Neuron[i]
    dW = -Pre_Neuron*Y + Post_Neuron*X
    res = [
        Time,
        (Pre_Neuron,Post_Neuron,X,Y,dW),
        {
            "label": ("Pre_Neuron", "Post_Neuron", "X", "Y", "dW"),
            "color": ("red", "blue", "red", "blue", "green")
        }
    ]
    return res
