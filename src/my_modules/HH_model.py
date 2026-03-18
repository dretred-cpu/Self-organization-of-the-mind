import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import time
import pandas as pd
import os



class HodgkinHuxley_model():
    def __init__(self):
        # Parametry symulacji
        self.dt = 0.01
        self.T = 50
        self.time_HH = np.arange(0, self.T, self.dt)
        self.V = np.zeros(len(self.time_HH))
        self.m = np.zeros(len(self.time_HH))
        self.h = np.zeros(len(self.time_HH))
        self.n = np.zeros(len(self.time_HH))
        # Parametry modelu
        self.Cm = 1.0       # uF/cm^2
        self.gNa = 120.0    # mS/cm^2
        self.gK = 36.0
        self.gL = 0.3

        self.ENa = 50.0    # mV
        self.EK = -77.0
        self.EL = -54.387
        self.V[0] = -65.0
        self.V_initial = self.V[0]


    # Funkcje alpha i beta
    def alpha_m(self, V):
        return 0.1 * (self.V + 40) / (1 - np.exp(-(self.V + 40) / 10))

    def beta_m(self, V):
        return 4.0 * np.exp(-(self.V + 65) / 18)
    def alpha_h(self, V):
        return 0.07 * np.exp(-(self.V + 65) / 20)

    def beta_h(self, V):
        return 1 / (1 + np.exp(-(self.V + 35) / 10))

    def alpha_n(self, V):
        return 0.01 * (self.V + 55) / (1 - np.exp(-(self.V + 55) / 10))

    def beta_n(self, V):
        return 0.125 * np.exp(-(self.V + 65) / 80)
    
    def initial_conditions(self, V):
    # Warunki początkowe (stan spoczynkowy)
        self.m[0] = self.alpha_m(V[0]) / (self.alpha_m(V[0]) + self.beta_m(V[0]))
        self.h[0] = self.alpha_h(V[0]) / (self.alpha_h(V[0]) + self.beta_h(V[0]))
        self.n[0] = self.alpha_n(V[0]) / (self.alpha_n(V[0]) + self.beta_n(V[0]))

    def I_ext(self, t):
        return 10.0 if 10 < t < 40 else 0.0


    # Główna pętla
    def run_simulation(self):
        Cm = self.Cm     # uF/cm^2
        gNa = self.gNa   # mS/cm^2
        gK = self.gK
        gL = self.gL

        ENa = self.ENa   # mV
        EK = self.EK
        EL =self.EL
        long = len(self.time_HH)


        self.V[i - 1] = self.V_initial
        
        for i in range(1, len(self.time_HH)):
            if i == 1:
                Vm = self.V_initial
            else:
                Vm = self.V[i-1]

            INa = self.gNa * (self.m[i-1]**3) * self.h[i-1] * (Vm - self.ENa)
            IK  = self.gK  * (self.n[i-1]**4) * (Vm - self.EK)
            IL  = self.gL * (Vm - self.EL)
            self.V[i] = Vm + self.dt * (self.I_ext(self.time_HH[i]) - INa - IK - IL) / self.Cm
            self.m[i] = self.m[i-1] + self.dt * (self.alpha_m(Vm)*(1-self.m[i-1]) - self.beta_m(Vm)*self.m[i-1])
            self.h[i] = self.h[i-1] + self.dt * (self.alpha_h(Vm)*(1-self.h[i-1]) - self.beta_h(Vm)*self.h[i-1])
            self.n[i] = self.n[i-1] + self.dt * (self.alpha_n(Vm)*(1-self.n[i-1]) - self.beta_n(Vm)*self.n[i-1])
    
        df_V = pd.DataFrame([self.V], columns=self.time_HH)
        print(df_V)
        print(df_V.shape)
        df_V.index = ["V"]
        df_V.to_csv("V_output.csv", index=True)