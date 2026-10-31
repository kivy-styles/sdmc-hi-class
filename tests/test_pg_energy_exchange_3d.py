"""Independent checks for 3D periodic Jordan geometric energy exchange."""
import unittest
from math import exp
import numpy as np
from sdmc.pg_energy_exchange import Epoch
from sdmc.pg_energy_exchange_3d import (
    averaged_energy_pressure_3d, trace_eom_3d,
    mean_jordan_exchange_3d, demo)
from sdmc.pg_linear_perturbations import trace_dynamic_derivatives

class PGThreeDimensionalExchange(unittest.TestCase):
    def setUp(self):
        self.n=24; self.L=2*np.pi
        ax=np.arange(self.n)*self.L/self.n
        self.x,self.y,self.z=np.meshgrid(ax,ax,ax,indexing="ij")
        self.ell=.2*np.cos(2*self.x)+.13*np.sin(3*self.y)
        self.vel=.07*np.sin(2*self.x)-.03*np.cos(4*self.z)
        self.u=.08*np.cos(self.x)+.04*np.sin(2*self.z)
        self.udot=.019*np.sin(self.x)-.011*np.cos(2*self.y)
        self.e=Epoch(.91,1.31,.043,.0645,1.3,.6,.21,.009,
                 Kdot=.017,c2dot=.013,mass2dot=.007,wetadot=-.0004)

    def test_nonzero_inhomogeneous_energy_with_zero_mean_trace(self):
        t=averaged_energy_pressure_3d(self.ell,self.vel,self.u,self.e,self.L)
        self.assertAlmostEqual(float(self.ell.mean()),0.,places=12)
        self.assertGreater(t["E"],0.)
        self.assertGreater(t["curvature"],0.)

    def test_cubic_3d_eom_reduces_to_fourier_trace_equation(self):
        e=self.e;k=3
        A=.4;B=-.1;U=.2
        ell=A*np.cos(k*self.x)
        vel=B*np.cos(k*self.x)
        u=U*np.cos(k*self.x)
        acc=trace_eom_3d(ell,vel,u,e,self.L)
        ref=trace_dynamic_derivatives(
            A,B,e.H_j,e.Kdot/e.K,-U,k,e.a_j,e.a_s,
            3*e.K,e.c2,e.mass2,e.weta)[1]
        self.assertAlmostEqual(float(acc[0,0,0]),ref,
                               delta=abs(ref)*1e-10+1e-10)

    def test_spatial_3d_energy_balance_matches_finite_difference(self):
        e=self.e
        acc=trace_eom_3d(self.ell,self.vel,self.u,e,self.L)
        present=averaged_energy_pressure_3d(
            self.ell,self.vel,self.u,e,self.L)
        J=mean_jordan_exchange_3d(
            self.ell,self.vel,self.u,self.udot,e,self.L)
        h=2e-5
        def epoch(t):
            return Epoch(
               a_j=e.a_j*exp(e.H_j*t),a_s=e.a_s*exp(e.H_s*t),
               H_j=e.H_j,H_s=e.H_s,
               K=e.K+e.Kdot*t,c2=e.c2+e.c2dot*t,
               mass2=e.mass2+e.mass2dot*t,
               weta=e.weta+e.wetadot*t,
               Kdot=e.Kdot,c2dot=e.c2dot,
               mass2dot=e.mass2dot,wetadot=e.wetadot)
        def E(t):
            return averaged_energy_pressure_3d(
                self.ell+self.vel*t+.5*acc*t*t,
                self.vel+acc*t,self.u+self.udot*t,
                epoch(t),self.L)["E"]
        derivative=(E(h)-E(-h))/(2*h)
        expected=J-3*e.H_j*(present["E"]+present["P_restricted"])
        self.assertAlmostEqual(derivative,expected,
                               delta=max(1e-8,abs(expected)*2e-6))

    def test_zero_curvature_constant_coeff_gives_zero_exchange(self):
        e=Epoch(1.,1.,.05,.075,1.2,.9,.3,0.)
        u=np.zeros_like(self.u)
        z=np.zeros_like(self.u)
        J=mean_jordan_exchange_3d(self.ell,self.vel,u,z,e,self.L)
        self.assertAlmostEqual(J,0.,places=12)

    def test_clock_curvature_exchange_is_negative_for_unforced_field(self):
        e=Epoch(1.,1.,.05,.075,1.2,.9,.3,.01)
        z=np.zeros_like(self.u)
        t=averaged_energy_pressure_3d(self.ell,z,z,e,self.L)
        J=mean_jordan_exchange_3d(self.ell,z,z,z,e,self.L)
        self.assertAlmostEqual(J,-4*e.H_s*t["curvature"],places=11)

    def test_driver_work_could_counteract_clock_exchange(self):
        e=Epoch(1.,1.,.05,.075,1.2,.9,.3,.01)
        z=np.zeros_like(self.u)
        without=mean_jordan_exchange_3d(self.ell,z,self.u,z,e,self.L)
        withdriver=mean_jordan_exchange_3d(self.ell,z,self.u,
                                         2*self.u,e,self.L)
        self.assertGreater(withdriver,without)

    def test_3d_demo_is_not_cosmological_dust_claim(self):
        d=demo()
        self.assertGreater(d["mean_energy"],0.)
        self.assertLess(d["mean_pressure_restricted"],0.)
        self.assertLess(d["exchange_residual"],0.)
        self.assertAlmostEqual(d["p_v"],1.5)
        self.assertIn("pressureless attractor",d["missing"])

if __name__=="__main__":
    unittest.main()
