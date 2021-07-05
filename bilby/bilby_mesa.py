#!/usr/bin/env python
import sys

import numpy as np
import bilby
from bilby.gw.likelihood import GravitationalWaveTransient
import memspectrum


class MESAGravitationalWaveTransient(GravitationalWaveTransient):

    def __init__(self, interferometers, waveform_generator,
                 priors=None, distance_marginalization=None,
                 phase_marginalization=None, time_marginalization=None):
        GravitationalWaveTransient.__init__(
            self, interferometers=interferometers, waveform_generator=waveform_generator, priors=priors,
            distance_marginalization=distance_marginalization,
            phase_marginalization=phase_marginalization,
            time_marginalization=time_marginalization)
        self._nll = np.nan
        self.m = memspectrum.MESA()

    def log_likelihood_ratio(self):
        return self.log_likelihood() - self.noise_log_likelihood()

    def log_likelihood(self):
        deltas = []
        for ifo in self.interferometers:
            psd_array = ifo.power_spectral_density._cache["psd_array"]
            frequency_mask = ifo.frequency_mask

            response = ifo.get_detector_response(
                self.waveform_generator.frequency_domain_strain(),
                self.parameters
            )
            signal_projected = bilby.core.utils.infft(
                response,
                self.waveform_generator.sampling_frequency
            )
            delta = ifo.time_domain_strain - signal_projected
            _ = self.m.solve(delta, m=5)
            _, psd = self.m.spectrum(
                dt=1 / ifo.sampling_frequency,
                onesided=True,
            )
            psd = np.concatenate([psd, [psd_array[-1]]])
            ifo.power_spectral_density._cache["psd_array"][frequency_mask] = psd[frequency_mask]

            deltas.append(-np.sum(np.log(2 * np.pi * ifo.power_spectral_density_array[
                ifo.frequency_mask]))
            )
            deltas.append(-2. / self.waveform_generator.duration * np.sum(
                abs(ifo.frequency_domain_strain[ifo.frequency_mask]) ** 2 /
                ifo.power_spectral_density_array[ifo.frequency_mask])
            )

        logl = GravitationalWaveTransient.log_likelihood_ratio(self) + np.sum(deltas)

        if np.isnan(logl):
            return -np.nan_to_num(np.inf)
        else:
            return logl

    def noise_log_likelihood(self):
        if np.isnan(self._nll):
            self._nll = 0
            for ifo in self.interferometers:
                self._nll += - sum(ifo.frequency_mask) / 2 - np.sum(
                    np.log(
                        2 * np.pi * abs(
                            ifo.frequency_domain_strain[ifo.frequency_mask]
                        )**2
                    )
                )
        return self._nll


duration = 4.
sampling_frequency = 512.
outdir = 'outdir'
np.random.seed(88170235)

injection_parameters = dict(
    chirp_mass=36., mass_ratio=0.9, a_1=0.0, a_2=0.0, tilt_1=0.0, tilt_2=0.0,
    phi_12=0.0, phi_jl=0.0, luminosity_distance=4000., theta_jn=0.4, psi=2.659,
    phase=1.3, geocent_time=1126259642.413, ra=1.375, dec=-1.2108)

waveform_arguments = dict(
    waveform_approximant='TaylorF2',
    reference_frequency=20.,
    minimum_frequency=20.,
    catch_waveform_errors=True,
)
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments)

ifos = bilby.gw.detector.InterferometerList(['H1'])
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency, duration=duration,
    start_time=injection_parameters['geocent_time'] - 3)
ifos.inject_signal(waveform_generator=waveform_generator,
                   parameters=injection_parameters)

priors = bilby.gw.prior.BBHPriorDict()
priors['geocent_time'] = bilby.core.prior.Uniform(
    minimum=injection_parameters['geocent_time'] - 0.1,
    maximum=injection_parameters['geocent_time'] + 0.1,
    name='geocent_time', latex_label='$t_c$', unit='$s$')
priors['chirp_mass'] = bilby.core.prior.Uniform(
    minimum=35, maximum=37, name="chirp_mass"
)

for key in bilby.gw.source.spin:
    if key in injection_parameters:
        priors[key] = injection_parameters[key]
for key in ['psi', 'ra', 'dec', 'geocent_time', 'phase', 'theta_jn', 'mass_ratio', 'luminosity_distance']:
    priors[key] = injection_parameters[key]

if "standard" in sys.argv:
    standard_likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        interferometers=ifos, waveform_generator=waveform_generator)

    result = bilby.run_sampler(
        likelihood=standard_likelihood, priors=priors, sampler='bilby_mcmc', nsamples=20000,
        check_point_delta_t=300, ntemps=1, thin_by_nact=0.2,
        use_ratio=False,
        injection_parameters=injection_parameters, outdir=outdir, label="standard")
else:
    mesa_likelihood = MESAGravitationalWaveTransient(
        interferometers=ifos, waveform_generator=waveform_generator)

    result = bilby.run_sampler(
        likelihood=mesa_likelihood, priors=priors, sampler='bilby_mcmc', nsamples=5000,
        check_point_delta_t=300, ntemps=1, thin_by_nact=0.2,
        use_ratio=False,
        injection_parameters=injection_parameters, outdir=outdir, label='mesa')

result.plot_corner()
